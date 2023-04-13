import operator
import re
from collections import defaultdict, Counter
from functools import reduce
from typing import Dict, Tuple

from django.conf import settings
from django.contrib.postgres.aggregates import StringAgg
from django.db.models import Count, Sum, Q

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models import VariantTranscriptAnnotation, AnnotationVersion
from classification.enums import ClinicalSignificance
from classification.models import Classification
from genes.models import GeneSymbol
from snpdb.models import Variant, VariantZygosityCountCollection


def get_method_summaries(variant, annotation_version, distance=None):
    if distance is None:
        distance = settings.VARIANT_DETAILS_NEARBY_RANGE

    if transcripts_and_codons := get_transcripts_and_codons(variant):
        codon_summary = " or ".join([f"Transcript: {t}, Codon: {codon}" for (t, codon) in transcripts_and_codons.items()])
    else:
        codon_summary = "Intergenic - no codon search performed"

    if transcripts_and_exons := get_transcript_and_exons(variant):
        exon_summary = " or ".join([f"Transcript: {t}, Exon: {e}" for t, e in transcripts_and_exons.items()])
    else:
        exon_summary = "Intergenic - no exon search performed."

    if transcript_and_domains := get_transcript_and_domains(variant):
        domains = []
        for t, domain_set in transcript_and_domains.items():
            domains.append(f"Transcript: {t} and domain contains: {', '.join(domain_set)}")
        domain_summary = " or ".join(domains)
    else:
        domain_summary = "Not in an annotated domain - no domain search performed"

    vav = annotation_version.variant_annotation_version
    if GeneSymbol.overlapping_variant(variant, vav).exists():
        gene_summary = f"{settings.ANNOTATION_VEP_DISTANCE}bp up or downstream of a gene"
    else:
        gene_summary = "Intergenic, no gene search performed"

    start = variant.start - distance
    end = variant.end + distance
    range_summary = f"{variant.locus.contig.name}: {start}-{end} ({distance} bases)"

    return {
        "codon": codon_summary,
        "exon": exon_summary,
        "domain": domain_summary,
        "genes": gene_summary,
        "range": range_summary,
    }


def get_nearby_qs(variant, annotation_version, distance=None):
    if distance is None:
        distance = settings.VARIANT_DETAILS_NEARBY_RANGE
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
    q = Variant.get_no_reference_q() & ~Q(pk=variant.pk)  # Exclude ref and self
    qs = qs.filter(q)

    qs_dict = {
        "codon": filter_variant_codon(qs, variant),
        "exon": filter_variant_exon(qs, variant),
        "domain": filter_variant_domain(qs, variant),
        "range": filter_variant_range(qs, variant, distance=distance),
    }
    if settings.VARIANT_DETAILS_NEARBY_SHOW_GENE:
        qs_dict["genes"] = get_gene_symbol_filters(qs, variant, annotation_version)
    return qs_dict


def get_nearby_summaries(user, variant, annotation_version, distance=None, clinical_significance=False):
    if distance is None:
        distance = settings.VARIANT_DETAILS_NEARBY_RANGE
    nearby_qs_dict = get_nearby_qs(variant, annotation_version, distance=distance)
    kwargs = {
        "user": user,
        "genome_build": annotation_version.genome_build,
        "clinical_significance": clinical_significance
    }

    def get_summary(qs, prefix=None):
        summary, tag_counts = interesting_summary(qs, **kwargs)
        if prefix:
            prefix_label = [prefix]
        else:
            prefix_label = []
        summary_label = "_".join(prefix_label + ["summary"])
        tag_counts_label = "_".join(prefix_label + ["tag_counts"])
        return {
            summary_label: summary,
            tag_counts_label: tag_counts,
        }

    gene_summaries = {}
    for gene_symbol, qs in nearby_qs_dict.pop("genes", {}).items():
        gene_summaries[gene_symbol] = get_summary(qs)

    nearby_summaries = {}
    for region, qs in nearby_qs_dict.items():
        nearby_summaries.update(get_summary(qs, prefix=region))

    nearby_summaries["genes"] = gene_summaries
    return nearby_summaries


def interesting_summary(qs, user, genome_build, total=True, clinvar=True, classifications=True,
                        clinical_significance=False) -> Tuple[str, Dict[str, int]]:
    """ returns a strin summary, and dict of tag counts (so you can format w/colors) """
    counts = interesting_counts(qs, user, genome_build, clinical_significance=clinical_significance)
    summary = ""
    tag_counts = {}

    if num_variants := counts['total']:
        classification_types = {}
        if classifications:
            classification_types["Classifications"] = "classification"

        if clinvar:
            classification_types["ClinVar"] = "clinvar"

        summaries = []
        for label, classification in classification_types.items():
            if classification_count := counts[f"{classification}_count"]:
                if clinical_significance:
                    cs_summaries = []
                    for cs, cs_label in ClinicalSignificance.SHORT_CHOICES:
                        field = f"{classification}_{cs}"
                        if c := counts[field]:
                            cs_summaries.append(f"{cs_label}: {c}")
                    classification_summary = ", ".join(cs_summaries)
                    classification_summary = f"({classification_summary})"
                else:
                    classification_summary = classification_count
                summaries.append(f"{label}: {classification_summary}")

        if tags := counts.get("tags"):
            counter = Counter(tags.split("|"))
            tag_counts = dict(counter.most_common())  # Sort highest -> lowest

        zygosity_counts = []
        for zygosity in ["REF", "HET", "HOM_ALT"]:
            if c := counts[zygosity]:
                zygosity_counts.append(f"{zygosity}: {c}")

        if zygosity_counts:
            db_summary = f"DB: ({', '.join(zygosity_counts)})"
            summaries.append(db_summary)

        optional_summary = ", ".join(summaries)

        if total:
            summary = f"{num_variants} variants. "
        else:
            summary = ""
        if optional_summary:
            summary += optional_summary
    return summary, tag_counts


def interesting_counts(qs, user, genome_build, clinical_significance=False):
    """ qs: Variant queryset annotated w/VariantZygosityCountCollection """

    agg_kwargs = {
        "total": Count("id", distinct=True),
        "REF": Sum("global_variant_zygosity__ref_count"),
        "HET": Sum("global_variant_zygosity__het_count"),
        "HOM_ALT": Sum("global_variant_zygosity__hom_count"),
        "tags": StringAgg("variantallele__allele__varianttag__tag", delimiter='|'),
    }

    clinical_significance_list = [c[0] for c in ClinicalSignificance.SHORT_CHOICES]
    classification_qs = Classification.get_classifications_qs(user,
                                                              clinical_significance_list=clinical_significance_list)
    classifications_ids = list(classification_qs.values_list("pk", flat=True))
    q_classification = Q(classification__in=classifications_ids)
    q_classification &= Classification.get_variant_q_from_classification_qs(classification_qs, genome_build)

    classifications = {
        "clinvar": ("clinvar", Q(clinvar__isnull=False), "highest_pathogenicity"),
        "classification": ("classification", q_classification, "clinical_significance")
    }

    for classification, (count_path, classification_q, clinical_significance_path) in classifications.items():
        agg_kwargs[f"{classification}_count"] = Count(count_path, filter=classification_q, distinct=True)
        if clinical_significance:
            for cs in clinical_significance_list:
                q_clinical_significance = Q(**{f"{classification}__{clinical_significance_path}": cs})
                if classification_q:
                    q_clinical_significance &= classification_q
                agg_kwargs[f"{classification}_{cs}"] = Count(count_path, filter=q_clinical_significance, distinct=True)

    return qs.aggregate(**agg_kwargs)


def filter_variant_range(qs, variant: Variant, distance):
    start = variant.start - distance
    end = variant.end + distance
    annotation_kwargs, q = Variant.get_overlap_annotate_and_q(variant.locus.contig, start, end)
    return qs.annotate(**annotation_kwargs).filter(q)


def get_transcripts_and_codons(variant: Variant):
    transcript_codons = {}
    transcript_qs = variant.varianttranscriptannotation_set.filter(exon__isnull=False, hgvs_c__isnull=False)
    for t, hgvs_c in transcript_qs.values_list("transcript_id", "hgvs_c"):
        if m := re.match(r".*(:c\.\d+)", hgvs_c):  # Pulls out e.g. ":c.1057"
            codon = m.group(1)
            transcript_codons[t] = codon
    return transcript_codons


def filter_variant_codon(qs, variant: Variant):
    q_or = []
    for transcript_id, codon in get_transcripts_and_codons(variant).items():
        codon_regex = codon + r"[^\d]"
        q_or.append(Q(varianttranscriptannotation__transcript_id=transcript_id,
                      varianttranscriptannotation__hgvs_c__regex=codon_regex))
    if q_or:
        q = reduce(operator.or_, q_or)
        qs = qs.filter(q).distinct()
    else:
        qs = qs.none()
    return qs


def get_transcript_and_exons(variant) -> Dict:
    transcript_qs = variant.varianttranscriptannotation_set.filter(exon__isnull=False)
    return dict(transcript_qs.values_list("transcript_id", "exon"))


def filter_variant_exon(qs, variant: Variant):
    q_or = []
    for t, e in get_transcript_and_exons(variant).items():
        q_or.append(Q(varianttranscriptannotation__transcript_id=t, varianttranscriptannotation__exon=e))
    if q_or:
        q = reduce(operator.or_, q_or)
        qs = qs.filter(q).distinct()
    else:
        qs = qs.none()
    return qs


def get_transcript_and_domains(variant) -> Dict[str, set]:
    transcript_and_domain = defaultdict(set)
    transcript_qs = variant.varianttranscriptannotation_set.filter(interpro_domain__isnull=False)
    for t, interpro_domain in transcript_qs.values_list("transcript_id", "interpro_domain"):
        transcript_and_domain[t].update(interpro_domain.split("&"))  # VEP separator
    return transcript_and_domain


def filter_variant_domain(qs, variant: Variant):
    q_or = []

    for t, domain_set in get_transcript_and_domains(variant).items():
        for domain in domain_set:
            q_or.append(Q(varianttranscriptannotation__transcript_id=t,
                          varianttranscriptannotation__interpro_domain__contains=domain))
    if q_or:
        q = reduce(operator.or_, q_or)
        qs = qs.filter(q).distinct()
    else:
        qs = qs.none()
    return qs


def get_gene_symbol_filters(qs, variant: Variant, annotation_version: AnnotationVersion) -> Dict:
    gene_qs_dict = {}
    vav = annotation_version.variant_annotation_version
    for gene_symbol in GeneSymbol.overlapping_variant(variant, vav).order_by("pk"):
        genes_qs = vav.gene_annotation_release.genes_for_symbol(gene_symbol)
        q = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, genes_qs)
        gene_qs_dict[str(gene_symbol)] = qs.filter(q)
    return gene_qs_dict
