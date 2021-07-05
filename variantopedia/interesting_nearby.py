import operator
import re
from collections import defaultdict
from functools import reduce
from typing import Dict

from django.conf import settings
from django.db.models import Count, Sum, Q

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from classification.enums import ClinicalSignificance
from classification.models import Classification
from snpdb.models import Variant


def get_method_summaries(variant, distance=None):
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

    start = variant.start - distance
    end = variant.end + distance
    range_summary = f"{variant.locus.contig.name}: {start}-{end} ({distance} bases)"

    return {
        "codon": codon_summary,
        "exon": exon_summary,
        "domain": domain_summary,
        "range": range_summary,
    }


def get_nearby_qs(variant, annotation_version, distance=None):
    if distance is None:
        distance = settings.VARIANT_DETAILS_NEARBY_RANGE
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    q = Variant.get_no_reference_q() & ~Q(pk=variant.pk)  # Exclude ref and self
    qs = qs.filter(q)

    return {
        "codon": filter_variant_codon(qs, variant),
        "exon": filter_variant_exon(qs, variant),
        "domain": filter_variant_domain(qs, variant),
        "range": filter_variant_range(qs, variant, distance=distance),
    }


def get_nearby_summaries(user, variant, annotation_version, distance=None, clinical_significance=False):
    if distance is None:
        distance = settings.VARIANT_DETAILS_NEARBY_RANGE
    nearby_qs = get_nearby_qs(variant, annotation_version, distance=distance)
    kwargs = {
        "user": user,
        "genome_build": annotation_version.genome_build,
        "clinical_significance": clinical_significance
    }
    return {f"{region}_summary": interesting_summary(qs, **kwargs) for region, qs in nearby_qs.items()}


def interesting_summary(qs, user, genome_build, total=True, clinvar=True, clinical_significance=False):
    counts = interesting_counts(qs, user, genome_build, clinical_significance=clinical_significance)
    # print(counts)
    summary = None
    if num_variants := counts['total']:
        classification_types = {
            "Classifications": "classification",
        }
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
    return summary


def interesting_counts(qs, user, genome_build, clinical_significance=False):
    """ qs: Variant queryset annotated w/VariantZygosityCountCollection """

    agg_kwargs = {
        "total": Count("id"),
        "REF": Sum("global_variant_zygosity__ref_count"),
        "HET": Sum("global_variant_zygosity__het_count"),
        "HOM_ALT": Sum("global_variant_zygosity__hom_count"),
    }

    clinical_significance_list = [c[0] for c in ClinicalSignificance.SHORT_CHOICES]
    q_classification = Classification.get_variant_q(user, genome_build,
                                                    clinical_significance_list=clinical_significance_list)
    classifications = {
        "clinvar": (Q(clinvar__isnull=False), "highest_pathogenicity"),
        "classification": (q_classification, "clinical_significance")
    }

    for classification, (classification_q, clinical_significance_path) in classifications.items():
        agg_kwargs[f"{classification}_count"] = Count("id", filter=classification_q)
        if clinical_significance:
            for cs in clinical_significance_list:
                q_clinical_significance = Q(**{f"{classification}__{clinical_significance_path}": cs})
                if classification_q:
                    q_clinical_significance &= classification_q
                agg_kwargs[f"{classification}_{cs}"] = Count("id", filter=q_clinical_significance)

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
        if m := re.match(r".*(:c\.\d+)", hgvs_c):  # Pulls out eg ":c.1057"
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
