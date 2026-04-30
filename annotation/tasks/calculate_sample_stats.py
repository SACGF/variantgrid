import logging
import time
from collections import defaultdict
from typing import Optional

import celery
import numpy as np
from django.conf import settings
from django.db.models.query_utils import Q

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models import AnnotationVersion, CohortGenotypeVariantAnnotationStats, \
    CohortGenotypeGeneAnnotationStats, CohortGenotypeClinVarAnnotationStats
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import InvalidAnnotationVersionError, VCFAnnotationStats
from eventlog.models import create_event
from library.django_utils import thread_safe_unique_together_get_or_create
from library.enums.log_level import LogLevel
from library.genomics.vcf_enums import VariantClass, VCFSymbolicAllele
from library.git import Git
from library.log_utils import get_traceback
from library.utils.json_utils import canonical_filter_key
from snpdb.models import Cohort, CohortGenotypeStats, ImportStatus, VCF, Variant, \
    SampleStatsCodeVersion, Sequence, VCFLengthStatsCollection, VCFLengthStats
from snpdb.models import Zygosity
from snpdb.models.models_genome import GenomeBuild


@celery.shared_task
def calculate_vcf_stats(vcf_id, annotation_version_id):
    vcf = VCF.objects.get(pk=vcf_id)
    try:
        annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
        msg = f"Build for annotation: {annotation_version.genome_build} and VCF: {vcf.genome_build} match"
        assert annotation_version.genome_build == vcf.genome_build, msg
        calculate_cohort_stats(vcf.cohort, annotation_version)
        _compute_vcf_specific_stats(vcf, annotation_version)
        vcf.sample_set.update(import_status=ImportStatus.SUCCESS)
    except:
        tb = get_traceback()
        create_event(None, 'calculate_vcf_stats', details=tb, severity=LogLevel.ERROR)
        vcf.sample_set.update(import_status=ImportStatus.ERROR)
        raise


def _get_sample_stats_code_version() -> SampleStatsCodeVersion:
    """
        Version | Notes on sample stats code changes
        0 - Set for historical versions created before code_version added
        1 - Existing code when code_version added
        2 - Added VCFLengthStats
        3 - Optimisations
    """
    code_version, _ = thread_safe_unique_together_get_or_create(SampleStatsCodeVersion,
                                                                name="SampleStats",
                                                                version=3,
                                                                code_git_hash=Git(settings.BASE_DIR).hash)
    return code_version


def _compute_vcf_specific_stats(vcf: VCF, annotation_version: AnnotationVersion):
    """ VCF-keyed stats — VEP-skipped count and variant-class length histograms.
        Per-sample / cohort-aggregate stats are written by calculate_cohort_stats. """
    start = time.time()
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    qs = qs.filter(Variant.get_no_reference_q())
    qs = vcf.get_variant_qs(qs)
    columns = [
        "locus__ref__seq",
        "alt__seq",
        "svlen",
        "variantannotation__vep_skipped_reason",
        "variantannotation__hgvs_g",
    ]
    vep_skipped_count = 0
    variant_class_lengths = defaultdict(list)

    for row in qs.values_list(*columns).iterator(chunk_size=10_000):
        ref, alt, svlen, vep_skipped, hgvs_g = row
        if vep_skipped:
            vep_skipped_count += 1
        ref_len = len(ref)
        alt_len = len(alt)
        if ref_len == 1 and alt_len == 1:
            continue  # SNP — no length histogram
        variant_class = _get_variant_class(alt, hgvs_g)
        variant_length = abs(svlen) if svlen else max(alt_len, ref_len)
        variant_class_lengths[variant_class].append(variant_length)

    if vep_skipped_count:
        stats, created = VCFAnnotationStats.objects.get_or_create(
            vcf=vcf, variant_annotation_version=annotation_version.variant_annotation_version,
            defaults={"vep_skipped_count": vep_skipped_count})
        if not created and vep_skipped_count > stats.vep_skipped_count:
            stats.vep_skipped_count = vep_skipped_count
            stats.save()

    if variant_class_lengths:
        code_version = _get_sample_stats_code_version()
        collection, created = VCFLengthStatsCollection.objects.update_or_create(
            vcf=vcf, defaults={"code_version": code_version})
        if not created:
            collection.vcflengthstats_set.all().delete()

        vcf_length_stats = []
        for variant_class, lengths in variant_class_lengths.items():
            lengths = np.array(lengths)
            is_log = False
            abs_lengths = np.abs(lengths)
            max_min_ratio = abs_lengths.max() / (abs_lengths.min() + 1)
            logging.info(f"Variant class {variant_class} has {abs_lengths.max()=} / {abs_lengths.min()=} = {max_min_ratio=}")
            if max_min_ratio > 1000:  # Spans 3 orders of magnitude
                is_log = True
                lengths = np.log10(lengths)

            histogram_counts, histogram_bin_edges = np.histogram(lengths, bins=50)
            vcf_length_stats.append(VCFLengthStats(
                collection=collection, variant_class=variant_class, is_log=is_log,
                histogram_counts=histogram_counts.tolist(),
                histogram_bin_edges=histogram_bin_edges.tolist()))
        if vcf_length_stats:
            VCFLengthStats.objects.bulk_create(vcf_length_stats)

    end = time.time()
    logging.info("VCF-specific stats for %s took %.2f seconds", vcf, end - start)


@celery.shared_task
def calculate_cohort_stats_task(cohort_id, annotation_version_id):
    cohort = Cohort.objects.get(pk=cohort_id)
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
    calculate_cohort_stats(cohort, annotation_version)


# Field indexes for the per-bucket counter arrays. Keep in lock-step with
# _COUNTER_FIELDS below — index N in the array is the field named at position N.
_COUNTER_FIELDS = (
    "variant_count",
    "snp_count",
    "insertions_count",
    "deletions_count",
    "ref_count",
    "het_count",
    "hom_count",
    "unk_count",
    "x_hom_count",
    "x_het_count",
    "x_unk_count",
    "variant_dbsnp_count",
    "snp_dbsnp_count",
    "insertions_dbsnp_count",
    "deletions_dbsnp_count",
    "ref_high_or_moderate_count",
    "het_high_or_moderate_count",
    "hom_high_or_moderate_count",
    "unk_high_or_moderate_count",
    "gene_count",
    "ref_omim_phenotype_count",
    "het_omim_phenotype_count",
    "hom_omim_phenotype_count",
    "unk_omim_phenotype_count",
    "clinvar_count",
    "ref_clinvar_pathogenic_count",
    "het_clinvar_pathogenic_count",
    "hom_clinvar_pathogenic_count",
    "unk_clinvar_pathogenic_count",
)
_FIELD_INDEX = {name: i for i, name in enumerate(_COUNTER_FIELDS)}
_NUM_FIELDS = len(_COUNTER_FIELDS)


def _new_counter():
    return [0] * _NUM_FIELDS


def _aggregate_zygosity_priority(zygosities: str) -> str:
    """ For aggregate rows, classify a variant by the highest-priority zygosity
        present across the cohort's samples — hom > het > ref > unk. This makes
        ref+het+hom+unk_count sum to variant_count for the aggregate row, the
        invariant SampleStats relies on. """
    if Zygosity.HOM_ALT in zygosities:
        return Zygosity.HOM_ALT
    if Zygosity.HET in zygosities:
        return Zygosity.HET
    if Zygosity.HOM_REF in zygosities:
        return Zygosity.HOM_REF
    return Zygosity.UNKNOWN_ZYGOSITY


def _trio_predicates(mum_zyg, dad_zyg, proband_zyg, contig, mother_affected, father_affected):
    """ Returns the set of inheritance modes the variant qualifies for under
        the trio's mom/dad/proband zygosities. Mirrors the predicates in
        analysis/models/nodes/sources/trio_node.py — keep in sync. """
    NO_VARIANT = (Zygosity.MISSING, Zygosity.HOM_REF)
    HAS_VARIANT = (Zygosity.HET, Zygosity.HOM_ALT)
    modes = set()
    # Denovo: mum & dad no variant; proband has variant
    if mum_zyg in NO_VARIANT and dad_zyg in NO_VARIANT and proband_zyg in HAS_VARIANT:
        modes.add("denovo")
    # Autosomal recessive: mum HET, dad HET, proband HOM_ALT
    if mum_zyg == Zygosity.HET and dad_zyg == Zygosity.HET and proband_zyg == Zygosity.HOM_ALT:
        modes.add("autosomal_recessive")
    # Autosomal dominant: each parent's expected zyg depends on affected status;
    # proband always has variant.
    mum_dominant = (mum_zyg in HAS_VARIANT) if mother_affected else (mum_zyg in NO_VARIANT)
    dad_dominant = (dad_zyg in HAS_VARIANT) if father_affected else (dad_zyg in NO_VARIANT)
    if mum_dominant and dad_dominant and proband_zyg in HAS_VARIANT:
        modes.add("autosomal_dominant")
    # X-linked recessive: mum HET, dad unconstrained, proband HOM_ALT, on chrX
    if contig == "X" and mum_zyg == Zygosity.HET and proband_zyg == Zygosity.HOM_ALT:
        modes.add("x_linked_recessive")
    return modes


def calculate_cohort_stats(cohort: Cohort, annotation_version: AnnotationVersion):
    """ Iterate the cohort's variant queryset once and populate the four
        CohortGenotype*Stats classes for:
          - per-sample rows (sample IS NOT NULL, filter_key NULL) — every cohort member
          - aggregate row (sample IS NULL, filter_key NULL)
          - filter-keyed aggregate rows (sample IS NULL, filter_key=...) for any keys
            returned by get_filter_keys_to_precompute_for_cohort(cohort)
        All rows are FK'd to cohort.cohort_genotype_collection (the current CGC).
        Idempotent: any existing rows for this CGC + the four annotation versions
        are deleted and rewritten. """
    from analysis.models.nodes.sources._stats_cache import get_filter_keys_to_precompute_for_cohort

    cgc = cohort.cohort_genotype_collection
    code_version = _get_sample_stats_code_version()
    cohort_samples = list(cohort.get_cohort_samples())
    sample_index_pairs = [(cs.sample, cs.cohort_genotype_packed_field_index) for cs in cohort_samples]
    has_filters = bool(_cohort_has_filters(cohort))

    # Filter-keyed buckets to precompute (always includes None for raw aggregate).
    precompute_keys = list(get_filter_keys_to_precompute_for_cohort(cohort))
    if None not in precompute_keys:
        precompute_keys = [None] + precompute_keys

    # Trio context, if present (used to evaluate inheritance predicates per variant).
    trio = cohort.trio_set.first()
    trio_indexes = None
    if trio is not None:
        cs_for = {cs.pk: cs for cs in cohort_samples}
        try:
            trio_indexes = (
                cs_for[trio.mother_id].cohort_genotype_packed_field_index,
                cs_for[trio.father_id].cohort_genotype_packed_field_index,
                cs_for[trio.proband_id].cohort_genotype_packed_field_index,
            )
        except KeyError:
            trio_indexes = None  # Trio not coherent with this cohort; skip inheritance bucketing

    # Accumulators — one counter array per bucket key.
    # Per-sample: (sample_pk, passing_filter_bool) → counter
    sample_counters: dict[tuple[int, bool], list[int]] = defaultdict(_new_counter)
    # Aggregate / filter-keyed: (filter_key_str_or_None, passing_filter_bool) → counter
    aggregate_counters: dict[tuple[Optional[str], bool], list[int]] = defaultdict(_new_counter)
    # Compound-het per-gene accumulator: (gene_id, passing_filter_bool) → {"mum", "dad"} sources seen
    chet_gene_sources: dict[tuple[int, bool], set] = defaultdict(set)
    # Per-variant chet contributions queued for post-processing once we know the genes:
    #   list of (gene_id, passing_filter_bool, proband_zyg, contig, ref_len, alt_len, alt, hgvs_g,
    #            dbsnp, impact, omim, clinvar_pathogenic, transcript)
    # We accumulate contributions per gene; only genes with both "mum" and "dad" sources
    # contribute to the compound_het aggregate counters.
    chet_pending: list[tuple] = []

    # Local bindings for hot loop
    HET, HOM_ALT, HOM_REF, MISSING, UNK = (
        Zygosity.HET, Zygosity.HOM_ALT, Zygosity.HOM_REF, Zygosity.MISSING, Zygosity.UNKNOWN_ZYGOSITY
    )
    HAS_VARIANT = (HET, HOM_ALT)
    NO_VARIANT = (MISSING, HOM_REF)
    pathogenicity_high_mod = {PathogenicityImpact.HIGH, PathogenicityImpact.MODERATE}
    symbolic_insertions = {VCFSymbolicAllele.DUP, VCFSymbolicAllele.INS}

    # Field indexes
    F_VARIANT = _FIELD_INDEX["variant_count"]
    F_SNP = _FIELD_INDEX["snp_count"]
    F_INS = _FIELD_INDEX["insertions_count"]
    F_DEL = _FIELD_INDEX["deletions_count"]
    F_REF = _FIELD_INDEX["ref_count"]
    F_HET = _FIELD_INDEX["het_count"]
    F_HOM = _FIELD_INDEX["hom_count"]
    F_UNK = _FIELD_INDEX["unk_count"]
    F_X_HOM = _FIELD_INDEX["x_hom_count"]
    F_X_HET = _FIELD_INDEX["x_het_count"]
    F_X_UNK = _FIELD_INDEX["x_unk_count"]
    F_VAR_DBSNP = _FIELD_INDEX["variant_dbsnp_count"]
    F_SNP_DBSNP = _FIELD_INDEX["snp_dbsnp_count"]
    F_INS_DBSNP = _FIELD_INDEX["insertions_dbsnp_count"]
    F_DEL_DBSNP = _FIELD_INDEX["deletions_dbsnp_count"]
    F_REF_HM = _FIELD_INDEX["ref_high_or_moderate_count"]
    F_HET_HM = _FIELD_INDEX["het_high_or_moderate_count"]
    F_HOM_HM = _FIELD_INDEX["hom_high_or_moderate_count"]
    F_UNK_HM = _FIELD_INDEX["unk_high_or_moderate_count"]
    F_GENE = _FIELD_INDEX["gene_count"]
    F_REF_OMIM = _FIELD_INDEX["ref_omim_phenotype_count"]
    F_HET_OMIM = _FIELD_INDEX["het_omim_phenotype_count"]
    F_HOM_OMIM = _FIELD_INDEX["hom_omim_phenotype_count"]
    F_UNK_OMIM = _FIELD_INDEX["unk_omim_phenotype_count"]
    F_CLINVAR = _FIELD_INDEX["clinvar_count"]
    F_REF_CV = _FIELD_INDEX["ref_clinvar_pathogenic_count"]
    F_HET_CV = _FIELD_INDEX["het_clinvar_pathogenic_count"]
    F_HOM_CV = _FIELD_INDEX["hom_clinvar_pathogenic_count"]
    F_UNK_CV = _FIELD_INDEX["unk_clinvar_pathogenic_count"]

    PRIORITY_TO_F = {
        HOM_ALT: (F_HOM, F_HOM_HM, F_HOM_OMIM, F_HOM_CV),
        HET: (F_HET, F_HET_HM, F_HET_OMIM, F_HET_CV),
        HOM_REF: (F_REF, F_REF_HM, F_REF_OMIM, F_REF_CV),
        UNK: (F_UNK, F_UNK_HM, F_UNK_OMIM, F_UNK_CV),
        MISSING: (F_UNK, F_UNK_HM, F_UNK_OMIM, F_UNK_CV),
    }

    PER_SAMPLE_TO_F = {
        HOM_ALT: (F_HOM, F_HOM_HM, F_HOM_OMIM, F_HOM_CV, F_X_HOM),
        HET: (F_HET, F_HET_HM, F_HET_OMIM, F_HET_CV, F_X_HET),
        HOM_REF: (F_REF, F_REF_HM, F_REF_OMIM, F_REF_CV, F_X_HOM),  # sample stats path: HOM_REF goes to X_HOM today
        UNK: (F_UNK, F_UNK_HM, F_UNK_OMIM, F_UNK_CV, F_X_UNK),
        MISSING: (F_UNK, F_UNK_HM, F_UNK_OMIM, F_UNK_CV, F_X_UNK),
    }

    qs = get_variant_queryset_for_annotation_version(annotation_version)
    qs = qs.filter(Variant.get_no_reference_q())
    qs = qs.annotate(**cgc.get_annotation_kwargs())
    qs = qs.filter(**{f"{cgc.cohortgenotype_alias}__isnull": False})

    samples_zygosity_column = f"{cgc.cohortgenotype_alias}__samples_zygosity"
    filters_column = f"{cgc.cohortgenotype_alias}__filters"
    columns = [
        "locus__contig__name",
        "locus__ref__seq",
        "alt__seq",
        "svlen",
        samples_zygosity_column,
        "variantannotation__dbsnp_rs_id",
        "variantannotation__impact",
        "variantannotation__transcript_id",
        "variantannotation__gene_id",
        "variantannotation__gene__geneannotation__omim_terms",
        "variantannotation__hgvs_g",
        "clinvar__highest_pathogenicity",
        filters_column,
    ]

    for row in qs.values_list(*columns).iterator(chunk_size=10_000):
        (chrom, ref, alt, svlen, samples_zygosity, dbsnp, impact, transcript, gene_id,
         omim, hgvs_g, clinvar, filters_value) = row
        clinvar_path = (clinvar or 0) >= 4
        ref_len = len(ref)
        alt_len = len(alt)
        is_snp = ref_len == 1 and alt_len == 1
        is_x = chrom == "X"
        impact_high_mod = impact in pathogenicity_high_mod
        has_omim = bool(omim)
        has_transcript = bool(transcript)

        # Variant-class field (snp / insertion / deletion)
        if is_snp:
            type_field = F_SNP
        elif Sequence.allele_is_symbolic(alt):
            if alt in symbolic_insertions:
                type_field = F_INS
            elif alt == VCFSymbolicAllele.DEL:
                type_field = F_DEL
            else:
                type_field = F_SNP  # fallback
        elif ref_len > alt_len:
            type_field = F_DEL
        else:
            type_field = F_INS

        # dbsnp variant-class field
        dbsnp_class_field = None
        if dbsnp:
            if is_snp:
                dbsnp_class_field = F_SNP_DBSNP
            elif ref_len > alt_len:
                dbsnp_class_field = F_DEL_DBSNP
            elif ref_len < alt_len:
                dbsnp_class_field = F_INS_DBSNP

        # ----- Per-sample accumulation (sample IS NOT NULL rows) -----
        for sample, sample_index in sample_index_pairs:
            zyg = samples_zygosity[sample_index]
            for pf_pass in (False, True):
                if pf_pass and not (has_filters and not filters_value):
                    continue
                counter = sample_counters[(sample.pk, pf_pass)]
                counter[F_VARIANT] += 1
                counter[type_field] += 1
                if dbsnp:
                    counter[F_VAR_DBSNP] += 1
                    if dbsnp_class_field is not None:
                        counter[dbsnp_class_field] += 1
                if has_transcript:
                    counter[F_GENE] += 1
                if clinvar:
                    counter[F_CLINVAR] += 1
                # Per-zygosity tally
                f_z, f_z_hm, f_z_omim, f_z_cv, f_z_x = PER_SAMPLE_TO_F.get(
                    zyg, PER_SAMPLE_TO_F[UNK])
                counter[f_z] += 1
                if impact_high_mod and zyg in (HET, HOM_ALT, HOM_REF):
                    counter[f_z_hm] += 1
                elif impact_high_mod and zyg == UNK:
                    counter[f_z_hm] += 1  # mirror today's behavior (unk also tallied)
                if has_omim:
                    counter[f_z_omim] += 1
                if clinvar_path:
                    counter[f_z_cv] += 1
                if is_x:
                    counter[f_z_x] += 1

        # ----- Aggregate accumulation (sample IS NULL, filter_key=None) -----
        agg_priority = _aggregate_zygosity_priority(samples_zygosity)
        for pf_pass in (False, True):
            if pf_pass and not (has_filters and not filters_value):
                continue
            agg_counter = aggregate_counters[(None, pf_pass)]
            agg_counter[F_VARIANT] += 1
            agg_counter[type_field] += 1
            if dbsnp:
                agg_counter[F_VAR_DBSNP] += 1
                if dbsnp_class_field is not None:
                    agg_counter[dbsnp_class_field] += 1
            if has_transcript:
                agg_counter[F_GENE] += 1
            if clinvar:
                agg_counter[F_CLINVAR] += 1
            f_z, f_z_hm, f_z_omim, f_z_cv = PRIORITY_TO_F[agg_priority]
            agg_counter[f_z] += 1
            if impact_high_mod:
                agg_counter[f_z_hm] += 1
            if has_omim:
                agg_counter[f_z_omim] += 1
            if clinvar_path:
                agg_counter[f_z_cv] += 1
            if is_x:
                if agg_priority == HOM_ALT or agg_priority == HOM_REF:
                    agg_counter[F_X_HOM] += 1
                elif agg_priority == HET:
                    agg_counter[F_X_HET] += 1
                else:
                    agg_counter[F_X_UNK] += 1

        # ----- Filter-keyed accumulation (sample IS NULL, filter_key=...) -----
        # Trio inheritance modes
        if trio_indexes is not None:
            mum_zyg = samples_zygosity[trio_indexes[0]]
            dad_zyg = samples_zygosity[trio_indexes[1]]
            proband_zyg = samples_zygosity[trio_indexes[2]]
            modes = _trio_predicates(mum_zyg, dad_zyg, proband_zyg, chrom,
                                     trio.mother_affected, trio.father_affected)
            # Compound het: defer; record gene-source so we can resolve once all variants seen
            chet_match = None
            if proband_zyg == HET:
                if mum_zyg == HET and dad_zyg in NO_VARIANT:
                    chet_match = "mum"
                elif mum_zyg in NO_VARIANT and dad_zyg == HET:
                    chet_match = "dad"
            for pf_pass in (False, True):
                if pf_pass and not (has_filters and not filters_value):
                    continue
                # Per-mode buckets — counted just like aggregate, classified by proband zyg
                proband_priority = (HOM_ALT if proband_zyg == HOM_ALT
                                    else HET if proband_zyg == HET
                                    else HOM_REF if proband_zyg == HOM_REF
                                    else UNK)
                f_z, f_z_hm, f_z_omim, f_z_cv = PRIORITY_TO_F[proband_priority]
                for mode in modes:
                    key = canonical_filter_key({"inheritance": mode})
                    if key not in precompute_keys:
                        continue
                    bucket = aggregate_counters[(key, pf_pass)]
                    bucket[F_VARIANT] += 1
                    bucket[type_field] += 1
                    if dbsnp:
                        bucket[F_VAR_DBSNP] += 1
                        if dbsnp_class_field is not None:
                            bucket[dbsnp_class_field] += 1
                    if has_transcript:
                        bucket[F_GENE] += 1
                    if clinvar:
                        bucket[F_CLINVAR] += 1
                    bucket[f_z] += 1
                    if impact_high_mod:
                        bucket[f_z_hm] += 1
                    if has_omim:
                        bucket[f_z_omim] += 1
                    if clinvar_path:
                        bucket[f_z_cv] += 1
                # Compound het: queue the contribution
                if chet_match is not None and gene_id is not None:
                    chet_gene_sources[(gene_id, pf_pass)].add(chet_match)
                    chet_pending.append((
                        gene_id, pf_pass, type_field, dbsnp, dbsnp_class_field,
                        has_transcript, clinvar, clinvar_path, impact_high_mod, has_omim,
                        proband_priority,
                    ))

    # Compound het post-processing — only genes with both "mum" and "dad"
    # sources contribute. Add to the compound_het filter_key bucket.
    chet_key = canonical_filter_key({"inheritance": "compound_het"})
    if chet_key in precompute_keys:
        for entry in chet_pending:
            (gene_id, pf_pass, type_field, dbsnp, dbsnp_class_field,
             has_transcript, clinvar, clinvar_path, impact_high_mod, has_omim,
             proband_priority) = entry
            if {"mum", "dad"} <= chet_gene_sources[(gene_id, pf_pass)]:
                bucket = aggregate_counters[(chet_key, pf_pass)]
                bucket[F_VARIANT] += 1
                bucket[type_field] += 1
                if dbsnp:
                    bucket[F_VAR_DBSNP] += 1
                    if dbsnp_class_field is not None:
                        bucket[dbsnp_class_field] += 1
                if has_transcript:
                    bucket[F_GENE] += 1
                if clinvar:
                    bucket[F_CLINVAR] += 1
                f_z, f_z_hm, f_z_omim, f_z_cv = PRIORITY_TO_F[proband_priority]
                bucket[f_z] += 1
                if impact_high_mod:
                    bucket[f_z_hm] += 1
                if has_omim:
                    bucket[f_z_omim] += 1
                if clinvar_path:
                    bucket[f_z_cv] += 1

    _persist_cohort_stats(cgc, annotation_version, code_version, has_filters,
                          sample_counters, aggregate_counters, precompute_keys,
                          [s for s, _ in sample_index_pairs])


def _cohort_has_filters(cohort: Cohort) -> bool:
    """ A cohort has 'filters' iff its underlying VCF (if any) does. Custom
        multi-VCF cohorts inherit filtering from their member VCFs but, for
        passing-filter stats purposes, we mirror the per-VCF flag of the
        cohort's base VCF if any. """
    vcf = cohort.get_vcf()
    return bool(vcf and vcf.has_filters)


def _persist_cohort_stats(cgc, annotation_version, code_version, has_filters,
                          sample_counters, aggregate_counters, precompute_keys, samples):
    """ Replace any existing rows for this CGC + annotation versions on the
        new models with the computed counters. """
    vav = annotation_version.variant_annotation_version
    gav = annotation_version.gene_annotation_version
    cv = annotation_version.clinvar_version

    # Wipe existing rows for this CGC scoped by annotation version where applicable
    CohortGenotypeStats.objects.filter(cohort_genotype_collection=cgc).delete()
    CohortGenotypeVariantAnnotationStats.objects.filter(
        cohort_genotype_collection=cgc, variant_annotation_version=vav).delete()
    CohortGenotypeGeneAnnotationStats.objects.filter(
        cohort_genotype_collection=cgc, gene_annotation_version=gav).delete()
    CohortGenotypeClinVarAnnotationStats.objects.filter(
        cohort_genotype_collection=cgc, clinvar_version=cv).delete()

    rows_genotype = []
    rows_variant = []
    rows_gene = []
    rows_clinvar = []

    def _populate(target, counter, fields):
        for f in fields:
            setattr(target, f, counter[_FIELD_INDEX[f]])

    GENOTYPE_FIELDS = (
        "variant_count", "snp_count", "insertions_count", "deletions_count",
        "ref_count", "het_count", "hom_count", "unk_count",
        "x_hom_count", "x_het_count", "x_unk_count",
    )
    VARIANT_ANN_FIELDS = (
        "variant_dbsnp_count", "snp_dbsnp_count", "insertions_dbsnp_count", "deletions_dbsnp_count",
        "ref_high_or_moderate_count", "het_high_or_moderate_count",
        "hom_high_or_moderate_count", "unk_high_or_moderate_count",
    )
    GENE_ANN_FIELDS = (
        "gene_count",
        "ref_omim_phenotype_count", "het_omim_phenotype_count",
        "hom_omim_phenotype_count", "unk_omim_phenotype_count",
    )
    CLINVAR_ANN_FIELDS = (
        "clinvar_count",
        "ref_clinvar_pathogenic_count", "het_clinvar_pathogenic_count",
        "hom_clinvar_pathogenic_count", "unk_clinvar_pathogenic_count",
    )

    sample_by_pk = {s.pk: s for s in samples}

    # Per-sample rows (filter_key NULL)
    for (sample_pk, pf_pass), counter in sample_counters.items():
        if pf_pass and not has_filters:
            continue
        sample = sample_by_pk[sample_pk]
        gs = CohortGenotypeStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
        )
        _populate(gs, counter, GENOTYPE_FIELDS)
        rows_genotype.append(gs)
        va = CohortGenotypeVariantAnnotationStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
        )
        _populate(va, counter, VARIANT_ANN_FIELDS)
        rows_variant.append(va)
        ga = CohortGenotypeGeneAnnotationStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
        )
        _populate(ga, counter, GENE_ANN_FIELDS)
        rows_gene.append(ga)
        ca = CohortGenotypeClinVarAnnotationStats(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
        )
        _populate(ca, counter, CLINVAR_ANN_FIELDS)
        rows_clinvar.append(ca)

    # Aggregate / filter-keyed rows (sample NULL)
    for (filter_key, pf_pass), counter in aggregate_counters.items():
        if pf_pass and not has_filters:
            continue
        gs = CohortGenotypeStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
        )
        _populate(gs, counter, GENOTYPE_FIELDS)
        rows_genotype.append(gs)
        va = CohortGenotypeVariantAnnotationStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
        )
        _populate(va, counter, VARIANT_ANN_FIELDS)
        rows_variant.append(va)
        ga = CohortGenotypeGeneAnnotationStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
        )
        _populate(ga, counter, GENE_ANN_FIELDS)
        rows_gene.append(ga)
        ca = CohortGenotypeClinVarAnnotationStats(
            cohort_genotype_collection=cgc, sample=None, filter_key=filter_key,
            passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
        )
        _populate(ca, counter, CLINVAR_ANN_FIELDS)
        rows_clinvar.append(ca)

    # Ensure every precomputed (key, pf_pass) bucket has a row even if no
    # variants matched it (so reads hit instead of falling back to live count).
    existing_agg_keys = {(fk, pf) for fk, pf in aggregate_counters.keys()}
    for fk in precompute_keys:
        for pf_pass in ((False, True) if has_filters else (False,)):
            if (fk, pf_pass) in existing_agg_keys:
                continue
            zero = _new_counter()
            rows_genotype.append(CohortGenotypeStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, import_status=ImportStatus.SUCCESS,
            ))
            rows_variant.append(CohortGenotypeVariantAnnotationStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, variant_annotation_version=vav,
            ))
            rows_gene.append(CohortGenotypeGeneAnnotationStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, gene_annotation_version=gav,
            ))
            rows_clinvar.append(CohortGenotypeClinVarAnnotationStats(
                cohort_genotype_collection=cgc, sample=None, filter_key=fk,
                passing_filter=pf_pass, code_version=code_version, clinvar_version=cv,
            ))
            del zero  # all-zero row, populated by defaults

    CohortGenotypeStats.objects.bulk_create(rows_genotype, batch_size=2000)
    CohortGenotypeVariantAnnotationStats.objects.bulk_create(rows_variant, batch_size=2000)
    CohortGenotypeGeneAnnotationStats.objects.bulk_create(rows_gene, batch_size=2000)
    CohortGenotypeClinVarAnnotationStats.objects.bulk_create(rows_clinvar, batch_size=2000)


def _get_variant_class(alt: str, hgvs_g: str) -> Optional[VariantClass]:
    """ The VEP VariantClass is not always accurate so basically DIY here """

    variant_type = None
    if Sequence.allele_is_symbolic(alt):
        if alt == VCFSymbolicAllele.INV:
            variant_type = VariantClass.INVERSION
        elif alt == VCFSymbolicAllele.DEL:
            variant_type = VariantClass.DELETION
        elif alt == VCFSymbolicAllele.DUP:
            variant_type = VariantClass.DUPLICATION
    elif hgvs_g:
        if "delins" in hgvs_g:
            variant_type = VariantClass.INDEL
        elif "inv" in hgvs_g:
            variant_type = VariantClass.INVERSION
        elif "ins" in hgvs_g:
            variant_type = VariantClass.INSERTION
        elif "del" in hgvs_g:
            variant_type = VariantClass.DELETION
        elif "dup" in hgvs_g:
            variant_type = VariantClass.DUPLICATION
    return variant_type


def calculate_needed_stats(run_async=False):
    """ Works out what needs to be (re)calculated and does so. A VCF needs
        recompute when its cohort's CGC has no genotype-level stats row, or
        when any of the three annotation versions on file no longer match
        AnnotationVersion.latest. """

    logging.info("Deleting CohortGenotypeStats with non-SUCCESS import_status (leftovers)")
    deleted = CohortGenotypeStats.objects.exclude(import_status=ImportStatus.SUCCESS).delete()
    logging.info(deleted)

    logging.info("Calculating cohort stats (run_async=%s)", run_async)
    for genome_build in GenomeBuild.builds_with_annotation():
        try:
            annotation_version = AnnotationVersion.latest(genome_build)
        except InvalidAnnotationVersionError:
            logging.info(f"Skipping calculating sample stats for incomplete annotation version for build {genome_build}")
            continue

        vav = annotation_version.variant_annotation_version
        gav = annotation_version.gene_annotation_version
        cv = annotation_version.clinvar_version

        # A VCF's stats are stale if any of the per-sample (sample IS NOT NULL,
        # filter_key NULL) rows for any of the four stat classes is missing
        # for the latest annotation versions.
        needs_stats = Q(cohort__cohortgenotypecollection__genotype_stats__isnull=True)
        needs_stats |= ~Q(
            cohort__cohortgenotypecollection__variant_annotation_stats__variant_annotation_version=vav)
        needs_stats |= ~Q(
            cohort__cohortgenotypecollection__gene_annotation_stats__gene_annotation_version=gav)
        needs_stats |= ~Q(
            cohort__cohortgenotypecollection__clinvar_annotation_stats__clinvar_version=cv)

        vcf_qs = VCF.objects.filter(needs_stats, genome_build=genome_build).distinct()
        logging.info("Build: %s VCFs needing stats: %d", genome_build, vcf_qs.count())
        for vcf in vcf_qs:
            task = calculate_vcf_stats.si(vcf.pk, annotation_version.pk)  # @UndefinedVariable
            if run_async:
                task.apply_async()
            else:
                result = task.apply()
                if result.successful():
                    logging.info("Successfully calculated stats for %s", vcf)
                else:
                    logging.error("Died for VCF %s: %s", vcf, result.result)
