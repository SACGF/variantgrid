import logging
import time
from collections import defaultdict
from typing import Optional

import celery
import numpy as np
from django.conf import settings
from django.db.models.query_utils import Q

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models import SampleVariantAnnotationStats, SampleGeneAnnotationStats, \
    SampleClinVarAnnotationStats, SampleVariantAnnotationStatsPassingFilter, \
    SampleGeneAnnotationStatsPassingFilter, SampleClinVarAnnotationStatsPassingFilter, \
    AnnotationVersion
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import InvalidAnnotationVersionError, VCFAnnotationStats
from eventlog.models import create_event
from library.django_utils import thread_safe_unique_together_get_or_create
from library.enums.log_level import LogLevel
from library.genomics.vcf_enums import VariantClass, VCFSymbolicAllele
from library.git import Git
from library.log_utils import get_traceback
from snpdb.models import Sample, SampleStats, ImportStatus, SampleStatsPassingFilter, VCF, Variant, \
    SampleStatsCodeVersion, Sequence, VCFLengthStatsCollection, VCFLengthStats
from snpdb.models import Zygosity
from snpdb.models.models_genome import GenomeBuild

SAMPLE_STATS = "sample_stats"
SAMPLE_ANNOTATION_STATS = "annotation_stats"
SAMPLE_GENE_STATS = "gene_stats"
SAMPLE_CLINVAR_STATS = "clinvar_stats"


@celery.shared_task
def calculate_vcf_stats(vcf_id, annotation_version_id):
    vcf = VCF.objects.get(pk=vcf_id)
    try:
        annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
        msg = f"Build for annotation: {annotation_version.genome_build} and VCF: {vcf.genome_build} match"
        assert annotation_version.genome_build == vcf.genome_build, msg
        _actually_calculate_vcf_stats(vcf, annotation_version)

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


def _create_stats_per_sample(vcf: VCF, annotation_version) -> tuple[dict, Optional[dict]]:
    vav_kwargs = {"variant_annotation_version": annotation_version.variant_annotation_version}
    ega_kwargs = {"gene_annotation_version": annotation_version.gene_annotation_version}
    cv_kwargs = {"clinvar_version": annotation_version.clinvar_version}

    code_version = _get_sample_stats_code_version()

    STATS_CLASSES = {
        SAMPLE_STATS: (SampleStats, SampleStatsPassingFilter, {}),
        SAMPLE_ANNOTATION_STATS: (SampleVariantAnnotationStats, SampleVariantAnnotationStatsPassingFilter, vav_kwargs),
        SAMPLE_GENE_STATS: (SampleGeneAnnotationStats, SampleGeneAnnotationStatsPassingFilter, ega_kwargs),
        SAMPLE_CLINVAR_STATS: (SampleClinVarAnnotationStats, SampleClinVarAnnotationStatsPassingFilter, cv_kwargs),
    }

    stats_per_sample = defaultdict(dict)
    stats_passing_filters_per_sample = defaultdict(dict)

    samples = list(vcf.sample_set.all())
    for st, (klass, pf_klass, extra_kwargs) in STATS_CLASSES.items():
        # We want to delete and recreate them so they reset to zero
        klass.objects.filter(sample__in=samples).delete()
        pf_klass.objects.filter(sample__in=samples).delete()
        for sample in samples:
            stats_per_sample[sample][st] = klass.objects.create(sample=sample, code_version=code_version,
                                                                **extra_kwargs)
            if vcf.has_filters:
                stats_passing_filters_per_sample[sample][st] = pf_klass.objects.create(sample=sample,
                                                                                       code_version=code_version,
                                                                                       **extra_kwargs)
    return stats_per_sample, stats_passing_filters_per_sample


def _actually_calculate_vcf_stats(vcf: VCF, annotation_version: AnnotationVersion):
    start = time.time()
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    qs = qs.filter(Variant.get_no_reference_q())
    qs = vcf.get_variant_qs(qs)
    cgc = vcf.cohort.cohort_genotype_collection
    samples_zygosity_column = f"{cgc.cohortgenotype_alias}__samples_zygosity"
    columns = [
        "locus__contig__name",
        "locus__ref__seq",
        "alt__seq",
        "svlen",
        samples_zygosity_column,
        "variantannotation__dbsnp_rs_id",
        "variantannotation__impact",
        "variantannotation__transcript_id",
        "variantannotation__gene__geneannotation__omim_terms",
        "variantannotation__vep_skipped_reason",
        "variantannotation__hgvs_g",
        "clinvar__highest_pathogenicity",
        f"{cgc.cohortgenotype_alias}__filters",
    ]
    cohort_samples = list(vcf.cohort.get_cohort_samples())
    vep_skipped_count = 0
    total_delta = 0
    variant_class_lengths = defaultdict(list)

    # bind constants and functions locally
    pathogenicity_high_mod = {PathogenicityImpact.HIGH, PathogenicityImpact.MODERATE}
    symbolic_insertions = {VCFSymbolicAllele.DUP, VCFSymbolicAllele.INS}
    HET, HOM_ALT, HOM_REF, UNK = Zygosity.HET, Zygosity.HOM_ALT, Zygosity.HOM_REF, Zygosity.UNKNOWN_ZYGOSITY
    sample_and_index = []
    for cohort_sample in cohort_samples:
        sample_and_index.append((cohort_sample.sample, cohort_sample.cohort_genotype_packed_field_index))

    # Use array indexes to avoid dict lookups
    F_CLINVAR_COUNT = 0
    F_DELETIONS_COUNT = 1
    F_DELETIONS_DBSNP_COUNT = 2
    F_GENE_COUNT = 3
    F_HET_CLINVAR_PATHOGENIC_COUNT = 4
    F_HET_COUNT = 5
    F_HET_HIGH_OR_MODERATE_COUNT = 6
    F_HET_OMIM_PHENOTYPE_COUNT = 7
    F_HOM_CLINVAR_PATHOGENIC_COUNT = 8
    F_HOM_COUNT = 9
    F_HOM_HIGH_OR_MODERATE_COUNT = 10
    F_HOM_OMIM_PHENOTYPE_COUNT = 11
    F_INSERTIONS_COUNT = 12
    F_INSERTIONS_DBSNP_COUNT = 13
    F_REF_CLINVAR_PATHOGENIC_COUNT = 14
    F_REF_COUNT = 15
    F_REF_HIGH_OR_MODERATE_COUNT = 16
    F_REF_OMIM_PHENOTYPE_COUNT = 17
    F_SNP_COUNT = 18
    F_SNP_DBSNP_COUNT = 19
    F_UNK_CLINVAR_PATHOGENIC_COUNT = 20
    F_UNK_COUNT = 21
    F_UNK_OMIM_PHENOTYPE_COUNT = 22
    F_VARIANT_COUNT = 23
    F_VARIANT_DBSNP_COUNT = 24
    F_X_HET_COUNT = 25
    F_X_HOM_COUNT = 26
    F_X_UNK_COUNT = 27
    NUM_FIELDS = 28

    # Use Counter instead of Django models as it has lower overhead
    sample_counts = defaultdict(lambda: [0] * NUM_FIELDS)
    sample_counts_passing_filters = defaultdict(lambda: [0] * NUM_FIELDS)

    for row in qs.values_list(*columns).iterator(chunk_size=10_000):
        (chrom, ref, alt, svlen, samples_zygosity, dbsnp, impact, transcript, omim, vep_skipped, hgvs_g, clinvar, filters) = row
        clinvar_highest_pathogenicity = clinvar or 0

        ref_len = len(ref)
        alt_len = len(alt)
        is_snp = ref_len == 1 and alt_len == 1
        impact_high_or_mod = impact in pathogenicity_high_mod

        if vep_skipped:
            vep_skipped_count += 1

        variant_type_field = F_SNP_COUNT
        if not is_snp:
            variant_class = _get_variant_class(alt, hgvs_g)
            if svlen:
                variant_length = abs(svlen)  # Need to be positive if we ever take log
            else:
                variant_length = max(alt_len, ref_len)
            variant_class_lengths[variant_class].append(variant_length)

            if Sequence.allele_is_symbolic(alt):
                if alt in symbolic_insertions:
                    variant_type_field = F_INSERTIONS_COUNT
                elif alt == VCFSymbolicAllele.DEL:
                    variant_type_field = F_DELETIONS_COUNT
            else:
                if ref_len > alt_len:
                    variant_type_field = F_DELETIONS_COUNT
                else:
                    variant_type_field = F_INSERTIONS_COUNT

        zygosity_updates = {
            HET: [F_HET_COUNT],
            HOM_ALT: [F_HOM_COUNT],
            HOM_REF: [F_REF_COUNT],
            UNK: [F_UNK_COUNT],
        }
        if chrom == 'X':
            zygosity_updates[HET].append(F_X_HET_COUNT)
            zygosity_updates[HOM_ALT].append(F_X_HOM_COUNT)
            zygosity_updates[HOM_REF].append(F_X_HOM_COUNT)
            zygosity_updates[UNK].append(F_X_UNK_COUNT)
        if impact_high_or_mod:
            zygosity_updates[HET].append(F_HET_HIGH_OR_MODERATE_COUNT)
            zygosity_updates[HOM_ALT].append(F_HOM_HIGH_OR_MODERATE_COUNT)
            zygosity_updates[HOM_REF].append(F_REF_HIGH_OR_MODERATE_COUNT)
        if omim:
            zygosity_updates[HET].append(F_HET_OMIM_PHENOTYPE_COUNT)
            zygosity_updates[HOM_ALT].append(F_HOM_OMIM_PHENOTYPE_COUNT)
            zygosity_updates[HOM_REF].append(F_REF_OMIM_PHENOTYPE_COUNT)
            zygosity_updates[UNK].append(F_UNK_OMIM_PHENOTYPE_COUNT)
        if clinvar_highest_pathogenicity >= 4:
            zygosity_updates[HET].append(F_HET_CLINVAR_PATHOGENIC_COUNT)
            zygosity_updates[HOM_ALT].append(F_HOM_CLINVAR_PATHOGENIC_COUNT)
            zygosity_updates[HOM_REF].append(F_REF_CLINVAR_PATHOGENIC_COUNT)
            zygosity_updates[UNK].append(F_UNK_CLINVAR_PATHOGENIC_COUNT)

        variant_updates = [F_VARIANT_COUNT, variant_type_field]
        if dbsnp:
            variant_updates.append(F_VARIANT_DBSNP_COUNT)

            if ref_len == 1 and alt_len == 1:
                variant_updates.append(F_SNP_DBSNP_COUNT)
            elif ref_len > alt_len:
                variant_updates.append(F_DELETIONS_DBSNP_COUNT)
            elif ref_len < alt_len:
                variant_updates.append(F_INSERTIONS_DBSNP_COUNT)

        if transcript:
            variant_updates.append(F_GENE_COUNT)

        if clinvar_highest_pathogenicity:
            variant_updates.append(F_CLINVAR_COUNT)

        for sample, sample_index in sample_and_index:
            zygosity = samples_zygosity[sample_index]

            update_indexes = variant_updates + zygosity_updates[zygosity]
            sc = sample_counts[sample]
            for i in update_indexes:
                sc[i] += 1
            # Update if passed
            if vcf.has_filters and not filters:
                sc_pf = sample_counts_passing_filters[sample]
                for i in update_indexes:
                    sc_pf[i] += 1
            # End of per-sample stuff

    copy_start = time.time()
    print(f"calc time: {copy_start - start}")

    stats_per_sample, stats_passing_filters_per_sample = _create_stats_per_sample(vcf, annotation_version)
    counter_and_models = [
        (sample_counts, stats_per_sample),
        (sample_counts_passing_filters, stats_passing_filters_per_sample),
    ]

    FIELD_LOOKUPS = {
        'clinvar_count': F_CLINVAR_COUNT,
        'deletions_count': F_DELETIONS_COUNT,
        'deletions_dbsnp_count': F_DELETIONS_DBSNP_COUNT,
        'gene_count': F_GENE_COUNT,
        'het_clinvar_pathogenic_count': F_HET_CLINVAR_PATHOGENIC_COUNT,
        'het_count': F_HET_COUNT,
        'het_high_or_moderate_count': F_HET_HIGH_OR_MODERATE_COUNT,
        'het_omim_phenotype_count': F_HET_OMIM_PHENOTYPE_COUNT,
        'hom_clinvar_pathogenic_count': F_HOM_CLINVAR_PATHOGENIC_COUNT,
        'hom_count': F_HOM_COUNT,
        'hom_high_or_moderate_count': F_HOM_HIGH_OR_MODERATE_COUNT,
        'hom_omim_phenotype_count': F_HOM_OMIM_PHENOTYPE_COUNT,
        'insertions_count': F_INSERTIONS_COUNT,
        'insertions_dbsnp_count': F_INSERTIONS_DBSNP_COUNT,
        'ref_clinvar_pathogenic_count': F_REF_CLINVAR_PATHOGENIC_COUNT,
        'ref_count': F_REF_COUNT,
        'ref_high_or_moderate_count': F_REF_HIGH_OR_MODERATE_COUNT,
        'ref_omim_phenotype_count': F_REF_OMIM_PHENOTYPE_COUNT,
        'snp_count': F_SNP_COUNT,
        'snp_dbsnp_count': F_SNP_DBSNP_COUNT,
        'unk_clinvar_pathogenic_count': F_UNK_CLINVAR_PATHOGENIC_COUNT,
        'unk_count': F_UNK_COUNT,
        'unk_omim_phenotype_count': F_UNK_OMIM_PHENOTYPE_COUNT,
        'variant_count': F_VARIANT_COUNT,
        'variant_dbsnp_count': F_VARIANT_DBSNP_COUNT,
        'x_het_count': F_X_HET_COUNT,
        'x_hom_count': F_X_HOM_COUNT,
        'x_unk_count': F_X_UNK_COUNT,
    }
    # Count is a big mix of everything, copy off into appropriate models
    for sc, sample_stats_models in counter_and_models:
        for sample, stats_models in sample_stats_models.items():
            stats_models[SAMPLE_STATS].import_status = ImportStatus.SUCCESS  # Only SampleStats has import_status

            for stats_model in stats_models.values():
                counts = sc[sample]
                for field_name, i in FIELD_LOOKUPS.items():
                    count = counts[i]
                    if hasattr(stats_model, field_name):
                        setattr(stats_model, field_name, count)

                stats_model.save()

    if vep_skipped_count:
        # We may be re-running a part of sample stats. Store the highest number skipped
        stats, created = VCFAnnotationStats.objects.get_or_create(vcf=vcf, variant_annotation_version=annotation_version.variant_annotation_version,
                                                                  defaults={"vep_skipped_count": vep_skipped_count})
        if not created:
            if vep_skipped_count > stats.vep_skipped_count:
                stats.vep_skipped_count = vep_skipped_count
                stats.save()

    if variant_class_lengths:
        code_version = _get_sample_stats_code_version()
        collection, created = VCFLengthStatsCollection.objects.update_or_create(vcf=vcf,
                                                                                defaults={"code_version": code_version})
        if not created:
            # Clear old ones (we're reloading now)
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

            vcf_ls = VCFLengthStats(collection=collection,
                                    variant_class=variant_class,
                                    is_log=is_log,
                                    histogram_counts=histogram_counts.tolist(),
                                    histogram_bin_edges=histogram_bin_edges.tolist())
            vcf_length_stats.append(vcf_ls)
        if vcf_length_stats:
            VCFLengthStats.objects.bulk_create(vcf_length_stats)

    end = time.time()
    logging.info("SampleStats for %s took %.2f seconds", vcf, end - start)


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
    """ Works out what needs to be calculated and does so """

    logging.info("Deleting Sample Stats where import_status != SUCCESS (leftovers)")
    deleted = SampleStats.objects.exclude(import_status=ImportStatus.SUCCESS).delete()
    logging.info(deleted)

    logging.info("Calculating Sample Stats (run_async=%s)", run_async)
    for genome_build in GenomeBuild.builds_with_annotation():
        try:
            annotation_version = AnnotationVersion.latest(genome_build)
        except InvalidAnnotationVersionError:
            logging.info(f"Skipping calculating sample stats for incomplete annotation version for build {genome_build}")
            continue

        needs_stats = Q(samplestats__isnull=True)
        needs_stats |= ~Q(samplevariantannotationstats__variant_annotation_version=annotation_version.variant_annotation_version)
        needs_stats |= ~Q(samplegeneannotationstats__gene_annotation_version=annotation_version.gene_annotation_version)
        needs_stats |= ~Q(sampleclinvarannotationstats__clinvar_version=annotation_version.clinvar_version)

        needs_stats_passing_filters = Q(samplestatspassingfilter__isnull=True)
        needs_stats_passing_filters |= ~Q(samplevariantannotationstats__variant_annotation_version=annotation_version.variant_annotation_version)
        needs_stats_passing_filters |= ~Q(samplegeneannotationstats__gene_annotation_version=annotation_version.gene_annotation_version)
        needs_stats_passing_filters |= ~Q(sampleclinvarannotationstats__clinvar_version=annotation_version.clinvar_version)

        qs_filter = needs_stats | (Q(vcf__vcffilter__isnull=False) & needs_stats_passing_filters)

        samples_qs = Sample.objects.filter(qs_filter, vcf__genome_build=genome_build).distinct()

        if run_async:
            logging.info("Launching sample stats jobs asynchronously")

        vcf_qs = VCF.objects.filter(sample__in=samples_qs).distinct()
        logging.info(f"Build: %s Samples: %d in %d VCFs", genome_build, samples_qs.count(), vcf_qs.count())
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
