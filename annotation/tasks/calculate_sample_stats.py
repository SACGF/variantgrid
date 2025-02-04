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
    """
    code_version, _ = thread_safe_unique_together_get_or_create(SampleStatsCodeVersion,
                                                                name="SampleStats",
                                                                version=2,
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
        "variantannotation__gene__geneannotation__omim_terms",
        "variantannotation__vep_skipped_reason",
        "variantannotation__hgvs_g",
        "clinvar__highest_pathogenicity",
    ]

    if vcf.has_filters:
        columns.append(filters_column)

    values_queryset = qs.values(*columns)
    cohort_samples = list(vcf.cohort.get_cohort_samples())
    stats_per_sample, stats_passing_filters_per_sample = _create_stats_per_sample(vcf, annotation_version)
    vep_skipped_count = 0
    variant_class_lengths = defaultdict(list)

    for vals in values_queryset.iterator():
        chrom = vals["locus__contig__name"]
        ref = vals["locus__ref__seq"]
        alt = vals["alt__seq"]
        svlen = vals["svlen"]
        samples_zygosity = vals[samples_zygosity_column]
        dbsnp = vals["variantannotation__dbsnp_rs_id"]
        impact = vals["variantannotation__impact"]
        transcript = vals["variantannotation__transcript_id"]
        omim = vals["variantannotation__gene__geneannotation__omim_terms"]
        vep_skipped = vals["variantannotation__vep_skipped_reason"]
        clinvar_highest_pathogenicity = vals["clinvar__highest_pathogenicity"]
        hgvs_g = vals["variantannotation__hgvs_g"]

        ref_len = len(ref)
        alt_len = len(alt)
        is_snp = ref_len == 1 and alt_len == 1
        impact_high_or_mod = impact in {PathogenicityImpact.HIGH, PathogenicityImpact.MODERATE}

        stats_list = [stats_per_sample]

        if vcf.has_filters:
            filters = vals[filters_column]
            if not filters:
                stats_list.append(stats_passing_filters_per_sample)

        if vep_skipped:
            vep_skipped_count += 1

        if not is_snp:
            variant_class = _get_variant_class(alt, hgvs_g)
            if svlen:
                variant_length = abs(svlen)  # Need to be positive if we ever take log
            else:
                variant_length = max(alt_len, ref_len)
            variant_class_lengths[variant_class].append(variant_length)

        # Sample related stuff

        for cohort_sample in cohort_samples:
            zygosity = samples_zygosity[cohort_sample.cohort_genotype_packed_field_index]
            sample_stats_list = [s[cohort_sample.sample] for s in stats_list]

            for ss in [sl[SAMPLE_STATS] for sl in sample_stats_list]:
                ss.variant_count += 1

                if Sequence.allele_is_symbolic(alt):
                    if alt in (VCFSymbolicAllele.DUP, VCFSymbolicAllele.INS):
                        ss.insertions_count += 1
                    elif alt == VCFSymbolicAllele.DEL:
                        ss.deletions_count += 1
                else:
                    if is_snp:
                        ss.snp_count += 1
                    elif ref_len > alt_len:
                        ss.deletions_count += 1
                    elif ref_len < alt_len:
                        ss.insertions_count += 1

                if zygosity == Zygosity.HET:
                    ss.het_count += 1
                elif zygosity == Zygosity.HOM_ALT:
                    ss.hom_count += 1
                elif zygosity == Zygosity.HOM_REF:
                    ss.ref_count += 1
                elif zygosity == Zygosity.UNKNOWN_ZYGOSITY:
                    ss.unk_count += 1

                if chrom == 'X':
                    if zygosity == Zygosity.HET:
                        ss.x_het_count += 1
                    elif zygosity in [Zygosity.HOM_REF, Zygosity.HOM_ALT]:
                        ss.x_hom_count += 1

            for ast in [sl[SAMPLE_ANNOTATION_STATS] for sl in sample_stats_list]:
                if dbsnp:
                    ast.variant_dbsnp_count += 1

                    if ref_len == 1 and alt_len == 1:
                        ast.snp_dbsnp_count += 1
                    elif ref_len > alt_len:
                        ast.deletions_dbsnp_count += 1
                    elif ref_len < alt_len:
                        ast.insertions_dbsnp_count += 1

                if impact_high_or_mod:
                    if zygosity == Zygosity.HET:
                        ast.het_high_or_moderate_count += 1
                    elif zygosity == Zygosity.HOM_ALT:
                        ast.hom_high_or_moderate_count += 1
                    elif zygosity == Zygosity.HOM_REF:
                        ast.ref_high_or_moderate_count += 1

            if transcript:
                for gs in [sl[SAMPLE_GENE_STATS] for sl in sample_stats_list]:
                    gs.gene_count += 1

                    if omim:
                        if zygosity == Zygosity.HET:
                            gs.het_omim_phenotype_count += 1
                        elif zygosity == Zygosity.HOM_ALT:
                            gs.hom_omim_phenotype_count += 1
                        elif zygosity == Zygosity.HOM_REF:
                            gs.ref_omim_phenotype_count += 1
                        elif zygosity == Zygosity.UNKNOWN_ZYGOSITY:
                            gs.unk_omim_phenotype_count += 1

            if clinvar_highest_pathogenicity:
                for cs in [sl[SAMPLE_CLINVAR_STATS] for sl in sample_stats_list]:
                    cs.clinvar_count += 1
                    if clinvar_highest_pathogenicity >= 4:
                        if zygosity == Zygosity.HET:
                            cs.het_clinvar_pathogenic_count += 1
                        elif zygosity == Zygosity.HOM_ALT:
                            cs.hom_clinvar_pathogenic_count += 1
                        elif zygosity == Zygosity.HOM_REF:
                            cs.ref_clinvar_pathogenic_count += 1
                        elif zygosity == Zygosity.UNKNOWN_ZYGOSITY:
                            cs.unk_clinvar_pathogenic_count += 1

    for stats_dict in [stats_per_sample, stats_passing_filters_per_sample or {}]:
        for stats in stats_dict.values():
            stats[SAMPLE_STATS].import_status = ImportStatus.SUCCESS
            for s in stats.values():
                s.save()

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
