import logging

from celery.canvas import Signature
from django.db.models import QuerySet
from django.dispatch import receiver

from classification.models import Classification, variants_classification_changed_signal
from classification.models.classification_import_run import ClassificationImportRun, \
    classification_imports_complete_signal
from snpdb.common_variants import get_classified_high_frequency_variants_qs, get_common_filters_in_use
from snpdb.models import CohortGenotypeCommonFilterVersion, Variant, VariantAllele

"""
Classifying a common variant means it can no longer be stored in the common partition, so it has to be moved
to the uncommon one - see snpdb.tasks.cohort_genotype_tasks.common_variant_classified_task
"""


def _launch_common_variant_classified_tasks(cgcfv: CohortGenotypeCommonFilterVersion, variant_qs: QuerySet[Variant]):
    """ common_variant_classified_task creates the CommonVariantClassified, so only call it for variants that
        don't have one yet """
    for variant in variant_qs.exclude(commonvariantclassified__common_filter=cgcfv):
        task_name = "snpdb.tasks.cohort_genotype_tasks.common_variant_classified_task"
        task = Signature(task_name, args=(variant.pk, cgcfv.pk), immutable=True)
        task.apply_async()


@receiver(variants_classification_changed_signal, sender=Classification)
def variants_classification_changed(sender, **kwargs):  # pylint: disable=unused-argument
    variants = kwargs['variants']
    genome_build = kwargs['genome_build']

    if ongoing_import := ClassificationImportRun.ongoing_imports():
        # Each check scans every classified high frequency variant in the build, which is far too expensive to
        # do per-classification during a bulk import - classification_imports_complete does one sweep instead
        logging.info("variants_classification_changed_signal delayed due to %s", ongoing_import)
        return

    logging.info("variants_classification_changed_signal!! %s, %s", genome_build, variants)

    # Look to see if any of these are used in a common filter (and not already handled)
    for cgcfv in get_common_filters_in_use(genome_build):
        va_qs = VariantAllele.objects.filter(variant__in=variants,
                                             genome_build=genome_build)
        va_qs = va_qs.exclude(variant__commonvariantclassified__common_filter=cgcfv)
        alleles = va_qs.values_list("allele")
        _launch_common_variant_classified_tasks(cgcfv, get_classified_high_frequency_variants_qs(cgcfv, alleles=alleles))


@receiver(classification_imports_complete_signal, sender=ClassificationImportRun)
def classification_imports_complete(sender, **kwargs):  # pylint: disable=unused-argument
    """ Catch up on everything variants_classification_changed skipped during the import. CommonVariantClassified
        already records what's been handled, so one sweep per filter covers it - no need to track which
        alleles changed """
    for cgcfv in get_common_filters_in_use():
        _launch_common_variant_classified_tasks(cgcfv, get_classified_high_frequency_variants_qs(cgcfv))
