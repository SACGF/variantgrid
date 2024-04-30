import logging

from django.dispatch import receiver

from classification.models import ImportedAlleleInfo
from library.log_utils import report_event
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models import liftover_run_complete_signal, LiftoverRun, AlleleConversionTool, Variant


@receiver(liftover_run_complete_signal, sender=LiftoverRun)
def liftover_run_complete_handler(sender, instance: LiftoverRun, **kwargs):

    logging.info("Classification liftover_run_complete_handler")

    allele_qs = instance.get_allele_qs()
    if instance.conversion_tool != AlleleConversionTool.CLINGEN_ALLELE_REGISTRY:
        build_variants = Variant.objects.filter(variantallele__genome_build=instance.genome_build,
                                                variantallele__allele__in=allele_qs)
        populate_clingen_alleles_for_variants(instance.genome_build, build_variants)

    ImportedAlleleInfo.relink_variants(liftover_run=instance)
    report_event('Completed import liftover',
                 extra_data={'liftover_id': instance.pk, 'allele_count': allele_qs.count()})
