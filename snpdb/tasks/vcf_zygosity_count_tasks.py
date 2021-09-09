import celery

from snpdb.models import VariantZygosityCountCollection, VCF
from snpdb.variant_zygosity_count import update_variant_zygosity_count_for_vcf


@celery.shared_task
def update_variant_zygosity_count_for_vcf_task(vzcc_id, vcf_id, operation):
    """ Fired when user clicks to manually change zygosity count """
    vcf = VCF.objects.get(pk=vcf_id)
    collection = VariantZygosityCountCollection.objects.get(pk=vzcc_id)
    update_variant_zygosity_count_for_vcf(collection, vcf, operation, manual_override=True)
