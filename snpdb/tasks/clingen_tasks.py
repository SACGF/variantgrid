import celery

from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models import VCF


@celery.shared_task
def populate_clingen_alleles_from_vcf(vcf_id, max_variants=0):
    vcf = VCF.objects.get(pk=vcf_id)
    variants_qs = vcf.get_variant_qs()
    if max_variants:
        # If we're going to limit it, make sure we only get useful ones....
        variants_qs = variants_qs.filter(variantallele__allele__clingen_allele__isnull=True)[:max_variants]
    populate_clingen_alleles_for_variants(vcf.genome_build, variants_qs)
