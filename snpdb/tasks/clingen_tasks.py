import celery

from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models import AlleleSource


@celery.task
def populate_clingen_alleles_from_allele_source(allele_source_id, max_variants=0):
    allele_source = AlleleSource.objects.get_subclass(pk=allele_source_id)
    variants_qs = allele_source.get_variant_qs()
    if max_variants:
        # If we're going to limit it, make sure we only get useful ones....
        variants_qs = variants_qs.filter(variantallele__allele__clingen_allele__isnull=True)[:max_variants]
    populate_clingen_alleles_for_variants(allele_source.get_genome_build(), variants_qs)
