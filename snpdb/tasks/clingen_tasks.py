import celery

from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models import AlleleSource


@celery.task
def populate_clingen_alleles_from_allele_source(allele_source_id):
    allele_source = AlleleSource.objects.get_subclass(pk=allele_source_id)
    populate_clingen_alleles_for_variants(allele_source.get_genome_build(), allele_source.get_variant_qs())
