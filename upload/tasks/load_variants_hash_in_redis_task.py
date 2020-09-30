import celery

from snpdb.variant_pk_lookup import VariantPKLookup


@celery.task
def load_variants_hash_in_redis():
    """ Make this a task so that we can run it in a queue and guarantee only 1 copy is doing this """

    variant_pk_lookup = VariantPKLookup(redis_check=False)
    variant_pk_lookup.update_redis_hashes()
