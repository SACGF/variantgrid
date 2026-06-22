import logging
from typing import Optional

import celery
from django.db.models.expressions import F, OuterRef, Subquery

from library.log_utils import log_traceback
from library.utils.database_utils import run_sql
from patients.models_enums import Zygosity
from snpdb.archive import DataArchivedError
from snpdb.models import Cohort, CohortGenotypeCollection, CohortVersion, SubCohortVariantCollection, \
    VariantCollection, ProcessingStatus


def enqueue_sub_cohort_any_sample_called_vc(cohort: Cohort, run_async: bool = True) -> Optional[str]:
    """ Sole entry point for triggers to schedule the any-sample-called VariantCollection build.
        Returns the celery task id, or None if the cohort is not a sub-cohort or has no samples.
        The task itself is idempotent (drops any pre-existing VC for this cohort version). """
    if not cohort.is_sub_cohort:
        return None
    if not cohort.cohortsample_set.exists():
        return None
    task = build_sub_cohort_any_sample_called_vc_task.si(cohort.pk)
    result = task.apply_async() if run_async else task.apply()
    return result.id


@celery.shared_task(ignore_result=False)
def build_sub_cohort_any_sample_called_vc_task(cohort_id):
    """ Pre-compute the set of variants where at least one sub-cohort sample is called (non-missing), storing
        it as a VariantCollection so the analysis EXCLUDE filter becomes a hash join instead of a
        regex seq-scan against samples_zygosity. Mirrors the build SQL shape in
        analysis.management.commands.profile_analysis_nodes (cohort_exclude_vc_join). Safe to run
        repeatedly - drops any prior collection for this cohort version first. @see issue #1551 """
    try:
        cohort = Cohort.objects.get(pk=cohort_id)
    except Cohort.DoesNotExist:
        return

    if not cohort.is_sub_cohort:
        logging.info("build_sub_cohort_any_sample_called_vc_task: %s is not a sub-cohort, skipping", cohort_id)
        return

    try:
        cgc = cohort.cohort_genotype_collection
    except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        logging.info("build_sub_cohort_any_sample_called_vc_task: %s has no usable CohortGenotypeCollection, skipping",
                     cohort_id)
        return

    cohort_version, _ = CohortVersion.objects.get_or_create(cohort=cohort, version=cohort.version)

    # Drop any pre-existing collection for this version (idempotent re-build). post_delete drops the partition.
    SubCohortVariantCollection.objects.filter(cohort_version=cohort_version).delete()

    # Variants kept by the EXCLUDE filter = those where NOT every sub-cohort sample is missing. This is
    # the same regex CohortGenotypeCollection.get_zygosity_q(exclude=True) uses at filter time.
    missing = [Zygosity.UNKNOWN_ZYGOSITY, Zygosity.MISSING]
    sample_zygosities_dict = {s: missing for s in cohort.get_samples()}
    regex_string = cgc.get_sample_zygosity_regex(sample_zygosities_dict, {})
    regex_excl = f"^((?!{regex_string}))"

    partitions = [cgc.get_partition_table()]
    collection_ids = [cgc.pk]
    if cgc.common_collection_id:
        partitions.append(cgc.common_collection.get_partition_table())
        collection_ids.append(cgc.common_collection_id)

    vc = VariantCollection.objects.create(name=f"sub_cohort_{cohort.pk}_v{cohort.version}_any_sample_called",
                                          status=ProcessingStatus.CREATED)
    vc_partition_table = vc.get_partition_table()
    union_sql_parts = [
        f'SELECT DISTINCT {vc.pk} AS variant_collection_id, variant_id '
        f'FROM "{partition}" '
        f'WHERE collection_id = ANY(%s) AND samples_zygosity ~ %s'
        for partition in partitions
    ]
    build_sql = (f'INSERT INTO "{vc_partition_table}" (variant_collection_id, variant_id) '
                 f'{" UNION ".join(union_sql_parts)};')
    # Each partition predicate needs its own (collection_id, regex) pair, matching its placeholders.
    params = []
    for collection_id in collection_ids:
        params.extend([[collection_id], regex_excl])

    try:
        _, rowcount = run_sql(build_sql, params)
        vc.count = rowcount
        vc.status = ProcessingStatus.SUCCESS
        vc.save()
        SubCohortVariantCollection.objects.create(cohort_version=cohort_version,
                                                  parent_cohort_genotype_collection=cgc,
                                                  variant_collection=vc)
        logging.info("build_sub_cohort_any_sample_called_vc_task: cohort %s built VC %s with %d variants",
                     cohort_id, vc.pk, rowcount)
    except Exception:
        vc.status = ProcessingStatus.ERROR
        vc.save()
        try:
            vc.delete_related_objects()
        except Exception:
            pass
        vc.delete()
        log_traceback()
        raise

    # Opportunistic cleanup of superseded version rows (cheap - small per cohort)
    delete_old_cohort_versions(cohort.pk)


@celery.shared_task(ignore_result=False)
def delete_old_cohort_versions(cohort_id=None):
    """ Drops every CohortVersion that isn't the current version of its Cohort. CASCADE removes any
        version-keyed caches FK'd to it (currently SubCohortVariantCollection, which drops its VC
        partition via post_delete). Mirrors analysis delete_analysis_old_node_versions. """
    qs = CohortVersion.objects.all()
    if cohort_id is not None:
        qs = qs.filter(cohort_id=cohort_id)
    latest = CohortVersion.objects.filter(cohort_id=OuterRef('cohort_id')).order_by('-version')
    sub_query = Subquery(latest.values('pk')[:1])
    qs.annotate(latest_pk=sub_query).exclude(pk=F('latest_pk')).delete()
