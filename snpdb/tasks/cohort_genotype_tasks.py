import logging
import time
from collections import defaultdict
from typing import Optional

import celery
from celery.result import AsyncResult
from django.db.models.query_utils import Q

from library.django_utils.django_postgres import pg_sql_array, model_to_insert_sql
from library.log_utils import log_traceback
from library.utils import single_quote
from library.utils.database_utils import run_sql
from patients.models_enums import Zygosity
from snpdb.common_variants import get_common_filter
from snpdb.models import Cohort, ImportStatus, CohortGenotypeCommonFilterVersion, Variant, CommonVariantClassified, \
    CohortGenotypeCollection, CohortGenotype, CohortGenotypeTaskVersion, CohortGenotypeCollectionType


def create_cohort_genotype_and_launch_task(cohort, run_async=True):
    logging.info("create_cohort_genotype_and_launch_task")
    cohort.delete_old_counts()

    q_not_sub_cohort = Q(parent_cohort__isnull=True)
    q_is_vcf = Q(vcf__isnull=False)  # Don't build off a cohort that could get deleted
    q_exclude_this_cohort = ~Q(pk=cohort.pk)
    q = q_not_sub_cohort & q_is_vcf & q_exclude_this_cohort

    containing_cohort = Cohort.get_cohort_containing_all_samples(cohort.get_sample_ids(), extra_q=q)
    status = None
    celery_task = None
    if containing_cohort:
        logging.info("Cohort %s can be fully contained within %s!", cohort, containing_cohort)
        cohort.parent_cohort = containing_cohort
        cohort.import_status = containing_cohort.import_status
        cohort.save()
        status = "SUCCESS"
    else:
        launch_task = False
        try:
            cgc = cohort.cohort_genotype_collection
            celery_task = cgc.celery_task
            logging.warning("This count task was already running")
            result = AsyncResult(celery_task)
            if (not result.successful()) or cohort.import_status == ImportStatus.ERROR:
                launch_task = True

                logging.info("Deleting existing and creating new partition")
                cgc.delete_related_objects()
                cgc.delete()
                cgc = create_cohort_genotype_collection(cohort)
        except CohortGenotypeCollection.DoesNotExist:
            cgc = create_cohort_genotype_collection(cohort)
            launch_task = True

        if launch_task:
            task = cohort_genotype_task.si(cgc.id)  # @UndefinedVariable
            if run_async:
                result = task.apply_async()
            else:
                result = task.apply()

            # Do as an update as other jobs may be modifying objects
            celery_task = result.id
            CohortGenotypeCollection.objects.filter(cohort_id=cohort.pk).update(celery_task=celery_task)
            Cohort.objects.filter(pk=cohort.pk).update(import_status=ImportStatus.IMPORTING)

    return status, celery_task


def create_cohort_genotype_collection(cohort):
    if not cohort.cohortsample_set.exists():
        msg = f"create_cohort_genotype_collection called for empty cohort {cohort.pk}: {cohort}"
        raise ValueError(msg)

    name = f"{cohort.name} ({cohort.pk}:{cohort.version})"
    num_samples = cohort.cohortsample_set.count()
    common_collection = None
    kwargs = {
        "cohort": cohort,
        "cohort_version": cohort.version,
        "num_samples": num_samples,
    }

    if common_filter := get_common_filter(cohort.genome_build):
        common_collection = CohortGenotypeCollection.objects.create(name=f"{name} common",
                                                                    common_filter=common_filter,
                                                                    collection_type=CohortGenotypeCollectionType.COMMON,
                                                                    **kwargs)
        logging.info(f"Created common collection: {common_collection}")

    collection = CohortGenotypeCollection.objects.create(name=name,
                                                         common_collection=common_collection,
                                                         collection_type=CohortGenotypeCollectionType.UNCOMMON,
                                                         **kwargs)

    logging.info(f"Created {collection}")
    return collection


def _get_sample_zygosity_count_sql(sample_value, zygosity):
    return f'CASE WHEN (({sample_value}) = \'{zygosity}\') THEN 1 ELSE 0 END'


def _get_left_outer_join_on_variant(partition_table):
    return f'LEFT OUTER JOIN "{partition_table}" ON ("snpdb_variant"."id" = "{partition_table}"."variant_id")'


def _get_insert_cohort_genotype_sql(cgc: CohortGenotypeCollection, common=False) -> Optional[str]:
    """ common: use common cohort genotype collection
        For old VCFs, imported before we had common filters - this can return None
        as nothing to do (everything in rare/default)
     """
    cohort = cgc.cohort
    samples = cohort.get_samples()

    zygosity_count_lists = defaultdict(list)
    column_pack_lists = defaultdict(list)
    joins = set()

    ZYGOSITY_COUNT_COLUMNS = {
        Zygosity.HOM_REF: "ref_count",
        Zygosity.HET: "het_count",
        Zygosity.HOM_ALT: "hom_count",
        Zygosity.UNKNOWN_ZYGOSITY: "unk_count",
    }

    for sample in samples:
        sample_cgc = sample.vcf.cohort.cohort_genotype_collection
        if common:
            if cc := sample_cgc.common_collection:
                sample_cgc = cc
            else:
                # No common - all in normal
                continue
        partition_table = sample_cgc.get_partition_table()
        joins.add(_get_left_outer_join_on_variant(partition_table))
        i = sample_cgc.get_sql_index_for_sample_id(sample.pk)

        for column, (is_array, empty_value) in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE.items():
            sample_field = f"{partition_table}.{column}"
            if is_array:
                sample_field = f"{sample_field}[{i}]"
            else:
                sample_field = f"SUBSTRING({sample_field}, {i}, 1)"
                empty_value = single_quote(empty_value)

            sample_value = f"coalesce({sample_field}, {empty_value})"
            column_pack_lists[column].append(sample_value)

            if column == "samples_zygosity":
                for zygosity in ZYGOSITY_COUNT_COLUMNS:
                    zygosity_count_lists[zygosity].append(_get_sample_zygosity_count_sql(sample_value, zygosity))

    empty_json_list = "'[]'::jsonb"
    columns = {
        "variant_id": '"snpdb_variant"."id"',
        "collection_id": f"{cgc.pk}",
        # TODO: These should be joined properly, only doing this to avoid not NULL constraint
        "format": empty_json_list,
        "info": empty_json_list,
    }

    for zygosity, c in ZYGOSITY_COUNT_COLUMNS.items():
        columns[c] = " + ".join(zygosity_count_lists[zygosity])

    for c, (is_array, _) in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE.items():
        if is_array:
            columns[c] = pg_sql_array(column_pack_lists[c])
        else:
            columns[c] = " || ".join(column_pack_lists[c])

    if not joins:
        return None

    joins = "\n".join(joins)
    insert_sql = f"""
insert into {cgc.get_partition_table()}
({','.join(columns)})
SELECT DISTINCT
{','.join(columns.values())}
FROM snpdb_variant
{joins}
WHERE
({columns['ref_count']} + {columns['het_count']} + {columns['hom_count']} + {columns['unk_count']}) > 0
"""
    # For the where clause, we could simplify (if running out of max SQL limit)
    # CASE WHEN COALESCE(SUBSTRING(samples_zygosity, 2, 1), '.') IN ('R', 'E', 'O', 'U') THEN 1 ELSE 0 END AS match_2
    # Though we are already calculating the columns in the select, so probably best to leave as it
    # Doesn't seem to have performance hit

    return insert_sql


@celery.shared_task(ignore_result=False)
def cohort_genotype_task(cohort_genotype_collection_id):
    """ This takes a cohort, performs a count on it, stores it to the database,
        and then saves the stored count to the cohort object  """

    cgc = CohortGenotypeCollection.objects.get(pk=cohort_genotype_collection_id)

    # CohortGenotypeTaskVersion versions:
    # 1. Legacy - inserted to tag legacy data before we added task_version
    # 2. cohort_genotype_task - 20200319 - Split hom_count into ref_count, hom_count
    # 3. cohort_genotype_task - 20200514 - Use code from grid_sample_columns
    # 4. cohort_genotype_task - 20241125 - Handle rare/common
    NAME = "cohort_genotype_task - 20241125 - packed fields (handle rare/common)"  # Change this if the data changes!
    task_version, _ = CohortGenotypeTaskVersion.objects.get_or_create(name=NAME)

    cohort = cgc.cohort
    cohort.import_status = ImportStatus.IMPORTING
    cohort.save()

    try:
        cohort_genotype_collection_list = [
            (cgc, False)
        ]
        if cc := cgc.common_collection:
            cohort_genotype_collection_list.append((cc, True))

        for cgc, common in cohort_genotype_collection_list:
            if insert_sql := _get_insert_cohort_genotype_sql(cgc, common=common):
                logging.info("cohort_genotype_task %s/common=%s SQL:", cgc, common)
                logging.info(insert_sql)

                start = time.time()
                run_sql(insert_sql)
                end = time.time()
                logging.info("SQL took %d secs to run", end - start)
                cgc.task_version = task_version
                cgc.save()

        # Do as an update as other jobs may be modifying object
        import_status = ImportStatus.SUCCESS
        Cohort.objects.filter(pk=cohort.pk).update(import_status=import_status)
    except:
        import_status = ImportStatus.ERROR
        Cohort.objects.filter(pk=cohort.pk).update(import_status=import_status)
        log_traceback()
        raise  # So it errors


@celery.shared_task(ignore_result=False)
def common_variant_classified_task(variant_id, common_filter_id):
    try:
        variant = Variant.objects.get(pk=variant_id)
        common_filter = CohortGenotypeCommonFilterVersion.objects.get(pk=common_filter_id)

        # We should only be called if CommonVariantClassified doesn't exist
        logging.info("common_variant_classified_task(%s, %s)", str(variant), str(common_filter))

        cohort_genotype_delete_ids = []
        for cgc in common_filter.cohortgenotypecollection_set.all():
            uncommon = cgc.uncommon
            for cg in cgc.cohortgenotype_set.filter(variant=variant):
                # Because they are in different partitions we can't just update them, need to re-insert and then delete
                # old ones
                cohort_genotype_delete_ids.append(cg.pk)
                cg.pk = None
                cg.collection = uncommon
                insert_sql = model_to_insert_sql([cg], ignore_fields=['id'],
                                                 db_table=uncommon.get_partition_table())[0]
                run_sql(insert_sql)

        # No errors, insert must have gone through ok
        logging.info("Deleting %d cohort genotype records", len(cohort_genotype_delete_ids))
        CohortGenotype.objects.filter(pk__in=cohort_genotype_delete_ids).delete()

        # When it's completed - we can create the record
        CommonVariantClassified.objects.create(variant=variant, common_filter=common_filter)
    except:
        log_traceback()
        raise
