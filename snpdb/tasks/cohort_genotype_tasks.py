import logging
import time
from collections import defaultdict

import celery
from celery.result import AsyncResult
from django.db.models.query_utils import Q

from library.database_utils import run_sql
from library.log_utils import log_traceback
from library.utils import single_quote
from patients.models_enums import Zygosity
from snpdb.grid_columns.grid_sample_columns import get_left_outer_join_on_variant
from snpdb.models import Cohort, ImportStatus
from snpdb.models import CohortGenotypeCollection, CohortGenotype, CohortGenotypeTaskVersion


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
            cgc = CohortGenotypeCollection.objects.get(cohort=cohort, cohort_version=cohort.version)
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
    common_collection = CohortGenotypeCollection.objects.create(name=f"{name} common",
                                                                cohort=cohort,
                                                                cohort_version=0,  # so it isn't retrieved
                                                                num_samples=cohort.cohortsample_set.count())

    collection = CohortGenotypeCollection.objects.create(name=f"{name} rare",
                                                         cohort=cohort,
                                                         cohort_version=cohort.version,
                                                         common_collection=common_collection,
                                                         num_samples=cohort.cohortsample_set.count())

    logging.info(f"Created {collection}")
    return collection


def _get_sample_zygosity_count_sql(sample_value, zygosity):
    return f'CASE WHEN (({sample_value}) = \'{zygosity}\') THEN 1 ELSE 0 END'


def pg_sql_array(values):
    return 'array[%s]' % ','.join(values)


def get_insert_cohort_genotype_sql(cgc: CohortGenotypeCollection):
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
        partition_table = sample_cgc.get_partition_table()
        joins.add(get_left_outer_join_on_variant(partition_table))
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

    columns = {
        "variant_id": '"snpdb_variant"."id"',
        "collection_id": f"{cgc.pk}",
    }

    for zygosity, c in ZYGOSITY_COUNT_COLUMNS.items():
        columns[c] = " + ".join(zygosity_count_lists[zygosity])

    for c, (is_array, _) in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE.items():
        if is_array:
            columns[c] = pg_sql_array(column_pack_lists[c])
        else:
            columns[c] = " || ".join(column_pack_lists[c])

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
    NAME = "cohort_genotype_task - 20200514 - packed fields"  # Change this if the data changes!
    task_version, _ = CohortGenotypeTaskVersion.objects.get_or_create(name=NAME)

    cohort = cgc.cohort
    cohort.import_status = ImportStatus.IMPORTING
    cohort.save()

    try:
        insert_sql = get_insert_cohort_genotype_sql(cgc)
        logging.info("cohort_genotype_task SQL:")
        logging.info(insert_sql)

        start = time.time()
        run_sql(insert_sql)
        end = time.time()
        logging.info("SQL took %d secs to run", end - start)

        # Do as an update as other jobs may be modifying object
        import_status = ImportStatus.SUCCESS
        Cohort.objects.filter(pk=cohort.pk).update(import_status=import_status)

        cgc.task_version = task_version
        cgc.save()
    except:
        import_status = ImportStatus.ERROR
        Cohort.objects.filter(pk=cohort.pk).update(import_status=import_status)
        log_traceback()
        raise  # So it errors
