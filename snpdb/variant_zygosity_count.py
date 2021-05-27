#
from typing import Tuple, List

from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Max
from django.utils import timezone
import logging

from eventlog.models import create_event
from library.database_utils import run_sql
from library.enums.log_level import LogLevel
from library.log_utils import get_traceback
from snpdb.models import VCF, VariantZygosityCountForVCF, VariantZygosityCountForSample, Sample, Variant, \
    VariantZygosityCount, VariantZygosityCountCollection


def create_variant_zygosity_counts():
    """ This needs to be run in a Queue w/1 thread so there are no race conditions """

    data = Variant.objects.filter(alt__isnull=False).aggregate(max_variant=Max("id"))
    max_variant_id = data["max_variant"] or 0

    for collection in VariantZygosityCountCollection.objects.all():
        data = VariantZygosityCount.objects.filter(collection=collection).aggregate(max_variant=Max("variant_id"))
        max_vzc_variant_id = data["max_variant"] or 0
        if max_variant_id > max_vzc_variant_id:
            logging.info("Current variant_id=%d, Need to create for collection %s up to variant_id=%d",
                         max_vzc_variant_id, collection, max_variant_id)

            params = {"table_name": collection.get_partition_table()}
            sql = """INSERT INTO %(table_name)s
                    (variant_id, collection_id, ref_count, het_count, hom_count, unk_count)\n""" % params
            sql += """select id, %s, 0, 0, 0, 0 from snpdb_variant where id > %s;"""
            run_sql(sql, [collection.pk, max_vzc_variant_id])


def check_valid_count_ops(operation):
    VALID_OPERATIONS = ['+', '-']
    if operation not in VALID_OPERATIONS:
        msg = f"Operation: {operation} not valid, must be one of {VALID_OPERATIONS}"
        raise ValueError(msg)


def _get_update_sql_and_params(collection: VariantZygosityCountCollection, vcf: VCF, operation, sample=None) -> Tuple[str, List]:
    extra_where_sql = ""
    cgc = vcf.cohort.cohort_genotype_collection
    if sample:
        sample_index = cgc.get_sql_index_for_sample_id(sample.pk)
        sample_zygosity = 'SUBSTRING("snpdb_cohortgenotype"."samples_zygosity", %d, 1)' % sample_index
        ref_count_delta = f"(case when {sample_zygosity} = 'R' then 1 else 0 end)"
        het_count_delta = f"(case when {sample_zygosity} = 'E' then 1 else 0 end)"
        hom_count_delta = f"(case when {sample_zygosity} = 'O' then 1 else 0 end)"
        extra_where_sql = f" AND {sample_zygosity} in ('R', 'E', 'O')"
    else:
        ref_count_delta = "snpdb_cohortgenotype.ref_count"
        het_count_delta = "snpdb_cohortgenotype.het_count"
        hom_count_delta = "snpdb_cohortgenotype.hom_count"

    update_params = {"cohort_genotype_collection_id": cgc.pk,
                     "table_name": collection.get_partition_table(),
                     "operation": operation,
                     "ref_count_delta": ref_count_delta,
                     "het_count_delta": het_count_delta,
                     "hom_count_delta": hom_count_delta}
    update_sql = """
        UPDATE %(table_name)s
        SET
            ref_count = %(table_name)s.ref_count %(operation)s %(ref_count_delta)s,
            het_count = %(table_name)s.het_count %(operation)s %(het_count_delta)s,
            hom_count = %(table_name)s.hom_count %(operation)s %(hom_count_delta)s """

    # This knows to join to the snpdb_cohortgenotype partition via where clause
    from_sql = """
        FROM snpdb_cohortgenotype
        WHERE
        %(table_name)s.variant_id = snpdb_cohortgenotype.variant_id
        AND snpdb_cohortgenotype.collection_id = %(cohort_genotype_collection_id)d
    """
    sql = (update_sql + from_sql + extra_where_sql) % update_params
    return sql, []


def update_all_variant_zygosity_counts_for_vcf(vcf: VCF, operation):
    """ This should run only once - as it could cause deadlocks if running simultaneously """

    for collection in VariantZygosityCountCollection.objects.all():
        update_variant_zygosity_count_for_vcf(collection, vcf, operation)


def update_variant_zygosity_count_for_vcf(collection: VariantZygosityCountCollection, vcf: VCF, operation,
                                          manual_override=False):
    UPDATE_VARIANT_ZYG_COUNT_EVENT = 'update_variant_zygosity_count_for_vcf'

    try:
        check_valid_count_ops(operation)
        logging.info("update_variant_zygosity_count_for_vcf(%s, %s, %s)", collection, vcf, operation)
        if not (vcf.has_genotype and (vcf.variant_zygosity_count or manual_override)):
            logging.info("VCF %s genotype=%s, variant_zygosity_count=%s - skipping VariantZygosityCount",
                         vcf, vcf.has_genotype, vcf.variant_zygosity_count)
            return

        try:
            cohort_genotype_collection = vcf.cohort.cohort_genotype_collection
        except ObjectDoesNotExist:
            logging.info("VCF %s has no cohort genotype collection", vcf)
            return

        # Want to be able to throw exceptions here
        if operation == '+':
            vzcv, created = VariantZygosityCountForVCF.objects.get_or_create(collection=collection, vcf=vcf)
            if created:
                if vzcv.deleted is None:
                    logging.warning("VCF pk=%d, collection=%s (Add) existing non-deleted VariantZygosityCountForVCF. Skipping",
                                    vcf.pk, collection.name)
                    return
                else:
                    vzcv.count_complete = None
                    vzcv.deleted = None
                    vzcv.save()
        else:
            try:
                vzcv = VariantZygosityCountForVCF.objects.get(collection=collection, vcf=vcf)
            except VariantZygosityCountForVCF.DoesNotExist:
                logging.warning("No VariantZygosityCountForVCF for VCF pk=%d (%s)",
                                vcf.pk, vcf.get_import_status_display())
                return  # no need to do anything

            vzcv.check_can_delete()

        use_cohort_genotype_collection = vzcv and not vzcv.is_split_to_sample_counts
        if use_cohort_genotype_collection:
            logging.info("Updating from Cohort Zygosity Collection %d", cohort_genotype_collection.pk)
            sql, params = _get_update_sql_and_params(collection, vcf, operation)
            run_sql(sql, params)

            now = timezone.now()
            if operation == '+':
                vzcv.count_complete = now
            else:
                vzcv.deleted = now
            vzcv.save()
        else:
            logging.info("Updating from sample...")

            for sample in vcf.sample_set.all():
                update_variant_zygosity_count_for_sample(collection, sample, operation)

    except:
        tb = get_traceback()
        details = f"vcf_id={vcf}, operation='{operation}', traceback: {tb}"
        create_event(None, UPDATE_VARIANT_ZYG_COUNT_EVENT, details=details, severity=LogLevel.ERROR)
        logging.error(details)
        raise


def update_all_variant_zygosity_counts_for_sample(sample, operation):
    for collection in VariantZygosityCountCollection.objects.all():
        update_variant_zygosity_count_for_sample(collection, sample, operation)


def update_variant_zygosity_count_for_sample(collection: VariantZygosityCountCollection, sample: Sample, operation):
    try:
        check_valid_count_ops(operation)
        logging.info("update_variant_zygosity_count_for_sample(%d, %s)", sample, operation)
        try:
            gvzcp = VariantZygosityCountForVCF.objects.get(vcf=sample.vcf)  # throws DoesNotExist if not there
            if operation == '+':
                msg = f"Attempting to count for sample ({sample}) which has already had vcf ({sample.vcf}) counted!"
                raise ValueError(msg)
            gvzcp.check_can_delete()
            if not gvzcp.is_split_to_sample_counts:
                gvzcp.split_to_sample_counts()

        except VariantZygosityCountForVCF.DoesNotExist:
            pass

        if operation == '+':
            gvzcs = VariantZygosityCountForSample.objects.create(collection=collection, sample=sample)
        else:
            gvzcs = VariantZygosityCountForSample.objects.get(collection=collection, sample=sample)
            gvzcs.check_can_delete()
            if gvzcs.count_complete is None:
                msg = f"VariantZygosityCountForSample ({sample}) count_complete was None (never completed) - aborting"
                raise ValueError(msg)

        logging.info("Updating from single sample %s...", sample)
        sql, params = _get_update_sql_and_params(collection, sample.vcf, operation, sample=sample)
        run_sql(sql, params)

        now = timezone.now()
        if operation == '+':
            gvzcs.count_complete = now
        else:
            gvzcs.deleted = now
        gvzcs.save()

    except:
        tb = get_traceback()
        details = f"sample_id={sample}, operation='{operation}', traceback: {tb}"
        create_event(None, 'update_variant_zygosity_count_for_sample', details=details, severity=LogLevel.ERROR)
        logging.error(details)
        raise
