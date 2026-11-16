"""
Unbinds VCF sample fields (genotype_field, read_depth_field...) bound to a FORMAT key the header declares but no record
uses, for VCFs imported before the import checked the records (#2833) - a DRAGEN TSO500 CNV VCF declares GT, AD, DP and
GQ but every record is only FC, so its calls were counted as nothing and the Sample node offered zygosity and depth
filters. Which keys the records used comes from the stored CohortGenotype.format JSON; a VCF whose rows predate that
column being filled is skipped. Sample stats are recalculated for a VCF that changed, so its contents show.
--dry-run reports without writing.
Registered as a deploy-time ManualOperation (upload/migrations/0005_one_off_clear_unused_vcf_format_fields.py).
"""
from django.core.management.base import BaseCommand
from django.db import connection

from annotation.models.models_annotation_version import AnnotationVersion
from annotation.tasks.calculate_sample_stats import calculate_vcf_stats
from snpdb.enums import ImportStatus
from snpdb.models import VCF, VCFSourceSettings
from upload.vcf.vcf_import import (
    clear_format_fields_not_in_records,
    get_format_fields_not_in_records,
)


def get_genotype_tables(vcf: VCF) -> list[str]:
    cgc = vcf.cohort.cohort_genotype_collection
    tables = [cgc.get_partition_table()]
    if cgc.common_collection_id:
        tables.append(cgc.get_partition_table(base_table_name="common"))
    return tables


def _exists(tables: list[str], where: str, params=None) -> bool:
    with connection.cursor() as cursor:
        for table in tables:
            cursor.execute(f'SELECT EXISTS(SELECT 1 FROM "{table}" WHERE {where})', params)
            if cursor.fetchone()[0]:
                return True
    return False


def get_record_format_fields(vcf: VCF, tables: list[str]) -> set[str]:
    """ Of the FORMAT keys the VCF has a field bound to, the ones some record carries. Every sample of a record shares
        its FORMAT, so the first sample's keys are the record's """
    bound_keys = {getattr(vcf, field) for field in VCFSourceSettings.OVERRIDABLE_SAMPLE_FIELDS} - {None}
    return {key for key in bound_keys if _exists(tables, "format->0 ? %s", [key])}


class Command(BaseCommand):
    help = "Unbind VCF sample fields whose FORMAT key the header declares but no record uses, and recalculate stats"
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would change, change nothing")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        existing_tables = set(connection.introspection.table_names())
        vcf_qs = VCF.objects.filter(genotype_samples__gt=0, cohort__isnull=False, data_archived_date__isnull=True,
                                    import_status=ImportStatus.SUCCESS).select_related("genome_build").order_by("pk")
        num_vcfs = vcf_qs.count()
        num_changed = 0
        num_no_format = 0
        for vcf in vcf_qs:
            tables = [t for t in get_genotype_tables(vcf) if t in existing_tables]
            if not tables:
                continue
            if not _exists(tables, "jsonb_array_length(format) > 0"):
                num_no_format += 1  # Imported before FORMAT was stored - nothing to tell used from unused
                continue

            record_format_fields = get_record_format_fields(vcf, tables)
            if dry_run:
                unbound = get_format_fields_not_in_records(vcf, record_format_fields)
            else:
                unbound = clear_format_fields_not_in_records(vcf, record_format_fields)
            if not unbound:
                continue

            num_changed += 1
            unused = ", ".join(sorted(set(unbound.values())))
            self.stdout.write(f"VCF {vcf.pk} ({vcf.name}): FORMAT {unused} declared but in no record - "
                              f"unbinding {', '.join(sorted(unbound))}")
            if not dry_run:
                annotation_version = AnnotationVersion.latest(vcf.genome_build)
                calculate_vcf_stats(vcf.pk, annotation_version.pk)

        mode = " (dry run)" if dry_run else ""
        self.stdout.write(f"Checked {num_vcfs} VCFs{mode}: {num_changed} with unused FORMAT fields, "
                          f"{num_no_format} skipped as imported before FORMAT was stored")
