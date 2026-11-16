"""
Unbinds VCF sample fields (genotype_field, read_depth_field...) bound to a FORMAT key the header declares but no record
uses, for VCFs imported before the import checked the records (#2833) - a DRAGEN TSO500 CNV VCF declares GT, AD, DP and
GQ but every record is only FC, so its calls were counted as nothing and the Sample node offered zygosity and depth
filters. Candidates come from sample stats alone: a GT field where every sample has no ref/het/hom call and some unknown
ones means no record had a GT. Only those VCFs (small CNV/SV files) have their genotype rows read, for which bound keys
the stored CohortGenotype.format JSON holds, so this never scans a large partition. A VCF whose records have GT but not
a declared depth field keeps it until reloaded; one without sample stats is skipped. Sample stats are recalculated for
a VCF that changed, so its contents show.
--dry-run reports without writing.
Registered as a deploy-time ManualOperation (upload/migrations/0005_one_off_clear_unused_vcf_format_fields.py).
"""
from django.core.management.base import BaseCommand
from django.db import connection
from django.db.models import F, Sum

from annotation.models.models_annotation_version import AnnotationVersion
from annotation.tasks.calculate_sample_stats import calculate_vcf_stats
from snpdb.enums import ImportStatus
from snpdb.models import VCF, VCFSourceSettings
from snpdb.models.models_cohort_stats import CohortGenotypeStats
from upload.vcf.vcf_import import (
    clear_format_fields_not_in_records,
    get_format_fields_not_in_records,
)


def get_no_genotype_call_vcf_ids() -> set[int]:
    """ VCFs with a GT field whose own cohort's per-sample stats show no ref/het/hom call at all """
    stats_qs = CohortGenotypeStats.objects.filter(
        sample__isnull=False, filter_key__isnull=True, passing_filter=False,
        sample__vcf__genotype_field__isnull=False,
        cohort_genotype_collection__cohort__vcf_id=F("sample__vcf_id"))
    vcf_calls_qs = stats_qs.values("sample__vcf_id").annotate(
        calls=Sum("ref_count") + Sum("het_count") + Sum("hom_count"), unknown=Sum("unk_count"))
    return set(vcf_calls_qs.filter(calls=0, unknown__gt=0).values_list("sample__vcf_id", flat=True))


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
        vcf_qs = VCF.objects.filter(pk__in=get_no_genotype_call_vcf_ids(), genotype_samples__gt=0, cohort__isnull=False,
                                    data_archived_date__isnull=True, import_status=ImportStatus.SUCCESS)
        num_changed = 0
        for vcf in vcf_qs.select_related("genome_build").order_by("pk"):
            tables = [t for t in get_genotype_tables(vcf) if t in existing_tables]
            if not (tables and _exists(tables, "jsonb_array_length(format) > 0")):
                continue  # Imported before FORMAT was stored - nothing to tell used from unused

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
        self.stdout.write(f"{num_changed} VCFs{mode} with a GT field but no genotype calls had unused FORMAT fields")
