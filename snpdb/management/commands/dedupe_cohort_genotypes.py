"""
Removes duplicate CohortGenotype rows - more than one row for a variant in a CohortGenotypeCollection - left by VCF
imports from before per-batch duplicate removal, where nearby multi-allelic records normalised to the same variant
(#2443). Keeps the lowest id, deletes the rest, then recounts the zygosity counts of the affected variants: the count
UPDATE joined to only one of the rows, so what was added can't be subtracted.

Each collection's partition is searched and fixed on its own, in its own transaction. Groups whose rows hold
different genotypes are listed, as keeping the lowest id there drops a call. --dry-run reports without writing.
Registered as a deploy-time ManualOperation (snpdb/migrations/0005_dedupe_cohort_genotypes.py).
"""
from dataclasses import dataclass

from django.core.management.base import BaseCommand
from django.db import connection, transaction

from library.django_utils.database_utils import run_sql
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.models.models_variant import Variant
from snpdb.variant_zygosity_count import recount_variant_zygosity_counts

GENOTYPE_COLUMNS = ["ref_count", "het_count", "hom_count", "unk_count", "samples_zygosity", "samples_allele_depth",
                    "samples_allele_frequency", "samples_read_depth", "samples_genotype_quality",
                    "samples_phred_likelihood", "samples_filters"]


@dataclass(frozen=True)
class DuplicateGenotypes:
    variant_id: int
    ids: list[int]  # ascending, so ids[0] is kept
    genotypes_differ: bool

    @property
    def delete_ids(self) -> list[int]:
        return self.ids[1:]


def find_duplicate_genotypes(cgc: CohortGenotypeCollection) -> list[DuplicateGenotypes]:
    genotype_row = ", ".join(GENOTYPE_COLUMNS)
    sql = f"""
        SELECT variant_id, array_agg(id ORDER BY id), count(DISTINCT ({genotype_row})) > 1
        FROM "{cgc.get_partition_table()}"
        GROUP BY variant_id
        HAVING count(*) > 1
        ORDER BY variant_id
    """
    with connection.cursor() as cursor:
        cursor.execute(sql)
        return [DuplicateGenotypes(*row) for row in cursor.fetchall()]


class Command(BaseCommand):
    help = "Delete duplicate CohortGenotype rows (keeping the lowest id) and recount their zygosity counts"
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would be deleted, change nothing")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        existing_tables = set(connection.introspection.table_names())
        cgc_qs = CohortGenotypeCollection.objects.filter(data_archived_date__isnull=True,
                                                         marked_for_deletion=False).order_by("pk")
        total_variants = total_rows = total_differ = 0
        for cgc in cgc_qs.select_related("cohort"):
            if cgc.get_partition_table() not in existing_tables:
                continue
            with transaction.atomic():
                duplicates = find_duplicate_genotypes(cgc)
                if not duplicates:
                    continue
                delete_ids = [pk for d in duplicates for pk in d.delete_ids]
                differ = [d for d in duplicates if d.genotypes_differ]
                if not dry_run:
                    run_sql(f'DELETE FROM "{cgc.get_partition_table()}" WHERE id = ANY(%s)', [delete_ids])
                    recount_variant_zygosity_counts([d.variant_id for d in duplicates])

            total_variants += len(duplicates)
            total_rows += len(delete_ids)
            total_differ += len(differ)
            verb = "would delete" if dry_run else "deleted"
            self.stdout.write(f"CohortGenotypeCollection {cgc.pk} ({cgc.cohort.name}): {len(duplicates)} duplicated "
                              f"variants, {verb} {len(delete_ids)} rows, {len(differ)} with differing genotypes")
            variants = Variant.objects.in_bulk([d.variant_id for d in differ])
            for d in differ:
                self.stdout.write(f"  differing: variant {d.variant_id} ({variants[d.variant_id]}) "
                                  f"keep {d.ids[0]}, {verb} {d.delete_ids}")

        mode = " (dry run)" if dry_run else ""
        self.stdout.write(f"Total{mode}: {total_variants} duplicated variants, {total_rows} rows, "
                          f"{total_differ} with differing genotypes")
