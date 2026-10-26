"""
Removes duplicate CohortGenotype rows - more than one row for a variant in a CohortGenotypeCollection - left by VCF
imports from before per-batch duplicate removal, where nearby multi-allelic records normalised to the same variant
(#2443). A row is deleted only when its genotype data is identical to the variant's lowest id row. Zygosity counts need
no recount: the count UPDATE ... FROM snpdb_cohortgenotype applies one joined row per variant, and identical rows add
the same values, so the counts already match a single row.

A row whose genotype differs from the lowest id row is left in place and reported as a warning: one VCF giving the
same variant two genotypes is an import / normalisation bug, and picking one would silently drop a call. Each
collection's partition is searched and fixed on its own, in its own transaction. --dry-run reports without writing.
Output is a summary at the end (big deployments have thousands of affected collections); --csv writes one row per
collection with differing duplicates.
Registered as a deploy-time ManualOperation (snpdb/migrations/0006_one_off_dedupe_cohort_genotypes.py).
"""
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from django.core.management.base import BaseCommand
from django.db import connection, transaction

from library.django_utils.database_utils import run_sql
from library.log_utils import report_message
from snpdb.models.models_cohort import CohortGenotypeCollection

GENOTYPE_COLUMNS = ["ref_count", "het_count", "hom_count", "unk_count", "samples_zygosity", "samples_allele_depth",
                    "samples_allele_frequency", "samples_read_depth", "samples_genotype_quality",
                    "samples_phred_likelihood", "samples_filters"]


@dataclass
class DuplicateGenotypes:
    variant_id: int
    keep_id: int
    identical_ids: list[int] = field(default_factory=list)
    differing_ids: list[int] = field(default_factory=list)
    calls_differ: bool = False  # samples_zygosity differs, not just the evidence (AD, DP, GQ...)


@dataclass
class CollectionDifferences:
    cgc: CohortGenotypeCollection
    variants: int
    calls_differ: int


def find_duplicate_genotypes(cgc: CohortGenotypeCollection) -> list[DuplicateGenotypes]:
    """ Every variant with more than one row, each extra row compared with the variant's lowest id row """
    table = cgc.get_partition_table()
    duplicate_row = ", ".join(f"duplicate.{c}" for c in GENOTYPE_COLUMNS)
    kept_row = ", ".join(f"kept.{c}" for c in GENOTYPE_COLUMNS)
    sql = f"""
        WITH duplicated AS (
            SELECT variant_id, min(id) AS keep_id FROM "{table}" GROUP BY variant_id HAVING count(*) > 1
        )
        SELECT duplicated.variant_id, duplicated.keep_id, duplicate.id,
               ({duplicate_row}) IS NOT DISTINCT FROM ({kept_row}),
               duplicate.samples_zygosity IS NOT DISTINCT FROM kept.samples_zygosity
        FROM duplicated
        JOIN "{table}" kept ON kept.id = duplicated.keep_id
        JOIN "{table}" duplicate ON duplicate.variant_id = duplicated.variant_id AND duplicate.id <> duplicated.keep_id
        ORDER BY duplicated.variant_id, duplicate.id
    """
    duplicates: dict[int, DuplicateGenotypes] = {}
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for variant_id, keep_id, duplicate_id, identical, same_calls in cursor.fetchall():
            if not (duplicate := duplicates.get(variant_id)):
                duplicate = duplicates[variant_id] = DuplicateGenotypes(variant_id, keep_id)
            (duplicate.identical_ids if identical else duplicate.differing_ids).append(duplicate_id)
            duplicate.calls_differ |= not same_calls
    return list(duplicates.values())


class Command(BaseCommand):
    help = ("Delete duplicate CohortGenotype rows identical to the variant's first row, and warn about duplicates "
            "with differing genotypes")
    category = "one-off"

    LARGEST_SHOWN = 20
    PROGRESS_EVERY = 1000

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would be deleted, change nothing")
        parser.add_argument('--csv', help="Write one row per collection with differing duplicates to this file")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        existing_tables = set(connection.introspection.table_names())
        cgc_qs = CohortGenotypeCollection.objects.filter(data_archived_date__isnull=True,
                                                         marked_for_deletion=False).order_by("pk")
        num_collections = cgc_qs.count()
        totals = defaultdict(int)
        differences: list[CollectionDifferences] = []
        for i, cgc in enumerate(cgc_qs.select_related("cohort__vcf"), start=1):
            if i % self.PROGRESS_EVERY == 0:
                self.stdout.write(f"Checked {i}/{num_collections} collections: {totals['deleted']} identical rows, "
                                  f"{totals['differing']} variants with differing genotypes so far")
            if cgc.get_partition_table() not in existing_tables:
                continue
            with transaction.atomic():
                duplicates = find_duplicate_genotypes(cgc)
                if not duplicates:
                    continue
                delete_ids = [pk for d in duplicates for pk in d.identical_ids]
                if delete_ids and not dry_run:
                    run_sql(f'DELETE FROM "{cgc.get_partition_table()}" WHERE id = ANY(%s)', [delete_ids])

            if delete_ids:
                totals["deleted"] += len(delete_ids)
                totals["deleted_collections"] += 1
            if differing := [d for d in duplicates if d.differing_ids]:
                totals["differing"] += len(differing)
                differences.append(CollectionDifferences(cgc, len(differing), sum(d.calls_differ for d in differing)))

        self._write_summary(num_collections, totals, differences, dry_run)
        if options["csv"]:
            self._write_csv(options["csv"], differences)
        if differences and not dry_run:
            report_message(self._differing_message(totals, differences),
                           extra_data={"differing_variants": totals["differing"],
                                       "collections": len(differences)})

    @staticmethod
    def _differing_message(totals, differences: list[CollectionDifferences]) -> str:
        return (f"{totals['differing']} variants in {len(differences)} CohortGenotypeCollections have more than one "
                "genotype row that disagree, left in place: the caller reported the same variant twice with different "
                "evidence (#2443). Reloading the VCF dedupes it, keeping the first record")

    def _write_summary(self, num_collections: int, totals, differences: list[CollectionDifferences], dry_run: bool):
        mode = " (dry run)" if dry_run else ""
        verb = "to delete" if dry_run else "deleted"
        self.stdout.write(f"\nChecked {num_collections} collections{mode}")
        self.stdout.write(f"Identical duplicate rows {verb}: {totals['deleted']} in {totals['deleted_collections']} "
                          "collections")
        if not differences:
            self.stdout.write("No duplicates with differing genotypes")
            return

        calls_differ = sum(d.calls_differ for d in differences)
        vcf_ids = {d.cgc.cohort.vcf_id for d in differences if d.cgc.cohort.vcf_id}
        self.stdout.write(f"WARNING: {self._differing_message(totals, differences)}")
        self.stdout.write(f"  {len(vcf_ids)} VCFs. Calls (samples_zygosity) differ for {calls_differ} variants, "
                          f"only the evidence (AD, DP, GQ...) for {totals['differing'] - calls_differ}")

        by_year = Counter()
        for d in differences:
            vcf = d.cgc.cohort.vcf
            by_year[vcf.date.year if vcf else "no VCF"] += d.variants
        self.stdout.write("  Variants by VCF import year: " + ", ".join(f"{year}: {n}" for year, n in
                                                                         sorted(by_year.items(), key=str)))

        largest = sorted(differences, key=lambda d: d.variants, reverse=True)[:self.LARGEST_SHOWN]
        self.stdout.write(f"  Largest {len(largest)} of {len(differences)} collections:")
        for d in largest:
            self.stdout.write(f"    {self._describe(d.cgc)}: {d.variants} variants, calls differ for {d.calls_differ}")

    @staticmethod
    def _describe(cgc: CohortGenotypeCollection) -> str:
        cohort = cgc.cohort
        if vcf := cohort.vcf:
            source = f"VCF {vcf.pk} '{vcf.name}' imported {vcf.date:%Y-%m-%d}"
        else:
            source = f"cohort {cohort.pk} '{cohort.name}'"
        return f"collection {cgc.pk} ({cgc.get_collection_type_display()}, {source})"

    def _write_csv(self, filename: str, differences: list[CollectionDifferences]):
        with open(filename, "w") as f:
            writer = csv.writer(f)
            writer.writerow(["collection_id", "collection_type", "cohort_id", "vcf_id", "vcf_name", "vcf_date",
                             "differing_variants", "calls_differ"])
            for d in differences:
                cohort = d.cgc.cohort
                vcf = cohort.vcf
                writer.writerow([d.cgc.pk, d.cgc.get_collection_type_display(), cohort.pk,
                                 vcf.pk if vcf else "", vcf.name if vcf else "", vcf.date.isoformat() if vcf else "",
                                 d.variants, d.calls_differ])
        self.stdout.write(f"Wrote {len(differences)} collections to {filename}")
