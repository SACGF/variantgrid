"""
Removes duplicate CohortGenotype rows - more than one row for a variant in a CohortGenotypeCollection - left by VCF
imports from before per-batch duplicate removal, where nearby multi-allelic records normalised to the same variant
(#2443). A row is deleted only when its genotype data is identical to the variant's lowest id row, then the zygosity
counts of those variants are recounted: the count UPDATE joined to only one of the rows, so what was added can't be
subtracted.

A row whose genotype differs from the lowest id row is left in place and reported as a warning: one VCF giving the
same variant two genotypes is an import / normalisation bug, and picking one would silently drop a call. Each
collection's partition is searched and fixed on its own, in its own transaction. --dry-run reports without writing.
Registered as a deploy-time ManualOperation (snpdb/migrations/0005_dedupe_cohort_genotypes.py).
"""
from collections import defaultdict
from dataclasses import dataclass, field

from django.core.management.base import BaseCommand
from django.db import connection, transaction
from django.db.models import F, Q

from library.django_utils.database_utils import run_sql
from library.log_utils import report_message
from snpdb.enums import CohortGenotypeCollectionType
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.models.models_zygosity_counts import (
    VariantZygosityCountCollection,
    VariantZygosityCountForSample,
    VariantZygosityCountForVCF,
)
from snpdb.variant_zygosity_count import (
    _get_update_sql_and_params,
    _sample_excluded_from_variant_zygosity_count,
)

GENOTYPE_COLUMNS = ["ref_count", "het_count", "hom_count", "unk_count", "samples_zygosity", "samples_allele_depth",
                    "samples_allele_frequency", "samples_read_depth", "samples_genotype_quality",
                    "samples_phred_likelihood", "samples_filters"]


@dataclass
class DuplicateGenotypes:
    variant_id: int
    keep_id: int
    identical_ids: list[int] = field(default_factory=list)
    differing_ids: list[int] = field(default_factory=list)


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
               ({duplicate_row}) IS NOT DISTINCT FROM ({kept_row})
        FROM duplicated
        JOIN "{table}" kept ON kept.id = duplicated.keep_id
        JOIN "{table}" duplicate ON duplicate.variant_id = duplicated.variant_id AND duplicate.id <> duplicated.keep_id
        ORDER BY duplicated.variant_id, duplicate.id
    """
    duplicates: dict[int, DuplicateGenotypes] = {}
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for variant_id, keep_id, duplicate_id, identical in cursor.fetchall():
            if not (duplicate := duplicates.get(variant_id)):
                duplicate = duplicates[variant_id] = DuplicateGenotypes(variant_id, keep_id)
            (duplicate.identical_ids if identical else duplicate.differing_ids).append(duplicate_id)
    return list(duplicates.values())


def _vcf_ids_with_variants(variant_ids: list[int]) -> set[int]:
    """ VCFs whose current genotype collection (Cohort.cohort_genotype_collection - what the count SQL reads)
        has a row for any of these variants """
    sql = "SELECT DISTINCT collection_id FROM snpdb_cohortgenotype WHERE variant_id = ANY(%s)"
    with connection.cursor() as cursor:
        cursor.execute(sql, [variant_ids])
        cgc_ids = [row[0] for row in cursor.fetchall()]
    cgc_qs = CohortGenotypeCollection.objects.filter(pk__in=cgc_ids,
                                                     collection_type=CohortGenotypeCollectionType.UNCOMMON,
                                                     cohort_version=F("cohort__version"),
                                                     data_archived_date__isnull=True,
                                                     cohort__vcf__isnull=False,
                                                     cohort__vcf__data_archived_date__isnull=True)
    return set(cgc_qs.values_list("cohort__vcf_id", flat=True))


def _add_count_for_variants(collection: VariantZygosityCountCollection, vcf, variant_ids: list[int], sample=None):
    """ The normal count UPDATE (which ends in its WHERE clause), restricted to these variants """
    sql, params = _get_update_sql_and_params(collection, vcf, '+', sample=sample)
    run_sql(sql + " AND snpdb_cohortgenotype.variant_id = ANY(%s)", [*params, variant_ids])


def recount_variant_zygosity_counts(variant_ids: list[int]):
    """ Recalculate these variants' counts from scratch: zero them, then re-add every VCF and sample counted in
        each collection, restricted to these variants. Like the other count updates, nothing else should be
        counting at the same time. """
    vcf_ids = _vcf_ids_with_variants(variant_ids)
    for collection in VariantZygosityCountCollection.objects.filter(data_archived_date__isnull=True):
        counted = Q(collection=collection, count_complete__isnull=False, deleted__isnull=True)
        counted_vcfs = VariantZygosityCountForVCF.objects.filter(counted, vcf__in=vcf_ids,
                                                                 is_split_to_sample_counts=False)
        counted_samples = VariantZygosityCountForSample.objects.filter(counted, sample__vcf__in=vcf_ids)
        with transaction.atomic():
            run_sql(f'UPDATE "{collection.get_partition_table()}" SET ref_count = 0, het_count = 0, hom_count = 0 '
                    'WHERE variant_id = ANY(%s)', [variant_ids])
            for vzc_vcf in counted_vcfs.select_related("vcf"):
                _add_count_for_variants(collection, vzc_vcf.vcf, variant_ids)
            for vzc_sample in counted_samples.select_related("sample__vcf"):
                sample = vzc_sample.sample
                if not _sample_excluded_from_variant_zygosity_count(sample):
                    _add_count_for_variants(collection, sample.vcf, variant_ids, sample=sample)
            collection.bump_data_version()


class Command(BaseCommand):
    help = ("Delete duplicate CohortGenotype rows identical to the variant's first row, recount their zygosity "
            "counts, and warn about duplicates with differing genotypes")
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would be deleted, change nothing")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        existing_tables = set(connection.introspection.table_names())
        cgc_qs = CohortGenotypeCollection.objects.filter(data_archived_date__isnull=True,
                                                         marked_for_deletion=False).order_by("pk")
        totals = defaultdict(int)
        for cgc in cgc_qs.select_related("cohort"):
            if cgc.get_partition_table() not in existing_tables:
                continue
            with transaction.atomic():
                duplicates = find_duplicate_genotypes(cgc)
                if not duplicates:
                    continue
                delete_ids = [pk for d in duplicates for pk in d.identical_ids]
                differing = [d for d in duplicates if d.differing_ids]
                if delete_ids and not dry_run:
                    run_sql(f'DELETE FROM "{cgc.get_partition_table()}" WHERE id = ANY(%s)', [delete_ids])
                    recount_variant_zygosity_counts([d.variant_id for d in duplicates if d.identical_ids])

            totals["variants"] += len(duplicates)
            totals["deleted"] += len(delete_ids)
            totals["differing"] += len(differing)
            verb = "would delete" if dry_run else "deleted"
            self.stdout.write(f"CohortGenotypeCollection {cgc.pk} ({cgc.cohort.name}): {len(duplicates)} duplicated "
                              f"variants, {verb} {len(delete_ids)} identical rows, {len(differing)} variants with "
                              "differing genotypes left in place")
            if differing:
                self._warn_differing(cgc, differing, dry_run)

        mode = " (dry run)" if dry_run else ""
        self.stdout.write(f"Total{mode}: {totals['variants']} duplicated variants, {totals['deleted']} identical rows "
                          f"{'to delete' if dry_run else 'deleted'}, {totals['differing']} variants with differing "
                          "genotypes left in place")

    def _warn_differing(self, cgc: CohortGenotypeCollection, differing: list[DuplicateGenotypes], dry_run: bool):
        message = (f"CohortGenotypeCollection {cgc.pk} has {len(differing)} variants with more than one genotype row "
                   "that disagree: the same variant with a different genotype inside one VCF is an import / "
                   "normalisation bug. Left in place, not deduplicated")
        self.stdout.write(f"  WARNING: {message}")
        for d in differing:
            self.stdout.write(f"    collection {cgc.pk} variant {d.variant_id}: rows {[d.keep_id, *d.differing_ids]}")
        if not dry_run:
            report_message(message, extra_data={"target": f"cohort_genotype_collection={cgc.pk}",
                                                "cohort": cgc.cohort.name,
                                                "differing_variants": len(differing)})
