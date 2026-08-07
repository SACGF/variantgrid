"""
Migrate existing 'common' CohortGenotype partitions so they are valid for multiple gnomAD versions.

@see https://github.com/SACGF/variantgrid/issues/1582

VCFs are split into 'common'/'uncommon' partitions using a gnomAD AF>5 reference. The 'common' partition
can be skipped when filtering for rare variants, but only when its CohortGenotypeCommonFilterVersion lists the
analysis annotation's gnomAD version (@see CohortGenotypeCollection.get_annotation_kwargs).

Existing GRCh38 data was partitioned against a single gnomAD version (eg 4.0). To make those partitions valid
for the combined filter (eg 4.0 + 4.1) we move any variant that is NOT common (AF > gnomad_af_min) in the newly
added version(s) out of the common partition and into the uncommon partition - leaving the intersection behind -
then repoint the partition to the combined filter.
"""
import logging

from django.conf import settings
from django.core.management import BaseCommand
from django.db import transaction

from annotation.models import VariantAnnotationVersion
from library.django_utils.django_postgres import model_to_insert_sql
from library.utils.database_utils import run_sql
from snpdb.common_variants import get_common_filter
from snpdb.models import GenomeBuild, CohortGenotype, CohortGenotypeCollection, CohortGenotypeCollectionType

BATCH_SIZE = 1000


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--genome-build', help="Only process this genome build (default: all in settings)")
        parser.add_argument('--dry-run', action='store_true',
                            help="Report what would move without changing anything")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        if build_name := options["genome_build"]:
            build_names = [build_name]
        else:
            build_names = list(settings.VCF_IMPORT_COMMON_FILTERS)

        for build_name in build_names:
            genome_build = GenomeBuild.get_name_or_alias(build_name)
            self._migrate_genome_build(genome_build, dry_run=dry_run)

    def _migrate_genome_build(self, genome_build: GenomeBuild, dry_run: bool):
        target_filter = get_common_filter(genome_build)
        if target_filter is None:
            logging.info("%s: no common filter configured - skipping", genome_build)
            return

        target_versions = target_filter.gnomad_versions
        logging.info("%s: target common filter is %s (versions: %s)",
                     genome_build, target_filter, sorted(target_versions))

        common_cgc_qs = CohortGenotypeCollection.objects.filter(
            collection_type=CohortGenotypeCollectionType.COMMON,
            common_filter__genome_build=genome_build,
        ).exclude(common_filter=target_filter)

        for common_cgc in common_cgc_qs.iterator():
            self._migrate_common_collection(common_cgc, target_filter, target_versions, dry_run=dry_run)

    def _migrate_common_collection(self, common_cgc: CohortGenotypeCollection, target_filter, target_versions,
                                   dry_run: bool):
        genome_build = common_cgc.common_filter.genome_build
        af_min = common_cgc.common_filter.gnomad_af_min
        versions_to_enforce = target_versions - common_cgc.common_filter.gnomad_versions

        uncommon_cgc = common_cgc.uncommon  # the rare/uncommon partner partition
        logging.info("%s: enforcing versions %s (af > %s)", common_cgc, sorted(versions_to_enforce), af_min)

        moved = 0
        for version in versions_to_enforce:
            vav = VariantAnnotationVersion.objects.filter(gnomad=version,
                                                          genome_build=genome_build).order_by("pk").last()
            if vav is None:
                raise VariantAnnotationVersion.DoesNotExist(
                    f"Can't enforce gnomAD {version} on {common_cgc}: no VariantAnnotationVersion for it. "
                    "Annotate against this version first so the AF data is available.")

            # Variants in the common partition that are NOT common (AF > af_min) in this version -> must be moved
            move_qs = CohortGenotype.objects.filter(collection=common_cgc).exclude(
                variant__variantannotation__version=vav,
                variant__variantannotation__gnomad_af__gt=af_min)
            moved += self._move_cohort_genotypes(move_qs, uncommon_cgc, dry_run=dry_run)

        if dry_run:
            logging.info("%s: [dry-run] would move %d records and repoint to %s", common_cgc, moved, target_filter)
        else:
            common_cgc.common_filter = target_filter
            common_cgc.save(update_fields=["common_filter"])
            logging.info("%s: moved %d records, repointed to %s", common_cgc, moved, target_filter)

    @staticmethod
    def _move_cohort_genotypes(move_qs, uncommon_cgc: CohortGenotypeCollection, dry_run: bool) -> int:
        """ Move CohortGenotype rows into the uncommon partition. They live in different partitions, so we re-insert
            then delete the originals (@see common_variant_classified_task). Returns the number moved. """
        if dry_run:
            return move_qs.count()

        uncommon_table = uncommon_cgc.get_partition_table()
        total = 0
        batch = []
        delete_ids = []
        for cg in move_qs.iterator():
            delete_ids.append(cg.pk)
            cg.pk = None
            cg.collection = uncommon_cgc
            batch.append(cg)
            if len(batch) >= BATCH_SIZE:
                total += Command._flush_batch(batch, delete_ids, uncommon_table)
                batch, delete_ids = [], []
        total += Command._flush_batch(batch, delete_ids, uncommon_table)
        return total

    @staticmethod
    def _flush_batch(batch, delete_ids, uncommon_table) -> int:
        if not batch:
            return 0
        with transaction.atomic():
            for insert_sql in model_to_insert_sql(batch, ignore_fields=['id'], db_table=uncommon_table):
                run_sql(insert_sql)
            CohortGenotype.objects.filter(pk__in=delete_ids).delete()
        return len(batch)
