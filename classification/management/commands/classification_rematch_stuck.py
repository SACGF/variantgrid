"""Re-match ImportedAlleleInfo records left in Processing by an import that died.

Variant matching for a batch runs as upload pipelines (@see classification.classification_import); if one
never completes, its records sit in Processing with nothing left to restart them - the FileUpload rows have
no link back to the ClassificationImport. Re-matching is the way back: it derives the coordinates again and
runs fresh pipelines. This is the admin's "Re-Match Soft" action over a queryset you can pick from the CLI.
"""
from datetime import timedelta

from django.core.management.base import BaseCommand
from django.utils import timezone

from classification.classification_import import reattempt_variant_matching
from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from library.guardian_utils import admin_bot


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--older-than-hours', type=int, default=24,
                            help="Only records that have been Processing this long (default 24)")
        parser.add_argument('--gene-level', action='store_true',
                            help="Only records whose coordinate is gene-level (fusion, splice, CNV)")
        parser.add_argument('--dry-run', action='store_true', help="List what would be re-matched")

    def handle(self, *args, **options):
        older_than_hours = options['older_than_hours']
        gene_level_only = options['gene_level']
        dry_run = options['dry_run']

        qs = ImportedAlleleInfo.objects.filter(status=ImportedAlleleInfoStatus.PROCESSING,
                                               modified__lte=timezone.now() - timedelta(hours=older_than_hours))
        if gene_level_only:
            gene_level_pks = []
            for allele_info in qs:
                variant_coordinate = allele_info.variant_coordinate_obj
                if variant_coordinate and variant_coordinate.is_gene_level:
                    gene_level_pks.append(allele_info.pk)
            qs = ImportedAlleleInfo.objects.filter(pk__in=gene_level_pks)

        if dry_run:
            for allele_info in qs.order_by("pk"):
                print(f"{allele_info.pk}\t{allele_info.imported_genome_build_patch_version}\t{allele_info.imported_c_hgvs}")
            print(f"Would re-match {qs.count()} records")
            return

        print(f"Re-matching {qs.count()} records")
        for allele_info in qs:
            allele_info.update_variant_coordinate()
            allele_info.refresh_and_save(force_update=True)
            allele_info.classification_import = None
            allele_info.status = ImportedAlleleInfoStatus.PROCESSING
            allele_info.save()

        queued = reattempt_variant_matching(admin_bot(), qs, False)
        print(f"Queued {queued} records for matching")
