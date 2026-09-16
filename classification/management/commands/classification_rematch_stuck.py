"""Re-match ImportedAlleleInfo records an import that died left behind.

Variant matching for a batch runs as upload pipelines (@see classification.classification_import). A pipeline
that dies leaves its records in Processing, and one that reaches the linking step without the variants it
expected marks them Failed - either way there is nothing to restart, as the FileUpload rows have no link back
to the ClassificationImport. Re-matching is the way back: it derives the coordinates again and runs fresh
pipelines. This is the admin's "Re-Match Soft" action over a queryset you can pick from the CLI.
"""
from datetime import timedelta

from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone

from classification.classification_import import reattempt_variant_matching
from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from library.guardian_utils import admin_bot


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--status', default=ImportedAlleleInfoStatus.PROCESSING,
                            help="Comma separated ImportedAlleleInfoStatus codes (default 'P' - Processing, "
                                 "'F' is Failed)")
        parser.add_argument('--older-than-hours', type=int, default=24,
                            help="Only records untouched for this long (default 24)")
        parser.add_argument('--gene-level', action='store_true',
                            help="Only records whose coordinate is gene-level (fusion, splice, CNV)")
        parser.add_argument('--dry-run', action='store_true', help="List what would be re-matched")

    def handle(self, *args, **options):
        statuses = options['status'].split(',')
        if unknown := set(statuses) - set(dict(ImportedAlleleInfoStatus.choices)):
            raise CommandError(f"Unknown --status {sorted(unknown)}")
        older_than_hours = options['older_than_hours']
        gene_level_only = options['gene_level']
        dry_run = options['dry_run']

        qs = ImportedAlleleInfo.objects.filter(status__in=statuses,
                                               modified__lte=timezone.now() - timedelta(hours=older_than_hours))
        if gene_level_only:
            pks = [ai.pk for ai in qs if (vc := ai.variant_coordinate_obj) and vc.is_gene_level]
        else:
            pks = list(qs.values_list("pk", flat=True))
        # By pk from here on - the loop below moves records out of the status the filter selected on
        qs = ImportedAlleleInfo.objects.filter(pk__in=pks)

        if dry_run:
            for allele_info in qs.order_by("pk"):
                print(f"{allele_info.pk}\t{allele_info.get_status_display()}\t"
                      f"{allele_info.imported_genome_build_patch_version}\t{allele_info.imported_c_hgvs}")
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
