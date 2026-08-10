from django.core.management.base import BaseCommand
from django.db.models import Q

from classification.models import ImportedAlleleInfo


class Command(BaseCommand):
    """ g.HGVS only submissions used to be validated as though a c.HGVS was expected, which left them with
        errors that excluded them from exports - see https://github.com/SACGF/variantgrid/issues/1063 """

    def handle(self, *args, **options):
        g_hgvs_only_qs = ImportedAlleleInfo.objects.filter(
            Q(imported_c_hgvs__isnull=True) | Q(imported_c_hgvs=""),
            Q(imported_transcript__isnull=True) | Q(imported_transcript=""),
        )
        total = g_hgvs_only_qs.count()
        self.stdout.write(f"Re-validating {total} g.HGVS only ImportedAlleleInfo records")

        changed = 0
        for allele_info in g_hgvs_only_qs.iterator():
            previous_validation = allele_info.latest_validation
            previous_tags = previous_validation.validation_tags if previous_validation else None
            previous_include = previous_validation.include if previous_validation else None

            allele_info.apply_validation(force_update=True)
            allele_info.save()

            latest_validation = allele_info.latest_validation
            if latest_validation.validation_tags != previous_tags or latest_validation.include != previous_include:
                changed += 1

        self.stdout.write(self.style.SUCCESS(f"Re-validation complete, {changed} of {total} records changed"))
