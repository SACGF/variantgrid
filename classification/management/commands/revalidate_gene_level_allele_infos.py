from django.core.management.base import BaseCommand

from classification.models import ImportedAlleleInfo
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME


class Command(BaseCommand):
    """ Gene-level submissions used to be validated as though a c.HGVS was expected, which tagged
        transcript_type_not_supported as fatal and excluded every one of them from exports -
        see https://github.com/SACGF/variantgrid/issues/1835 """
    category = "one-off"

    def handle(self, *args, **options):
        gene_level_qs = ImportedAlleleInfo.objects.filter(
            variant_coordinate__startswith=f"{GENE_LEVEL_CONTIG_NAME}:")
        total = gene_level_qs.count()
        self.stdout.write(f"Re-validating {total} gene-level ImportedAlleleInfo records")

        changed = 0
        for allele_info in gene_level_qs.iterator():
            previous_validation = allele_info.latest_validation
            previous_tags = previous_validation.validation_tags if previous_validation else None
            previous_include = previous_validation.include if previous_validation else None

            allele_info.apply_validation(force_update=True)
            allele_info.save()

            latest_validation = allele_info.latest_validation
            if latest_validation.validation_tags != previous_tags or latest_validation.include != previous_include:
                changed += 1

        self.stdout.write(self.style.SUCCESS(f"Re-validation complete, {changed} of {total} records changed"))
