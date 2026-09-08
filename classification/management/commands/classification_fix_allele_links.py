"""
Re-home classifications whose clinical context or grouping sits under a different Allele than the
classification does - what an allele merge or the #1361 VariantAllele dedupe leaves behind.

Registered as a ManualOperation by snpdb/migrations/0259_variantallele_unique_variant_build.py
"""
from django.core.management import BaseCommand

from classification.models.clinical_context_utils import (
    classifications_needing_rehoming,
    rehome_classifications,
)


class Command(BaseCommand):

    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would be re-homed and stop")

    def handle(self, *args, **options):
        classification_qs = classifications_needing_rehoming()
        num_mismatched = classification_qs.count()
        if not num_mismatched:
            self.stdout.write("All classifications are filed under their own allele")
            return

        if options["dry_run"]:
            for classification in classification_qs:
                self.stdout.write(f"{classification.cr_lab_id}: allele {classification.allele_id}, "
                                  f"clinical context allele {classification.clinical_context.allele_id if classification.clinical_context else None}")
            self.stdout.write(f"{num_mismatched} classifications would be re-homed")
            return

        num_rehomed = rehome_classifications(classification_qs, force_recalc_text="allele merged (#1361)")
        self.stdout.write(f"Re-homed {num_rehomed} classifications")
