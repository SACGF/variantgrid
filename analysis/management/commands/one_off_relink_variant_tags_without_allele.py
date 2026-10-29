"""
Links variant tags left with no allele by the ClinGen request-level error crash fixed in #2354: tag liftover crashed
before assigning the allele, and nothing goes back to tags once their variant later gets one, so the tag never shows
on the same variant in other builds. Each tag goes through the tag liftover path again, which uses the variant's
existing VariantAllele (or calls ClinGen to make one) and queues liftover. --dry-run reports without writing.
Registered as a deploy-time ManualOperation (analysis/migrations/0004_one_off_relink_variant_tags_without_allele.py).
"""
from django.core.management.base import BaseCommand

from analysis.models.models_variant_tag import VariantTag
from analysis.tasks.variant_tag_tasks import _liftover_variant_tag


class Command(BaseCommand):
    help = "Link variant tags with no allele (left by the ClinGen error crash, #2354) and lift them over"
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report the tags without changing them")

    def handle(self, *args, **options):
        tags_qs = VariantTag.objects.filter(allele__isnull=True).select_related("variant", "genome_build")
        num_tags = tags_qs.count()
        if options["dry_run"]:
            self.stdout.write(f"{num_tags} variant tags have no allele")
            return

        for variant_tag in tags_qs.order_by("pk"):
            _liftover_variant_tag(variant_tag)
        self.stdout.write(f"Linked {num_tags} variant tags to their allele")
