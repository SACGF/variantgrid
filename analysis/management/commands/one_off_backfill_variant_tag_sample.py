"""
Fill in VariantTag.sample for taggings made before the field existed.

The sample is what puts a tagging into a case's classify queue (@see the Classify & Report tab). Tags made
in an analysis can usually say which sample they were about without anyone being asked - the node had a single
sample, or only one of the analysis's samples carries the variant. Anything still ambiguous is left null and
resolved at classification time, where a sample dropdown already exists.

@see https://github.com/SACGF/variantgrid_sapath/issues/246
"""
import logging

from django.core.management import BaseCommand

from analysis.models import VariantTag
from analysis.variant_tag_operations import get_sample_for_variant_tag


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="Report what would be set without changing anything")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        qs = VariantTag.objects.filter(sample__isnull=True, analysis__isnull=False) \
            .select_related("analysis", "node", "variant", "allele").order_by("pk")

        resolved = 0
        ambiguous = 0
        for variant_tag in qs.iterator():
            sample = get_sample_for_variant_tag(variant_tag)
            if sample is None:
                ambiguous += 1
                continue
            resolved += 1
            logging.info("VariantTag %s (%s) -> sample %s", variant_tag.pk, variant_tag.tag_id, sample)
            if not dry_run:
                VariantTag.objects.filter(pk=variant_tag.pk).update(sample=sample)

        prefix = "Would set" if dry_run else "Set"
        self.stdout.write(f"{prefix} sample on {resolved} taggings, left {ambiguous} ambiguous")
