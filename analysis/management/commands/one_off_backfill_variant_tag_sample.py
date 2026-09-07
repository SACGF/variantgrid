"""
Fill in VariantTag.sample for taggings made before the field existed.

The sample is what puts a tagging into a case's classify queue (@see the Classify & Report tab). A tag made in
an analysis can usually say which sample it was about without anyone being asked - the node it was made in
knows the study's proband. Anything still ambiguous is left null and picked at classification time, where a
sample dropdown already exists.

Taggings are done an analysis at a time so the node graph - which is what the answer actually comes from - is
loaded once for the thousands of taggings that share it. Only taggings without a sample are selected, so an
interrupted run picks up where it stopped.

@see https://github.com/SACGF/variantgrid_sapath/issues/246
"""
import logging

from django.core.management import BaseCommand

from analysis.models import Analysis, VariantTag
from analysis.variant_tag_operations import get_proband_sample_by_node_id

BATCH_SIZE = 1000


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="Report what would be set without changing anything")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        base_qs = VariantTag.objects.filter(sample__isnull=True, analysis__isnull=False, node__isnull=False)
        analysis_ids = list(base_qs.order_by("analysis_id").values_list("analysis_id", flat=True).distinct())

        filled = 0
        ambiguous = 0
        for i, analysis_id in enumerate(analysis_ids, start=1):
            analysis = Analysis.objects.get(pk=analysis_id)
            proband_sample_by_node_id = get_proband_sample_by_node_id(analysis)

            to_update = []
            for variant_tag in base_qs.filter(analysis=analysis).order_by("pk").iterator():
                sample = proband_sample_by_node_id.get(variant_tag.node_id)
                if sample is None:
                    ambiguous += 1
                    continue
                variant_tag.sample = sample
                to_update.append(variant_tag)
                if len(to_update) >= BATCH_SIZE:
                    filled += self._write(to_update, dry_run)
                    to_update = []
            filled += self._write(to_update, dry_run)

            logging.info("Analysis %s (%d/%d): %d filled, %d ambiguous so far",
                         analysis_id, i, len(analysis_ids), filled, ambiguous)

        prefix = "Would set" if dry_run else "Set"
        self.stdout.write(f"{prefix} sample on {filled} taggings, left {ambiguous} ambiguous")

    @staticmethod
    def _write(variant_tags: list[VariantTag], dry_run: bool) -> int:
        if variant_tags and not dry_run:
            VariantTag.objects.bulk_update(variant_tags, ["sample"])
        return len(variant_tags)
