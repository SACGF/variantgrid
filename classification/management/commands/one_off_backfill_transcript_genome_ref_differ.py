from collections import defaultdict

from django.core.management.base import BaseCommand
from django.db.models import Q

from classification.models.classification_variant_info_models import ImportedAlleleInfo
from classification.tasks.transcript_genome_ref_differ_backfill_task import (
    backfill_transcript_genome_ref_differ,
    backfill_transcript_genome_ref_differ_task,
)
from genes.hgvs.biocommons_hgvs.data_provider import is_contig
from genes.hgvs.hgvs import HGVSComponents
from genes.hgvs.hgvs_converter import HGVSConverterType


class Command(BaseCommand):
    """ ImportedAlleleInfo.transcript_genome_ref_differ is set during matching, so records matched before the check
        existed are null. Only biocommons resolution can compare the two sequences, so those records are grouped by
        transcript and build and queued a transcript at a time (#1266) """
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("--sync", action="store_true", help="Run in this process rather than queuing celery tasks")

    def handle(self, *args, **options):
        allele_info_qs = ImportedAlleleInfo.objects.filter(
            transcript_genome_ref_differ__isnull=True,
            variant_coordinate__isnull=False,
            hgvs_converter_version__used_converter_type=HGVSConverterType.BIOCOMMONS_HGVS.name,
        ).exclude(Q(imported_c_hgvs__isnull=True) | Q(imported_c_hgvs=""))

        ids_by_transcript_and_build: dict[tuple[str, int], list[int]] = defaultdict(list)
        values = allele_info_qs.values_list("pk", "imported_transcript", "imported_c_hgvs",
                                            "imported_genome_build_patch_version__genome_build_id")
        for pk, imported_transcript, imported_c_hgvs, genome_build_id in values.iterator():
            transcript_accession = imported_transcript or HGVSComponents(imported_c_hgvs).transcript
            if transcript_accession and not is_contig(transcript_accession):
                ids_by_transcript_and_build[(transcript_accession, genome_build_id)].append(pk)

        total = sum(len(ids) for ids in ids_by_transcript_and_build.values())
        self.stdout.write(f"{total} allele infos across {len(ids_by_transcript_and_build)} transcript/build groups")

        if options["sync"]:
            checked = differ = 0
            for (transcript_accession, _), allele_info_ids in ids_by_transcript_and_build.items():
                group_checked, group_differ = backfill_transcript_genome_ref_differ(transcript_accession, allele_info_ids)
                checked += group_checked
                differ += group_differ
            self.stdout.write(f"Checked {checked}, {differ} differ")
        else:
            for (transcript_accession, _), allele_info_ids in ids_by_transcript_and_build.items():
                backfill_transcript_genome_ref_differ_task.delay(transcript_accession, allele_info_ids)
            self.stdout.write(f"Queued {len(ids_by_transcript_and_build)} tasks on db_workers")
