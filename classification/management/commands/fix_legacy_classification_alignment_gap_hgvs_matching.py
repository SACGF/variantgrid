from collections import defaultdict

from django.core.management import BaseCommand

from classification.classification_import import reattempt_variant_matching
from classification.models import Classification
from genes.models import TranscriptVersion
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):

        # Get a list of transcripts/genome build used by classifications
        build_transcript_classification_ids = defaultdict(lambda: defaultdict(list))
        for c in Classification.objects.all():
            if transcript := c.transcript:
                build_transcript_classification_ids[c.get_genome_build()][transcript].append(c.pk)

        # We have to deal with bumping transcripts etc, so probably the best way to go is find the ones that would be
        # OK and then anything else is assumed to need re-matching
        # This will also get LRGs
        classification_ids_to_rematch = []
        for genome_build, transcript_classification_ids in build_transcript_classification_ids.items():
            # Things that have gaps, NR etc will not be returned so will forced to rematch
            previously_good_tv_qs = TranscriptVersion.objects.filter(genome_build=genome_build,
                                                                     data__chrom__isnull=False,
                                                                     data__cdna_match__isnull=True,
                                                                     data__partial__isnull=True,
                                                                     transcript__identifier__startswith='NM_')
            previously_good_transcript_accessions = set()
            for tv in previously_good_tv_qs:
                if not tv.alignment_gap:
                    previously_good_transcript_accessions.add(tv.accession)

            for t, classifications in transcript_classification_ids.items():
                if t not in previously_good_transcript_accessions:
                    classification_ids_to_rematch.extend(classifications)

        classification_qs = Classification.objects.filter(pk__in=classification_ids_to_rematch)
        num_classifications = classification_qs.count()
        print(f"Rematching {num_classifications} classifications")

        valid_record_count, invalid_record_count = reattempt_variant_matching(admin_bot(), classification_qs)
        print(f"{valid_record_count=}, {invalid_record_count=}")
