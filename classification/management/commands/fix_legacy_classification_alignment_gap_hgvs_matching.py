from django.core.management import BaseCommand
from django.db.models import Q

from classification.classification_import import reattempt_variant_matching
from classification.models import Classification
from genes.models import TranscriptVersion
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):

        # Do transcript as well (doesn't matter if we have false positives will just spend a bit longer rematching)
        transcripts_with_gaps = set()
        for tv in TranscriptVersion.objects.filter(alignment_gap=True):
            transcripts_with_gaps.add(tv.accession)
            transcripts_with_gaps.add(tv.transcript_id)

        gene_symbols = set()
        classification_to_rematch = set()
        for c in Classification.objects.filter(variant__isnull=False):
            if c.transcript in transcripts_with_gaps:
                classification_to_rematch.add(c.pk)
                if gene_symbol := c.get("gene_symbol"):
                    gene_symbols.add(gene_symbol)

        gene_symbols_str = ", ".join(sorted(gene_symbols))
        print(f"Rematching {len(classification_to_rematch)} classifications in transcripts with gaps. Symbols: {gene_symbols_str}")

        classification_qs = Classification.objects.filter(Q(variant__isnull=True) | Q(pk__in=classification_to_rematch))
        num_classifications = classification_qs.count()
        print(f"Rematching {num_classifications} classifications total (also includes variant=None)")

        valid_record_count, invalid_record_count = reattempt_variant_matching(admin_bot(), classification_qs)
        print(f"{valid_record_count=}, {invalid_record_count=}")
