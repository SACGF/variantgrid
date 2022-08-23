from collections import Counter
from traceback import print_exc

import pandas as pd
from django.core.management import BaseCommand
from django.db.models import Max
from pyhgvs import HGVSName

from classification.models import Classification
from snpdb.models import Variant


class Command(BaseCommand):
    """
        Classifications are matched to variants, and adding/changing transcript data (cdot), or matching algorithms
        can alter which transcript/variant is used.

        This can be re-triggered, but that may alter the classification's variant. This tool allows checking without
        altering the classification
    """

    def _get_last_modified(self):
        qs = Classification.objects.all()
        data = qs.aggregate(Max("modified"))
        return data["modified__max"]

    def handle(self, *args, **options):
        start_last_modified = self._get_last_modified()

        variant_diff_count = Counter()
        transcript_diff_count = Counter()
        diff_rows = []
        for c in Classification.objects.all().iterator(chunk_size=100):
            try:
                if v := c.variant:
                    existing_tuple = v.as_tuple()
                else:
                    existing_tuple = tuple()

                genome_build = c.get_genome_build()
                if genome_build.name == "GRCh37":
                    c_hgvs_name = c.chgvs_grch37_full
                elif genome_build.name == "GRCh38":
                    c_hgvs_name = c.chgvs_grch38_full
                else:
                    raise ValueError(f"Don't know how to get out cached chgvs from {c} ({genome_build=}) ")

                if c_hgvs_name:
                    existing_resolved_transcript = HGVSName(c_hgvs_name).transcript
                else:
                    existing_resolved_transcript = None

                vcfe = c._get_variant_coordinates_from_evidence()
                provided_transcript = c.transcript
                variant_tuple = vcfe.variant_coordinate
                transcript_accession = vcfe.transcript_accession

                if variant_tuple == existing_tuple:
                    variant_diff = ""
                elif not existing_tuple:
                    variant_diff = "gained variant match"
                elif not variant_tuple:
                    variant_diff = "lost variant match"
                else:
                    variant_diff = "variant matched changed"

                if transcript_accession == existing_resolved_transcript:
                    transcript_diff = ""
                elif not existing_resolved_transcript:
                    transcript_diff = "gained transcript"
                elif not transcript_accession:
                    transcript_diff = "lost transcript"
                else:
                    transcript_diff = "transcript changed"
                    if provided_transcript == transcript_accession:
                        transcript_diff += ": matched exact"
                    elif provided_transcript == existing_resolved_transcript:
                        transcript_diff += ": LOST EXACT MATCH!"

                if variant_diff or transcript_diff:
                    count_key = ", ".join([x for x in [variant_diff, transcript_diff] if x])

                    def _format_tuple(t):
                        if t:
                            return Variant.format_tuple(*t)
                        else:
                            return "."

                    #print(f"{existing_tuple=}")
                    diff_rows.append({
                        "classification": c.get_absolute_url(),
                        "provided_transcript": provided_transcript,
                        "old_resolved_transcript": existing_resolved_transcript,
                        "new_resolved_transcript": transcript_accession,
                        "transcript_diff": transcript_diff,
                        "old_variant": _format_tuple(existing_tuple),
                        "new_variant": _format_tuple(variant_tuple),
                        "variant_diff": variant_diff,
                    })

                if not variant_diff:
                    variant_diff = "No change"
                variant_diff_count[variant_diff] += 1

                if not transcript_diff:
                    transcript_diff = "No change"
                transcript_diff_count[transcript_diff] += 1
            except ValueError:
                print_exc()

        print("==== Variant matching ====")
        print(variant_diff_count)
        print("==== Transcripts: ==== ")
        print(transcript_diff_count)

        if diff_rows:
            filename = "classification_match_diff.csv"
            print(f"Writing diffs to '{filename}'")
            df = pd.DataFrame.from_records(diff_rows)
            df.to_csv(filename, index=False)

        end_last_modified = self._get_last_modified()
        if start_last_modified != end_last_modified:
            print(f"Beware - looks like someone edited a classification!!! {start_last_modified=} vs {end_last_modified=}")

