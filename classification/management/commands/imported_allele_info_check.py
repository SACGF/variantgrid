from collections import Counter
from traceback import print_exc

import pandas as pd
from django.core.management import BaseCommand
from django.db.models import Max

from classification.models import ImportedAlleleInfo
from genes.hgvs import HGVSMatcher, HGVSConverterType
from snpdb.models import Variant, GenomeBuild


class Command(BaseCommand):
    """
        Classifications are matched to variants, and adding/changing transcript data (cdot), or matching algorithms
        can alter which transcript/variant is used.

        This can be re-triggered, but that may alter the classification's variant. This tool allows checking without
        altering the classification
    """

    def _get_last_modified(self):
        qs = ImportedAlleleInfo.objects.all()
        data = qs.aggregate(Max("modified"))
        return data["modified__max"]

    def handle(self, *args, **options):
        start_last_modified = self._get_last_modified()

        variant_diff_count = Counter()
        transcript_diff_count = Counter()
        diff_rows = []
        iai: ImportedAlleleInfo
        hgvs_matchers_by_build = {}
        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build, hgvs_converter_type=HGVSConverterType.COMBO)
            hgvs_matchers_by_build[genome_build] = matcher

        for iai in ImportedAlleleInfo.objects.filter(imported_c_hgvs__isnull=False).iterator(chunk_size=100):
            try:
                genome_build = iai.imported_genome_build
                matcher = hgvs_matchers_by_build[genome_build]
                provided_transcript_accession = matcher.get_transcript_accession(iai.imported_c_hgvs)
                vcd = matcher.get_variant_tuple_used_transcript_kind_method_and_matches_reference(iai.imported_c_hgvs)
                current_variant_coordinate, current_transcript_accession, _kind, _method, _matches_ref = vcd

                rvi = iai[iai.imported_genome_build]
                if rvi.variant:
                    existing_variant_coordinate = rvi.variant.coordinate
                else:
                    existing_variant_coordinate = None

                resolved_transcript_accession = rvi.transcript_version.accession
                # resolved_c_hgvs_name = rvi.c_hgvs_full

                if current_variant_coordinate == existing_variant_coordinate:
                    variant_diff = ""
                elif not existing_variant_coordinate:
                    variant_diff = "gained variant match"
                elif not current_variant_coordinate:
                    variant_diff = "lost variant match"
                else:
                    variant_diff = "variant matched changed"

                if provided_transcript_accession == resolved_transcript_accession:
                    transcript_diff = ""
                elif not resolved_transcript_accession:
                    transcript_diff = "gained transcript"
                elif not provided_transcript_accession:
                    transcript_diff = "lost transcript"
                else:
                    transcript_diff = "transcript changed"
                    if provided_transcript_accession == current_transcript_accession:
                        transcript_diff += ": matched exact"
                    elif provided_transcript_accession == resolved_transcript_accession:
                        transcript_diff += ": LOST EXACT MATCH!"

                if variant_diff or transcript_diff:
                    def _format_tuple(t):
                        if t:
                            return Variant.format_tuple(*t)
                        else:
                            return "."

                    #print(f"{existing_tuple=}")
                    diff_rows.append({
                        "imported_allele_info": iai.get_absolute_url(),
                        "provided_transcript_accession": provided_transcript_accession,
                        "old_resolved_transcript": resolved_transcript_accession,
                        "new_resolved_transcript": current_transcript_accession,
                        "transcript_diff": transcript_diff,
                        "old_variant": _format_tuple(existing_variant_coordinate),
                        "new_variant": _format_tuple(current_variant_coordinate),
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
