from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Optional
from django.core.management import BaseCommand
from django.db.models import Max

from classification.models import ImportedAlleleInfo
from genes.hgvs import HGVSMatcher, HGVSConverterType
from genes.models import TranscriptVersion, TranscriptParts
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row
from snpdb.models import Variant, GenomeBuild, VariantCoordinate


class VariantChange(int, Enum):
    same = 0
    gained = 1
    lost = 2
    changed = 3

    @staticmethod
    def change_for(resolved: Optional[VariantCoordinate], updated: Optional[VariantCoordinate]):
        if resolved == updated:
            return VariantChange.same
        elif (not resolved) and updated:
            return VariantChange.gained
        elif resolved and (not updated):
            return VariantChange.lost
        else:
            return VariantChange.changed

    def __str__(self):
        if self == VariantChange.same:
            return "same"
        elif self == VariantChange.gained:
            return "match gained"
        elif self == VariantChange.lost:
            return "match lost"
        elif self == VariantChange.changed:
            return "changed"

    def __repr__(self):
        return str(self)


class TranscriptChange(int, Enum):
    same = 0
    gained = 1
    lost = 2
    changed = 3
    changed_from_exact = 4
    changed_to_exact = 5

    @staticmethod
    def change_for(provided: Optional[TranscriptParts], resolved: Optional[TranscriptParts], updated: Optional[TranscriptParts]):
        if resolved == updated:
            return TranscriptChange.same
        elif (not resolved) and updated:
            return TranscriptChange.gained
        elif resolved and (not updated):
            return TranscriptChange.lost
        else:
            if provided:
                if provided == resolved:
                    return TranscriptChange.changed_from_exact
                if provided == updated:
                    return TranscriptChange.changed_to_exact
            return TranscriptChange.changed

    def __str__(self):
        if self == TranscriptChange.same:
            return "same"
        elif self == TranscriptChange.gained:
            return "match gained"
        elif self == TranscriptChange.lost:
            return "match lost"
        elif self == TranscriptChange.changed:
            return "changed"
        elif self == TranscriptChange.changed_from_exact:
            return "changed - lost exact"
        elif self == TranscriptChange.changed_to_exact:
            return "changed - to exact"

    def __repr__(self):
        return str(self)

@dataclass
class HgvsSummary(ExportRow):
    transcript: Optional[TranscriptParts] = None
    variant_coordinate: Optional[VariantCoordinate] = None
    error_str: Optional[str] = None

    @export_column("variant_coordinate")
    def _variant_coordinate(self):
        if self.variant_coordinate:
            return Variant.format_tuple(*self.variant_coordinate)

    @export_column("transcript")
    def _transcript_export(self):
        return self.transcript

    @export_column("error")
    def _error_str(self):
        return self.error_str


@dataclass
class ChgvsDiff(ExportRow):

    imported_allele_info: ImportedAlleleInfo
    provided: HgvsSummary
    resolved: HgvsSummary
    updated: HgvsSummary

    @property
    def has_difference(self):
        if self.resolved.variant_coordinate != self.updated.variant_coordinate:
            return True
        if self.resolved.transcript != self.resolved.transcript:
            return True
        return False

    @export_column("$site_name URL")
    def _allele_id(self):
        return get_url_from_view_path(self.imported_allele_info.get_absolute_url())

    @export_column()
    def variant_change(self):
        return VariantChange.change_for(self.resolved.variant_coordinate, self.updated.variant_coordinate)

    @export_column()
    def transcript_change(self):
        return TranscriptChange.change_for(self.provided.transcript, self.resolved.transcript, self.updated.transcript)

    @export_column()
    def _imported_c_hgvs(self):
        return self.imported_allele_info.imported_c_hgvs

    @export_column()
    def _imported_genome_build(self):
        return self.imported_allele_info.imported_genome_build

    @export_column("provided", sub_data=HgvsSummary)
    def _provided(self):
        return self.provided

    @export_column("resolved", sub_data=HgvsSummary)
    def _resolved(self):
        return self.resolved

    @export_column("updated", sub_data=HgvsSummary)
    def _updated(self):
        return self.updated


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
        iai: ImportedAlleleInfo
        hgvs_matchers_by_build = {}
        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build, hgvs_converter_type=HGVSConverterType.COMBO)
            hgvs_matchers_by_build[genome_build] = matcher

        filename = "classification_match_diff.csv"
        with open(filename, "w") as out:
            out.write(delimited_row(ChgvsDiff.csv_header()))
            iai_qs = ImportedAlleleInfo.objects.filter(imported_c_hgvs__isnull=False)
            # Skip these as they won't match with HGVSMatcher
            # iai_qs = iai_qs.exclude(message="HGVS matched by 'ClinGen Allele Registry'")
            for iai in iai_qs.iterator(chunk_size=100):
                genome_build = iai.imported_genome_build
                matcher = hgvs_matchers_by_build[genome_build]
                imported_c_hgvs = iai.imported_c_hgvs

                provided = HgvsSummary()
                resolved = HgvsSummary()
                updated = HgvsSummary()
                diff = ChgvsDiff(
                    imported_allele_info=iai,
                    provided=provided,
                    resolved=resolved,
                    updated=updated
                )

                # PROVIDED
                try:
                    provided.full_c_hgvs = imported_c_hgvs
                    provided.transcript = matcher.get_transcript_parts(imported_c_hgvs)
                except Exception as ex:
                    provided.error_str = str(ex)

                # RESOLVED
                try:
                    if variant_info := iai[iai.imported_genome_build]:
                        resolved.full_c_hgvs = variant_info.c_hgvs
                        if variant := variant_info.variant:
                            resolved.variant_coordinate = variant.coordinate
                            if transcript_version := variant_info.transcript_version:
                                resolved.transcript = transcript_version.as_parts
                        if error := variant_info.error:
                            resolved.error_str = error
                except Exception as ex:
                    resolved.error_str = str(ex)

                # UPDATED
                try:
                    vcd = matcher.get_variant_tuple_used_transcript_kind_method_and_matches_reference(imported_c_hgvs)
                    updated.variant_coordinate = vcd.variant_coordinate
                    updated.transcript = TranscriptVersion.get_transcript_id_and_version(vcd.transcript_accession)
                except Exception as ex:
                    updated.error_str = str(ex)

                if diff.has_difference:
                    out.write(delimited_row(diff.to_csv()))

                variant_diff_count[diff.variant_change()] += 1
                transcript_diff_count[diff.transcript_change()] += 1

        print("==== Variant matching ====")
        print(dict(variant_diff_count))
        print("==== Transcripts: ==== ")
        print(dict(transcript_diff_count))

        end_last_modified = self._get_last_modified()
        if start_last_modified != end_last_modified:
            print(f"Beware - looks like someone edited a classification during this comparison!!! {start_last_modified=} vs {end_last_modified=}")
