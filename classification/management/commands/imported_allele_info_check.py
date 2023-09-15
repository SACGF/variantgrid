from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Any, List

from django.core.management import BaseCommand
from django.db.models import Max
from django.urls import reverse
from django.utils.http import urlencode

from classification.models import ImportedAlleleInfo
from genes.hgvs import HGVSMatcher, HGVSConverterType
from genes.models import TranscriptVersion, TranscriptParts
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row
from snpdb.models import Variant, GenomeBuild, VariantCoordinate


class Change(int, Enum):
    same = 0
    gained = 1
    lost = 2
    changed = 3
    changed_from_exact = 4
    changed_to_exact = 5

    @staticmethod
    def compare_3(provided: Any, resolved: Any, updated: Any) -> 'Change':
        if resolved == updated:
            return Change.same
        elif (not resolved) and updated:
            return Change.gained
        elif resolved and (not updated):
            return Change.lost
        else:
            if provided:
                if provided == resolved:
                    return Change.changed_from_exact
                if provided == updated:
                    return Change.changed_to_exact
            return Change.changed

    @staticmethod
    def compare_2(resolved: Any, updated: Any) -> 'Change':
        if resolved == updated:
            return Change.same
        elif (not resolved) and updated:
            return Change.gained
        elif resolved and (not updated):
            return Change.lost
        else:
            return Change.changed

    def __str__(self):
        if self == Change.same:
            return "same"
        elif self == Change.gained:
            return "gained value"
        elif self == Change.lost:
            return "lost value"
        elif self == Change.changed:
            return "changed"
        elif self == Change.changed_from_exact:
            return "changed - lost exact match"
        elif self == Change.changed_to_exact:
            return "changed - gained exact match"

    def __repr__(self):
        return str(self)


@dataclass
class HgvsSummary(ExportRow):
    transcript: Optional[TranscriptParts] = None
    variant_coordinate: Optional[VariantCoordinate] = None
    c_hgvs: Optional[str] = None
    errors: List[str] = None
    error_str: Optional[str] = None

    @property
    def variant_coordinate_str(self):
        if self.variant_coordinate:
            return Variant.format_tuple(*self.variant_coordinate, abbreviate=True)
        else:
            return None

    @property
    def c_hgvs_str(self):
        if self.c_hgvs:
            return str(self.c_hgvs)
        else:
            return None


@dataclass
class ChgvsDiff(ExportRow):

    imported_allele_info: ImportedAlleleInfo
    provided: HgvsSummary
    resolved: HgvsSummary
    updated: HgvsSummary

    @property
    def has_difference(self):
        return self.variant_coordinate_change() or self.transcript_change() or self.c_hgvs_change()

    # @export_column("$site_name URL")
    # def _allele_id(self):
    #     return get_url_from_view_path(self.imported_allele_info.get_absolute_url())

    @export_column("$site_name Compare URL")
    def _compare_url(self):
        return get_url_from_view_path(reverse('hgvs_resolution_tool')) + "?" + urlencode(
            {
                "genome_build": self.imported_allele_info.imported_genome_build.pk,
                "hgvs": self.imported_allele_info.imported_c_hgvs
            }
        )

    # @export_column("Labs")
    # def _labs(self):
    #     parts = []
    #     for classification in self.imported_allele_info.classification_set.select_related('lab', 'lab__organization'):
    #         parts.append(str(classification.lab))
    #     return ", ".join(parts)

    @export_column("Imported Genome Build")
    def _imported_genome_build(self):
        return self.imported_allele_info.imported_genome_build

    @export_column("c.HGVS Imported")
    def _c_hgvs_imported(self):
        return self.provided.c_hgvs_str

    @export_column("c.HGVS/VariantCoordinate/Transcript Change")
    def _all_changes(self):
        return str(self.c_hgvs_change()) + "/" + str(self.variant_coordinate_change()) + "/" + str(self.transcript_change())

    @export_column("Previous Error")
    def _error_previous(self):
        return self.resolved.error_str

    @export_column("New Error")
    def _error_new(self):
        return self.updated.error_str

    @export_column("c.HGVS Change")
    def c_hgvs_change(self):
        return Change.compare_3(self.provided.c_hgvs_str, self.resolved.c_hgvs_str, self.updated.c_hgvs_str)

    @export_column("c.HGVS Previous")
    def _c_hgvs_previous(self):
        return self.resolved.c_hgvs_str

    @export_column("c.HGVS New")
    def _c_hgvs_new(self):
        return self.updated.c_hgvs_str

    @export_column("Variant Coordinate Change")
    def variant_coordinate_change(self):
        genome_build = self._imported_genome_build
        resolved = self.resolved.variant_coordinate
        if resolved:
            resolved = resolved.as_internal_symbolic()

        updated = self.updated.variant_coordinate
        if updated:
            updated = updated.as_internal_symbolic()

        return Change.compare_2(resolved, updated)

    @export_column("Variant Coordinate Old")
    def _variant_coordinate_old(self):
        return self.resolved.variant_coordinate_str

    @export_column("Variant Coordinate New")
    def _variant_coordinate_new(self):
        return self.updated.variant_coordinate_str

    @export_column("Transcript Change")
    def transcript_change(self):
        return Change.compare_3(self.provided.transcript, self.resolved.transcript, self.updated.transcript)

    @export_column("Transcript Imported")
    def _transcript_imported(self):
        return self.provided.transcript

    @export_column("Transcript Previous")
    def _transcript_old(self):
        return self.resolved.transcript

    @export_column("Transcript New")
    def _transcript_updated(self):
        return self.updated.transcript


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
        c_hgvs_diff_count = Counter()

        iai: ImportedAlleleInfo
        hgvs_matchers_by_build = {}
        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build, hgvs_converter_type=HGVSConverterType.BIOCOMMONS_HGVS)
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
                    provided.c_hgvs = imported_c_hgvs
                    provided.transcript = matcher.get_transcript_parts(imported_c_hgvs)
                except Exception as ex:
                    provided.error_str = ex.__class__.__name__ + ": " + str(ex)

                # RESOLVED
                try:
                    if variant_info := iai[iai.imported_genome_build]:
                        resolved.c_hgvs = variant_info.c_hgvs
                        if variant := variant_info.variant:
                            resolved.variant_coordinate = variant.coordinate
                            if transcript_version := variant_info.transcript_version:
                                resolved.transcript = transcript_version.as_parts
                        if error := variant_info.error:
                            resolved.error_str = error
                except Exception as ex:
                    resolved.error_str = ex.__class__.__name__ + ": " + str(ex)

                # UPDATED
                stage = "Getting Coordinates"
                try:
                    vcd = matcher.get_variant_coordinate_used_transcript_kind_method_and_matches_reference(imported_c_hgvs)
                    updated.variant_coordinate = vcd.variant_coordinate

                    stage = "Getting Transcript"
                    updated.transcript = TranscriptVersion.get_transcript_id_and_version(vcd.transcript_accession)

                    if updated.variant_coordinate and updated.transcript:
                        stage = "Resolving c.HGVS"
                        if hgvs_variant := matcher.variant_coordinate_to_hgvs_variant(updated.variant_coordinate, str(updated.transcript)):
                            updated.c_hgvs = hgvs_variant.format()

                except Exception as ex:
                    updated.error_str = stage + ": " + ex.__class__.__name__ + ": " + str(ex)

                if diff.has_difference:
                    out.write(delimited_row(diff.to_csv()))

                variant_diff_count[diff.variant_coordinate_change()] += 1
                transcript_diff_count[diff.transcript_change()] += 1
                c_hgvs_diff_count[diff.c_hgvs_change()] += 1

        print("==== Variant matching ====")
        print(dict(variant_diff_count))
        print("==== Transcripts: ==== ")
        print(dict(transcript_diff_count))
        print("==== c.HGVS: ==== ")
        print(dict(c_hgvs_diff_count))

        end_last_modified = self._get_last_modified()
        if start_last_modified != end_last_modified:
            print(f"Beware - looks like someone edited a classification during this comparison!!! {start_last_modified=} vs {end_last_modified=}")
