from dataclasses import dataclass
from typing import Optional, List

from django import forms
from django.http import HttpRequest
from django.shortcuts import render

from classification.models import ImportedAlleleInfo
from genes.hgvs import HGVSConverterType, HGVSMatcher
from genes.models import TranscriptVersion, TranscriptParts
from library.django_utils import require_superuser
from library.utils import all_equal
from snpdb.models import GenomeBuild, VariantCoordinate


class HgvsResolutionForm(forms.Form):
    genome_build = forms.ChoiceField(widget=forms.RadioSelect)
    hgvs = forms.CharField(label="HGVS")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['genome_build'].choices = GenomeBuild.get_choices()


@dataclass
class MatcherOutput:
    matcher_name: str
    imported_allele_info: Optional[ImportedAlleleInfo] = None
    variant_coordinate: Optional[VariantCoordinate] = None
    transcript_version: Optional[TranscriptParts] = None
    hgvs: Optional[str] = None
    message: Optional[str] = None

    @property
    def explicit_variant_coordinate(self):
        if vc := self.variant_coordinate:
            return vc.as_internal_symbolic()

    @property
    def is_error(self):
        return not self.hgvs


class MatcherOutputs:

    def __init__(self, genome_build: GenomeBuild, hgvs: str):
        self.all_output: List[MatcherOutput] = []
        self.genome_build = genome_build
        self.hgvs = hgvs

    def append(self, output: MatcherOutput):
        self.all_output.append(output)

    def __iter__(self):
        return iter(self.all_output)

    def __len__(self):
        return len(self.all_output)

    def all_equal_variant_coordinate(self) -> bool:
        return all_equal(x.explicit_variant_coordinate for x in self.all_output)

    def all_equal_transcript_version(self) -> bool:
        return all_equal(x.transcript_version for x in self.all_output)

    def all_equal_hgvs(self) -> bool:
        return all_equal(x.hgvs for x in self.all_output)


@require_superuser
def hgvs_resolution_tool(request: HttpRequest):

    hgvs_form = HgvsResolutionForm(request.GET)

    context = {"form": hgvs_form}

    if hgvs_form.is_valid():
        data = hgvs_form.clean()
        genome_build = GenomeBuild.get_name_or_alias(data.get("genome_build"))
        hgvs_str = data.get("hgvs")

        all_output = MatcherOutputs(genome_build=genome_build, hgvs=hgvs_str)
        context["outputs"] = all_output

        # check for existing
        for iai in ImportedAlleleInfo.objects.filter(imported_c_hgvs=hgvs_str, imported_genome_build_patch_version__genome_build=genome_build):
            output = MatcherOutput(matcher_name="Previously resolved", imported_allele_info=iai)
            all_output.append(output)
            if resolved_variant := iai[genome_build]:
                output.hgvs = resolved_variant.c_hgvs
                if tv := resolved_variant.transcript_version:
                    output.transcript_version = tv.as_parts
                if v := resolved_variant.variant:
                    output.variant_coordinate = v.coordinate
            output.message = iai.message

        for matcher_id in [HGVSConverterType.PYHGVS, HGVSConverterType.BIOCOMMONS_HGVS]:
            matcher = HGVSMatcher(genome_build, hgvs_converter_type=matcher_id)

            output = MatcherOutput(matcher_name=matcher.hgvs_converter.description())
            all_output.append(output)

            try:
                variant_coordinate: Optional[VariantCoordinate] = None
                if vcd := matcher.get_variant_coordinate_used_transcript_kind_method_and_matches_reference(hgvs_str):
                    variant_coordinate = vcd.variant_coordinate
                    output.variant_coordinate = variant_coordinate

                if vcd.transcript_accession:
                    output.transcript_version = TranscriptVersion.transcript_parts(vcd.transcript_accession)

                    if variant_coordinate:
                        if variant_details := matcher.variant_coordinate_to_hgvs_variant(variant_coordinate,
                                                                                         vcd.transcript_accession):
                            output.hgvs = variant_details.format()

            except Exception as ex:
                output.message = str(ex)

    return render(request, "classification/hgvs_resolution_tool.html", context)
