from dataclasses import dataclass
from typing import Optional, List

from django.http import HttpRequest
from django.shortcuts import render

from classification.models import ImportedAlleleInfo
from genes.hgvs import HGVSConverterType, HGVSMatcher
from genes.models import TranscriptVersion, TranscriptParts
from library.django_utils import require_superuser
from library.utils import all_equal
from snpdb.models import GenomeBuild, VariantCoordinate
from django import forms


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
    error: Optional[str] = None
    error_stage: Optional[str] = None

    @property
    def explicit_variant_coordinate(self):
        if vc := self.variant_coordinate:
            return vc.explicit_reference()

    @property
    def error_str(self):
        if error := self.error:
            return f"{self.error_stage}: {error}"


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

    genome_build: Optional[GenomeBuild] = None
    hgvs_str: Optional[str] = None

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
                    iai.variant_coordinate = v.coordinate

        use_matchers = [(HGVSConverterType.BIOCOMMONS_HGVS, "biocommons"), (HGVSConverterType.PYHGVS, "pyhgvs")]
        for matcher_id, matcher_name in use_matchers:
            matcher = HGVSMatcher(genome_build, hgvs_converter_type=matcher_id)

            output = MatcherOutput(matcher_name=matcher_name)
            all_output.append(output)

            stage = "Basic Parsing"
            try:
                vcd = matcher.get_variant_tuple_used_transcript_kind_method_and_matches_reference(hgvs_str)
                output.variant_coordinate = vcd.variant_coordinate

                stage = "Getting Transcript"
                output.transcript = TranscriptVersion.get_transcript_id_and_version(vcd.transcript_accession)

                if output.variant_coordinate and output.transcript:
                    stage = "Resolving c.HGVS"
                    if variant_details := matcher.variant_coordinate_to_c_hgvs_variant(output.variant_coordinate,
                                                                                       str(output.transcript)):
                        output.hgvs = variant_details.format()

            except Exception as ex:
                output.error = str(ex)
                output.error_stage = stage

    return render(request, "classification/hgvs_resolution_tool.html", context)

