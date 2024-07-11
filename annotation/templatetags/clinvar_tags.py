from typing import Optional, Union

import pydantic
from django.template import Library
from more_itertools import first
from pydantic import ConfigDict

from annotation.models import ClinVar, AnnotationVersion, VariantAnnotation
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from library.log_utils import report_exc_info
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Allele, Variant, GenomeBuild

register = Library()


@register.inclusion_tag("annotation/tags/clinvar_stars.html")
def clinvar_stars(stars, review_status: Optional[str] = None):
    MAX_STARS = 4
    star_bools = ([True] * stars) + ([False] * (MAX_STARS - stars))
    return {"stars": star_bools, "review_status": review_status}


class ClinVarDetails(pydantic.BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True, strict=True)

    """
    All the data that we have from ClinVar based on the pre-loaded annotations.
    If we don't have a ClinVar instance, then we have access to the g_hgvs so we can prompt the user to search
    """

    clinvar: Optional[ClinVar]
    is_desired_build: bool
    """
    When asking for ClinVar details on an allele, set this based on if it was the genome build that was originally
    requested for
    """

    genome_build: GenomeBuild
    annotation_version: AnnotationVersion
    g_hgvs: Optional[str]

    @property
    def is_expert_panel_or_greater(self) -> bool:
        return self.clinvar and self.clinvar.is_expert_panel_or_greater

    @staticmethod
    def instance_from(
            allele: Optional[Union[int, Allele]] = None,
            variant: Optional[Union[int, Variant]] = None,
            genome_build: Optional[GenomeBuild] = None,
            annotation_version: Optional[AnnotationVersion] = None) -> Optional['ClinVarDetails']:

        if not allele and not variant:
            raise ValueError("One of allele or variant must be provided")

        if not genome_build:
            if annotation_version:
                genome_build = annotation_version.genome_build
            else:
                genome_build = GenomeBuildManager.get_current_genome_build()

        is_desired_build = True

        if not variant:
            variant: Variant
            if isinstance(allele, int):
                allele = Allele.objects.get(pk=allele)
            variant = allele.variant_for_build_optional(genome_build)
            if not variant:
                if variant := allele.variants.first():
                    genome_build = first(variant.genome_builds)
                is_desired_build = False

        if not annotation_version:
            annotation_version = AnnotationVersion.latest(genome_build)

        clinvar_record: Optional[ClinVar] = None

        g_hgvs: Optional[str] = None
        if variant:
            if variant.can_have_annotation:
                clinvar_qs = ClinVar.objects.filter(variant=variant, version=annotation_version.clinvar_version)
                try:
                    clinvar_record = clinvar_qs.get()
                except ClinVar.MultipleObjectsReturned:
                    # Report this - but carry on for the user
                    report_exc_info({"target": f"Variant {variant.pk}), Annotation Version {annotation_version.pk}"})
                    clinvar_record = clinvar_qs.first()
                except ClinVar.DoesNotExist:
                    pass

            if not clinvar_record and variant.can_make_g_hgvs:
                g_hgvs = VariantAnnotation.get_hgvs_g(variant)

        return ClinVarDetails(
            clinvar=clinvar_record,
            is_desired_build=is_desired_build,
            genome_build=genome_build,
            annotation_version=annotation_version,
            g_hgvs=g_hgvs
        )


@register.inclusion_tag("annotation/tags/clinvar_tag.html")
def clinvar(
        allele: Optional[Union[int, Allele]] = None,
        variant: Optional[Union[int, Variant]] = None,
        genome_build: Optional[GenomeBuild] = None,
        annotation_version: Optional[AnnotationVersion] = None,
        expert_panel_only=False):
    """
    :param allele: In most cases can just provide allele and nothing else
    :param variant: If retrieving just for a single genome build (provide in combination with genome build)
    :param genome_build: Defaults to genome build manager
    :param annotation_version: Defaults to the latest for a genome build
    :param expert_panel_only: Set to true if we only care about expert panels
    """

    data = ClinVarDetails.instance_from(allele=allele, variant=variant, genome_build=genome_build, annotation_version=annotation_version)

    return {
        "data": data,
        "expert_panel_only": expert_panel_only,
        "min_stars": CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE if expert_panel_only else 0
    }
