from dataclasses import dataclass, asdict
from functools import cached_property
from typing import Optional, Union, List

from django.template import Library
from more_itertools import first

from annotation.clinvar_xml_parser import ClinVarParser
from annotation.models import ClinVar, AnnotationVersion
from genes.hgvs import HGVSMatcher
from library.log_utils import report_exc_info
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Allele, Variant, GenomeBuild

register = Library()


@register.inclusion_tag("annotation/tags/clinvar_stars.html")
def clinvar_stars(stars):
    MAX_STARS = 4
    stars = ([True] * stars) + ([False] * (MAX_STARS - stars))
    return {"stars": stars}


@dataclass
class ClinVarDetails:
    clinvar: ClinVar
    is_desired_build: bool
    genome_build: GenomeBuild
    annotation_version: AnnotationVersion
    g_hgvs: str
    clinvar_citations: Optional[List[str]]
    num_clinvar_citations: Optional[int]

    def instance_from(allele: Optional[Union[int, Allele]] = None, variant: Optional[Union[int, Variant]] = None, genome_build: Optional[GenomeBuild] = None, annotation_version: Optional[AnnotationVersion] = None):
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
                variant = allele.variants.first()
                genome_build = first(variant.genome_builds)
                is_desired_build = False

        if not annotation_version:
            annotation_version = AnnotationVersion.latest(genome_build)

        clinvar_record: Optional[ClinVar] = None

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

        clinvar_citations: List[str] = None
        num_clinvar_citations: Optional[int] = None

        g_hgvs: Optional[str] = None
        if not clinvar_record:
            g_hgvs = HGVSMatcher(genome_build).variant_to_g_hgvs(variant)
        else:
            clinvar_citations = clinvar_record.citation_ids
            num_clinvar_citations = len(clinvar_citations)

        return ClinVarDetails(
            clinvar=clinvar_record,
            is_desired_build=is_desired_build,
            genome_build=genome_build,
            annotation_version=annotation_version,
            g_hgvs=g_hgvs,
            clinvar_citations=clinvar_citations,
            num_clinvar_citations=num_clinvar_citations
        )


@register.inclusion_tag("annotation/tags/clinvar_tag.html")
def clinvar(allele: Optional[Union[int, Allele]] = None, variant: Optional[Union[int, Variant]] = None, genome_build: Optional[GenomeBuild] = None, annotation_version: Optional[AnnotationVersion] = None):
    data = ClinVarDetails.instance_from(allele=allele, variant=variant, genome_build=genome_build, annotation_version=annotation_version)
    return dict(asdict(data).items())