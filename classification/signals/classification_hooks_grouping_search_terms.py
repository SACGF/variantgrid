from collections import defaultdict
from collections.abc import Iterable
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Optional

from django.dispatch import receiver

from annotation.models import AnnotationVersion, VariantAnnotation
from classification.models import (
    ClassificationGrouping,
    ClassificationGroupingSearchTermBuilder,
    ClassificationGroupingSearchTermStub,
    ClassificationGroupingSearchTermType,
    ClinVarExport,
    ImportedAlleleInfo,
    ResolvedVariantInfo,
    classification_grouping_search_term_signal,
)
from genes.models import GeneSymbol, GeneSymbolAlias, GeneVersion
from ontology.models import OntologyTerm
from snpdb.models import GenomeBuild

_latest_annotation_versions: ContextVar[Optional[dict[GenomeBuild, Optional[AnnotationVersion]]]] = \
    ContextVar("latest_annotation_versions", default=None)


@contextmanager
def latest_annotation_versions_cached():
    """ For rebuilding many groupings' search terms at once: each build's latest AnnotationVersion is looked up once """
    token = _latest_annotation_versions.set({})
    try:
        yield
    finally:
        _latest_annotation_versions.reset(token)


def _latest_annotation_version(genome_build: GenomeBuild, grouping: ClassificationGrouping) -> Optional[AnnotationVersion]:
    annotation_versions = _latest_annotation_versions.get()
    if annotation_versions is not None and genome_build in annotation_versions:
        return annotation_versions[genome_build]
    av = AnnotationVersion.latest_or_none(genome_build, context=f"ClassificationGrouping {grouping.pk} search terms")
    if annotation_versions is not None:
        annotation_versions[genome_build] = av
    return av


@receiver(classification_grouping_search_term_signal)
def _gene_symbol_search_for(grouping: ClassificationGrouping, **kwargs) -> Optional[Iterable[ClassificationGroupingSearchTermStub]]:
    imported_allele_infos = ImportedAlleleInfo.objects.filter(
        pk__in=grouping.classificationgroupingentry_set.all().values_list(
            'classification__allele_info_id', flat=True)
    ).select_related("imported_genome_build_patch_version__genome_build")

    all_transcript_version_ids: set[int] = set()
    genome_build_to_variant_ids: dict[GenomeBuild, set[int]] = defaultdict(set)
    gene_symbol_stubs: dict[str, ClassificationGroupingSearchTermBuilder] = {}

    def stub_for(gene_symbol: GeneSymbol | str):
        gene_symbol = str(gene_symbol).upper()
        if existing := gene_symbol_stubs.get(gene_symbol):
            return existing
        entry = ClassificationGroupingSearchTermBuilder(term=gene_symbol, term_type=ClassificationGroupingSearchTermType.GENE_SYMBOL)
        gene_symbol_stubs[gene_symbol] = entry
        return entry

    # GeneSymbol's primary key is the symbol
    for rb in ResolvedVariantInfo.objects.filter(allele_info__in=imported_allele_infos).select_related("genome_build"):
        genome_build_to_variant_ids[rb.genome_build].add(rb.variant_id)
        all_transcript_version_ids.add(rb.transcript_version_id)
        if gene_symbol_id := rb.gene_symbol_id:
            stub_for(gene_symbol_id).extra["from_normalized"] = True

    for imported_allele_info in imported_allele_infos:
        if c_hgvs := imported_allele_info.imported_c_hgvs_obj:
            if gene_symbol := GeneSymbol.objects.filter(symbol=c_hgvs.gene_symbol).first():
                stub_for(gene_symbol).extra["from_imported"] = True
            else:
                for alias in GeneSymbolAlias.objects.filter(alias=c_hgvs.gene_symbol).select_related("gene_symbol"):
                    stub_for(alias.gene_symbol).extra["from_imported"] = True

    # gene symbols from transcripts
    if all_transcript_version_ids:
        for gene_symbol_id in GeneVersion.objects.filter(transcriptversion__in=all_transcript_version_ids, gene_symbol__isnull=False).values_list(
                "gene_symbol_id", flat=True):
            stub_for(gene_symbol_id).extra["from_transcript"] = True

    # gene symbols from allele
    for genome_build, variant_ids in genome_build_to_variant_ids.items():
        av = _latest_annotation_version(genome_build, grouping)
        if av is None:
            continue
        vav = av.variant_annotation_version
        gene_ids = VariantAnnotation.objects.filter(variant__in=variant_ids, version=vav, gene__isnull=False).values_list("gene_id")
        for gene_symbol_id in GeneVersion.objects.filter(gene__in=gene_ids, gene_symbol__isnull=False).values_list("gene_symbol_id", flat=True).distinct():
            stub_for(gene_symbol_id).extra["from_allele"] = True

    return [term_builder.as_stub() for term_builder in gene_symbol_stubs.values()]


@receiver(classification_grouping_search_term_signal)
def _condition_terms(grouping: ClassificationGrouping, **kwargs) -> Optional[Iterable[ClassificationGroupingSearchTermStub]]:
    all_terms: set[OntologyTerm] = set()
    all_stubs: list[ClassificationGroupingSearchTermStub] = []

    # redundantly duplicate a lot for the condition free text
    # this means that this needs to stay in sync with ontology
    for modification in grouping.classification_modifications:
        if condition := modification.classification.condition_resolution_obj:
            all_terms |= set(condition.terms)

    for term in all_terms:
        all_stubs.append(ClassificationGroupingSearchTermStub(
            term_type=ClassificationGroupingSearchTermType.CONDITION_ID,
            term=term.id.upper()
        ))

    return all_stubs


# @receiver(classification_grouping_search_term_signal)
# def _patient_sample_ids(grouping: ClassificationGrouping, **kwargs) -> Optional[Iterable[ClassificationGroupingSearchTermStub]]:
#     # TODO consider only checking this if we've got
#     stubs = []
#     for cm in grouping.classification_modifications:
#         if patient_id := cm.get(SpecialEKeys.PATIENT_ID):
#             stubs.append(ClassificationGroupingSearchTermStub(
#                 term_type=ClassificationGroupingSearchTermType.PATIENT_ID,
#                 term=patient_id
#             ))
#         if sample_id := cm.get(SpecialEKeys.SAMPLE_ID):
#             stubs.append(ClassificationGroupingSearchTermStub(
#                 term_type=ClassificationGroupingSearchTermType.PATIENT_ID,
#                 term=sample_id
#             ))
#     return stubs


@receiver(classification_grouping_search_term_signal)
def _clinvar_scv(grouping: ClassificationGrouping, **kwargs) -> Optional[Iterable[ClassificationGroupingSearchTermStub]]:
    stubs = []
    for scv in ClinVarExport.objects.filter(classification_based_on__classification__in=[cm.classification_id for cm in grouping.classification_modifications]).values_list("scv", flat=True):
        if scv:
            stubs.append(ClassificationGroupingSearchTermStub(
                term_type=ClassificationGroupingSearchTermType.CLINVAR_SCV,
                term=scv
            ))
    return stubs
