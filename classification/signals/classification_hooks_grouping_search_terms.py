from collections import defaultdict
from typing import Iterable, Optional, Set

from django.dispatch import receiver

from annotation.models import AnnotationVersion, VariantAnnotation
from classification.models import ClassificationGroupingSearchTermStub, ClassificationGroupingSearchTermType, \
    ClassificationGrouping, ImportedAlleleInfo, classification_grouping_search_term_signal, \
    ClassificationGroupingSearchTermBuilder, ClinVarExport
from genes.models import GeneSymbol, GeneSymbolAlias, GeneVersion, TranscriptVersion
from ontology.models import OntologyTerm
from snpdb.models import GenomeBuild, Variant


@receiver(classification_grouping_search_term_signal)
def _gene_symbol_search_for(grouping: ClassificationGrouping, **kwargs) -> Optional[Iterable[ClassificationGroupingSearchTermStub]]:
    imported_allele_infos = ImportedAlleleInfo.objects.filter(
        pk__in=grouping.classificationgroupingentry_set.all().values_list(
            'classification__allele_info_id', flat=True)
    )

    all_transcript_versions: set[TranscriptVersion] = set()
    genome_build_to_variants: dict[GenomeBuild, Variant] = defaultdict(set)
    gene_symbol_stubs: dict[str, ClassificationGroupingSearchTermBuilder] = {}

    def stub_for(gene_symbol: GeneSymbol | str):
        gene_symbol = str(gene_symbol).upper()
        if existing := gene_symbol_stubs.get(gene_symbol):
            return existing
        entry = ClassificationGroupingSearchTermBuilder(term=gene_symbol, term_type=ClassificationGroupingSearchTermType.GENE_SYMBOL)
        gene_symbol_stubs[gene_symbol] = entry
        return entry

    for imported_allele_info in imported_allele_infos:
        for rb in imported_allele_info.resolved_builds:
            if variant := rb.variant:
                genome_build_to_variants[rb.genome_build].add(variant)

            all_transcript_versions.add(rb.transcript_version)
            if gene_symbol := rb.gene_symbol:
                stub_for(gene_symbol).extra["from_normalized"] = True

        if c_hgvs := imported_allele_info.imported_c_hgvs_obj:
            if gene_symbol := GeneSymbol.objects.filter(symbol=c_hgvs.gene_symbol).first():
                stub_for(gene_symbol).extra["from_imported"] = True
            else:
                for alias in GeneSymbolAlias.objects.filter(alias=c_hgvs.gene_symbol).select_related("gene_symbol"):
                    stub_for(alias.gene_symbol).extra["from_imported"] = True

    # gene symbols from transcripts
    if all_transcript_versions:
        for gv in GeneVersion.objects.filter(transcriptversion__in=all_transcript_versions).select_related(
                "gene_symbol"):
            stub_for(gv.gene_symbol).extra["from_transcript"] = True

    # gene symbols from allele
    for genome_build, variants in genome_build_to_variants.items():
        vav = AnnotationVersion.latest(genome_build=genome_build).variant_annotation_version
        for va in VariantAnnotation.objects.filter(variant__in=variants, version=vav).select_related("gene"):
            if gene := va.gene:
                for gene_symbol in gene.get_symbols():
                    stub_for(gene_symbol).extra["from_allele"] = True

    return [term_builder.as_stub() for term_builder in gene_symbol_stubs.values()]


@receiver(classification_grouping_search_term_signal)
def _condition_terms(grouping: ClassificationGrouping, **kwargs) -> Optional[Iterable[ClassificationGroupingSearchTermStub]]:
    all_terms: Set[OntologyTerm] = set()
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
