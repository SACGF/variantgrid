from collections import defaultdict
from typing import Iterable, Optional, Set
from django.dispatch import receiver
from annotation.models import AnnotationVersion, VariantAnnotation
from classification.enums import SpecialEKeys
from classification.models import ClassificationGroupingSearchTermStub, ClassificationGroupingSearchTermType, \
    ClassificationGrouping, ImportedAlleleInfo, classification_grouping_search_term_signal, \
    ClassificationGroupingSearchTermBuilder
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
    all_free_text_conditions: Set[str] = set()

    all_stubs: list[ClassificationGroupingSearchTermStub] = []

    # this is duplicated from the code that populates the grouped classification, but still seems cleanest to do so
    for modification in grouping.classification_modifications:
        if condition := modification.classification.condition_resolution_obj:
            all_terms |= set(condition.terms)
            if plain_text := condition.plain_text:
                all_free_text_conditions.add(plain_text.upper())
        elif condition_text := modification.get(SpecialEKeys.CONDITION):
            all_free_text_conditions.add(condition_text.upper())

    for term in all_terms:
        all_stubs.append(ClassificationGroupingSearchTermStub(
            term_type=ClassificationGroupingSearchTermType.CONDITION_ID,
            term=term.id.upper()
        ))
    for free_text in all_free_text_conditions:
        all_stubs.append(ClassificationGroupingSearchTermStub(
            term_type=ClassificationGroupingSearchTermType.CONDITION_TEXT,
            term=free_text
        ))
    return all_stubs