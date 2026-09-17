from collections import defaultdict
from collections.abc import Iterable

from django.conf import settings
from django.db.models.query_utils import Q
from django.shortcuts import get_object_or_404, render

from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.forms import CustomGeneListForm, GeneSymbolForm, UserGeneListForm
from genes.models import (
    CustomTextGeneList,
    GeneCoverageCanonicalTranscript,
    GeneCoverageCollection,
    GeneList,
    GeneListCategory,
    GeneSymbol,
)
from seqauto.models import EnrichmentKit
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings


def gene_coverage_graphs(request, genome_build, gene_symbols: Iterable[str]):
    NON_KIT_NAME = "Other / No enrichment kit"
    fields = ("mean", "percent_20x")
    enrichment_kits = EnrichmentKit.get_enrichment_kits(settings.SEQAUTO_COVERAGE_ENRICHMENT_KITS)
    # field_enrichment_kit_gene_json = { field_name : {'enrichment_kit' : {"gene1" : list, "gene2" : list}} }
    field_enrichment_kit_gene_json = defaultdict(lambda: defaultdict(dict))
    enrichment_kit_names = list(map(str, enrichment_kits))
    has_non_kit_coverage = GeneCoverageCollection.objects.filter(qcgenecoverage__isnull=True).exists()
    has_coverage = has_non_kit_coverage

    if has_non_kit_coverage:
        enrichment_kit_names.append(NON_KIT_NAME)

    for gene_symbol in gene_symbols:
        base_gene_coverage_qs = GeneCoverageCanonicalTranscript.get_for_symbol(genome_build, gene_symbol)
        has_coverage = has_coverage or base_gene_coverage_qs.exists()

        for enrichment_kit in enrichment_kits:
            filter_q = Q(gene_coverage_collection__qcgenecoverage__qc__bam_file__sequencing_sample__enrichment_kit=enrichment_kit)
            enrichment_kit_data = get_coverage_stats(base_gene_coverage_qs, filter_q, fields)
            enrichment_kit_name = str(enrichment_kit)
            for field_name in fields:
                field_enrichment_kit_gene_json[field_name][enrichment_kit_name][gene_symbol] = enrichment_kit_data.get(field_name, [])

        if has_non_kit_coverage:
            filter_q = Q(gene_coverage_collection__qcgenecoverage__qc__isnull=True)
            other_data = get_coverage_stats(base_gene_coverage_qs, filter_q, fields)
            for field_name in fields:
                field_enrichment_kit_gene_json[field_name][NON_KIT_NAME][gene_symbol] = other_data.get(field_name, [])

    context = {'has_coverage': has_coverage,
               'fields': fields,
               'enrichment_kits_list': enrichment_kit_names,
               'gene_symbols': gene_symbols,
               'field_enrichment_kit_gene': field_enrichment_kit_gene_json}
    return render(request, 'genes/coverage/gene_coverage_graphs.html', context)


def get_coverage_stats(base_gene_coverage_qs, filter_q, fields):
    gene_coverage_qs = base_gene_coverage_qs.filter(filter_q)

    values = defaultdict(list)
    for data in gene_coverage_qs.values(*fields):
        for k, v in data.items():
            values[k].append(v)
    return values


def qc_coverage(request, genome_build_name=None):
    SPECIAL_COVERAGE_CUSTOM_GENE_LIST = f"__QC_COVERAGE_CUSTOM_GENE_LIST__{request.user}"
    custom_text_gene_list, _ = CustomTextGeneList.objects.get_or_create(name=SPECIAL_COVERAGE_CUSTOM_GENE_LIST)

    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    custom_gene_list_form = CustomGeneListForm(request.POST or None,
                                               initial={"custom_gene_list_text": custom_text_gene_list.text})
    if custom_gene_list_form.is_valid():
        custom_text_gene_list.text = custom_gene_list_form.cleaned_data['custom_gene_list_text']
        custom_text_gene_list.save()
        create_custom_text_gene_list(custom_text_gene_list, request.user, GeneListCategory.QC_COVERAGE_CUSTOM_TEXT, hidden=True)
        gene_list_id = custom_text_gene_list.gene_list.pk
    else:
        gene_list_id = None

    context = {"genome_build": genome_build,
               'gene_symbol_form': GeneSymbolForm(),
               'gene_list_id': gene_list_id,
               'gene_list_form': UserGeneListForm(),
               'custom_gene_list_form': custom_gene_list_form}
    return render(request, 'genes/coverage/qc_coverage.html', context)


def gene_coverage_collection_graphs(request, genome_build_name, gene_symbol):
    gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    return gene_coverage_graphs(request, genome_build, [gene_symbol.symbol])


def qc_gene_list_coverage_graphs(request, genome_build_name, gene_list_id):
    gene_list = GeneList.get_for_user(request.user, gene_list_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    gene_symbols = list(gene_list.get_gene_names())
    return gene_coverage_graphs(request, genome_build, gene_symbols)
