from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Any, Optional

from django.conf import settings
from django.contrib import messages
from django.db.models.aggregates import Count
from django.http.response import JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils.datastructures import OrderedSet
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST

from annotation.models.models import AnnotationVersion
from classification.views.exports import ClassificationExportFormatterCSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_csv import FormatDetailsCSV
from genes.forms import (
    GeneAnnotationReleaseGenomeBuildForm,
    GeneListForm,
    GeneSymbolForm,
    NamedCustomGeneListForm,
    UserGeneListForm,
)
from genes.gene_symbol_view_info import GeneSymbolViewInfo
from genes.graphs.gene_list_chromosome_graph import GeneListChromosomeGraph
from genes.hgvs import HGVSMatcher
from genes.models import (
    CanonicalTranscriptCollection,
    Gene,
    GeneList,
    GeneListCategory,
    GeneSymbol,
    NoTranscript,
    SampleGeneList,
    Transcript,
    TranscriptVersion,
    TranscriptVersionSequenceInfo,
)
from genes.serializers import SampleGeneListSerializer
from library.constants import WEEK_SECS
from library.django_utils import add_save_message, get_field_counts
from library.utils import LazyAttribute, full_class_name
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.graphs import graphcache
from snpdb.models import Sample
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings


def genes(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    av = AnnotationVersion.latest(genome_build)
    gene_annotation_release = av.gene_annotation_version.gene_annotation_release

    gene_annotation_release_form = GeneAnnotationReleaseGenomeBuildForm(genome_build=genome_build,
                                                                        initial={'release': gene_annotation_release})

    context = {"genome_build": genome_build,
               "gene_symbol_form": GeneSymbolForm(),
               "gene_annotation_release_form": gene_annotation_release_form}
    return render(request, 'genes/genes.html', context)


def view_gene(request, gene_id):
    gene = get_object_or_404(Gene, pk=gene_id)

    gene_versions_by_build = defaultdict(dict)

    versions = set()
    for gv in gene.geneversion_set.order_by("version"):
        genome_build_id = gv.genome_build.pk
        gene_versions_by_build[genome_build_id][gv.version] = gv
        versions.add(gv.version)

    gene_genome_build_ids = sorted(gene_versions_by_build.keys())
    gene_versions = []
    for version in sorted(versions):
        data = [version]
        for genome_build_id in gene_genome_build_ids:
            gv = gene_versions_by_build.get(genome_build_id, {}).get(version)
            data.append(gv)
        gene_versions.append(data)

    transcript_versions_by_id_and_build = defaultdict(lambda: defaultdict(OrderedSet))
    transcript_genome_build_ids = set()
    transcript_biotype = defaultdict(set)
    transcript_chroms = set()
    for tv in TranscriptVersion.objects.filter(gene_version__gene=gene).order_by("transcript_id", "version"):
        transcript_genome_build_ids.add(tv.genome_build_id)
        transcript_versions_by_id_and_build[tv.transcript_id][tv.genome_build_id].add(tv)
        transcript_biotype[tv.transcript_id].add(tv.biotype)
        transcript_chroms.update(tv.get_chromosomes())

    if len(transcript_chroms) > 1:
        other_chroms_msg = f"Gene has transcripts that maps to multiple chromosomes ({', '.join(transcript_chroms)})"
        messages.add_message(request, messages.WARNING, other_chroms_msg)

    transcript_genome_build_ids = sorted(transcript_genome_build_ids)
    transcript_versions = []  # ID, biotype, [[GRCh37 versions...], [GRCh38 versions]]
    has_biotype = False
    for transcript_id, biotype_set in transcript_biotype.items():
        biotype = ", ".join(sorted(b for b in biotype_set if b is not None))
        if biotype:
            has_biotype = True
        build_versions = []
        for genome_build_id in transcript_genome_build_ids:
            tv_set = transcript_versions_by_id_and_build[transcript_id].get(genome_build_id, [])
            build_versions.append(tv_set)
        transcript_versions.append([transcript_id, biotype, build_versions])

    context = {
        "gene": gene,
        "gene_genome_build_ids": gene_genome_build_ids,
        "gene_versions": gene_versions,
        "has_biotype": has_biotype,
        "transcript_versions": transcript_versions,
        "transcript_genome_build_ids": transcript_genome_build_ids,
    }
    return render(request, "genes/view_gene.html", context)



def export_classifications_gene_symbol(request, gene_symbol: str, genome_build_name: str):
    genome_build = GenomeBuild.get_from_fuzzy_string(genome_build_name)
    gene_symbol_info = GeneSymbolViewInfo(
        gene_symbol=get_object_or_404(GeneSymbol, pk=gene_symbol),
        desired_genome_build=genome_build,
        user=request.user
    )
    return ClassificationExportFormatterCSV(
        ClassificationFilter(
            user=request.user,
            genome_build=genome_build,
            file_prefix=f"classifications_{gene_symbol}",
            starting_query=gene_symbol_info.classifications
        ),
        FormatDetailsCSV()
    ).serve()


def view_gene_symbol(request, gene_symbol: str, genome_build_name: Optional[str] = None):

    view_info = GeneSymbolViewInfo(
        gene_symbol=get_object_or_404(GeneSymbol, pk=gene_symbol),
        desired_genome_build=GenomeBuildManager.get_current_genome_build(),
        user=request.user)

    for warning in view_info.warnings():
        messages.add_message(request, messages.WARNING, warning)

    context = LazyAttribute.lazy_context(
        view_info,
        [
            "annotation_description",
            "consortium_genes_and_aliases",
            "citations_ids",
            "dbnsfp_gene_annotation",
            "gene_symbol",
            "gene_external_urls",
            "gene_constraint",
            "gene_in_gene_lists",
            "gene_summary",
            "genome_build",
            "has_classified_variants",
            "has_gene_coverage",
            "has_samples_in_other_builds",
            "hgnc",
            "panel_app_servers",
            "omim_and_hpo_for_gene",
            "has_variants",
            "show_classifications_hotspot_graph",
            "show_hotspot_graph",
            "unmatched_classifications",
            "unmatched_classifications_title"
        ]
    )
    context["show_wiki"] = settings.VIEW_GENE_WIKI
    context["show_annotation"] = settings.VARIANT_DETAILS_SHOW_ANNOTATION

    return render(request, "genes/view_gene_symbol.html", context)


def view_classifications(request, gene_symbol: str, genome_build_name: str):

    genome_build = GenomeBuild.get_from_fuzzy_string(genome_build_name)
    gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)

    view_info = GeneSymbolViewInfo(
        gene_symbol=gene_symbol,
        desired_genome_build=genome_build,
        user=request.user)

    return render(request, "genes/view_gene_symbol_classifications.html", {
        "classifications": view_info.classifications,
        "gene_symbol": gene_symbol,
        "genome_build": genome_build
    })


@dataclass(frozen=True)
class GenomeBuildGenes:
    genome_build: GenomeBuild
    genes: list[Gene]


@dataclass(frozen=True)
class TranscriptVersionDetails:
    tv: Optional[TranscriptVersion]
    version: int  # redundant to tv if provided
    genome_build: GenomeBuild  # redundant to tv if provided
    hgvs_method: Any


def view_transcript(request, transcript_id):
    transcript = get_object_or_404(Transcript, pk=transcript_id)

    gene_by_build: dict[GenomeBuild, set[Gene]] = defaultdict(set)
    transcripts_versions_by_build = defaultdict(dict)

    versions: set[int] = set()
    transcript_chroms = set()
    for tv in transcript.transcriptversion_set.order_by("version"):
        gene_by_build[tv.genome_build].add(tv.gene)
        version = tv.version or 0  # 0 = unknown
        transcripts_versions_by_build[tv.genome_build][version] = tv
        versions.add(version)
        transcript_chroms.update(tv.get_chromosomes())

    genome_builds = sorted(gene_by_build.keys())
    genome_build_genes = [GenomeBuildGenes(genome_build, sorted(gene_by_build.get(genome_build))) for genome_build in genome_builds]
    transcript_version_details: list[TranscriptVersionDetails] = []

    for version in sorted(versions):
        transcript_accession = f"{transcript}.{version}"
        for genome_build in genome_builds:
            tv = transcripts_versions_by_build.get(genome_build, {}).get(version)
            matcher = HGVSMatcher.instance(genome_build)
            hgvs_method = matcher.filter_best_transcripts_and_converter_type_by_accession(transcript_accession)

            transcript_version_details.append(
                TranscriptVersionDetails(
                    tv=tv,
                    genome_build=genome_build,
                    version=version,
                    hgvs_method=hgvs_method
                )
            )

    if len(transcript_chroms) > 1:
        other_chroms_msg = f"Transcript maps to multiple chromosomes ({', '.join(transcript_chroms)})"
        messages.add_message(request, messages.WARNING, other_chroms_msg)

    context = {
        "transcript": transcript,
        "genome_build_genes": genome_build_genes,
        "transcript_version_details": transcript_version_details
    }
    return render(request, "genes/view_transcript.html", context)


def view_transcript_version(request, transcript_id, version):
    transcript = Transcript.objects.filter(pk=transcript_id).first()
    if not transcript:
        return render(request, "genes/view_transcript.html", {'transcript_id': transcript_id})

    accession = TranscriptVersion.get_accession(transcript_id, version)
    no_transcript_message = ""
    try:
        # Call this before retrieving TranscriptVersions - as it will retrieve it and set alignment_gap
        # if lengths are different
        tv_sequence_info = TranscriptVersionSequenceInfo.get(accession)
    except NoTranscript as e:
        tv_sequence_info = None
        no_transcript_message = str(e)

    tv_set = transcript.transcriptversion_set.filter(version=version)
    tv: TranscriptVersion = tv_set.first()
    data = transcript.transcriptversion_set.aggregate(Count("version", distinct=True))
    version_count = data["version__count"]

    context = {"accession": accession,
               "transcript": transcript,
               "tv_sequence_info": tv_sequence_info,
               "no_transcript_message": no_transcript_message,
               "version_count": version_count}

    poly_a_tail = 0
    if tv:
        transcript_versions_by_build = {}
        builds_missing_data = set()
        alignment_gap = False
        transcript_chroms = set()

        for tv in tv_set.order_by("genome_build__name"):
            if tv_sequence_info:
                poly_a_tail = max(poly_a_tail, tv.sequence_poly_a_tail)
            genome_build_id = tv.genome_build.pk
            alignment_gap |= tv.alignment_gap
            transcript_chroms.update(tv.get_chromosomes())
            transcript_versions_by_build[genome_build_id] = tv
            if not tv.has_valid_data:
                builds_missing_data.add(tv.genome_build)

        if len(transcript_chroms) > 1:
            other_chroms_msg = f"Transcript version maps to multiple chromosomes ({', '.join(transcript_chroms)})"
            messages.add_message(request, messages.WARNING, other_chroms_msg)

        differences = []
        if builds_missing_data:
            builds = ', '.join([str(b) for b in builds_missing_data])
            msg = f"Transcripts in builds {builds} missing data, no difference comparison possible"
            messages.add_message(request, messages.WARNING, msg)
        else:
            for a, b in combinations(transcript_versions_by_build.keys(), 2):
                t_a = transcript_versions_by_build[a]
                t_b = transcript_versions_by_build[b]
                diff = t_a.get_differences(t_b)
                if diff:
                    differences.append(((a, b), diff))

        context = {**context, **{"accession": accession,
                                 "transcript_versions_by_build": transcript_versions_by_build,
                                 "differences": differences,
                                 "alignment_gap": alignment_gap}}

    if tv_sequence_info:
        if poly_a_tail:
            sequence_length = f"{tv_sequence_info.length - poly_a_tail} (w/{poly_a_tail}bp polyA tail)"
        else:
            sequence_length = tv_sequence_info.length
        context["sequence_length"] = sequence_length
    return render(request, "genes/view_transcript_version.html", context)


def view_transcript_accession(request, transcript_accession):
    """ When you don't know whether it has a version or not """

    transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
    if version:
        return view_transcript_version(request, transcript_id, version)
    return view_transcript(request, transcript_id)


@cache_page(WEEK_SECS)
def gene_symbol_info_tab(request, gene_symbol, tool_tips):
    """ Condensed info about gene symbol, loaded into variants page """

    # We pass tool tips and don't use user/build so that it can be globally cached
    tool_tips = tool_tips == "True"
    view_info = GeneSymbolViewInfo(
        gene_symbol=get_object_or_404(GeneSymbol, pk=gene_symbol),
        desired_genome_build=None,
        user=None,
        tool_tips=tool_tips)

    context = LazyAttribute.lazy_context(
        view_info,
        [
            "annotation_description",
            "gene_summary",
            "gene_symbol",
            "hgnc",
            "panel_app_servers",
        ]
    )
    return render(request, 'genes/gene_symbol_info_tab.html', context=context)


def gene_lists(request):
    context = {
        "create_gene_list_form": NamedCustomGeneListForm(username=request.user),
        "gene_list_form": UserGeneListForm(),
    }
    return render(request, 'genes/gene_lists.html', context=context)


def add_gene_list_unmatched_genes_message(request, gene_list, instructions=None):
    if unmatched_symbols := list(gene_list.unmatched_gene_symbols):
        if instructions:
            messages.add_message(request, messages.WARNING, instructions)

        for glg in unmatched_symbols:
            msg = f"Unmatched gene symbol: {glg.original_name}"
            messages.add_message(request, messages.WARNING, msg)


def view_gene_list(request, gene_list_id):
    gene_list = GeneList.get_for_user(request.user, gene_list_id, success_only=False)
    gl_form = GeneListForm(request.POST or None, instance=gene_list)
    if request.method == "POST":
        gene_list.check_can_write(request.user)
        valid = gl_form.is_valid()
        if valid:
            gene_list = gl_form.save()

        add_save_message(request, valid, "GeneList")

    add_gene_list_unmatched_genes_message(request, gene_list)

    gene_symbols = sorted(gene_list.get_gene_names())
    context = {'gene_list': gene_list,
               'gene_list_form': gl_form,
               'gene_symbols_text': "\n".join(gene_symbols),
               'num_gene_symbols': len(gene_symbols),
               'has_write_permission': gene_list.can_write(request.user)}
    return render(request, 'genes/view_gene_list.html', context)


@require_POST
def create_custom_gene_list(request):
    form = NamedCustomGeneListForm(request.POST, username=request.user)
    if form.is_valid():
        custom_text_gene_list = form.save()
        gene_list_id = custom_text_gene_list.gene_list.pk
    else:
        gene_list_id = None

    return JsonResponse({"gene_list_id": gene_list_id})


def view_canonical_transcript_collection(request, pk):
    canonical_transcript_collection = get_object_or_404(CanonicalTranscriptCollection, pk=pk)

    summary = None
    qs = canonical_transcript_collection.genecoveragecanonicaltranscript_set.all()
    if qs.exists():
        summary = get_field_counts(qs, "gene_coverage_collection__qcgenecoverage__qc__bam_file__sequencing_sample__enrichment_kit__name")

    is_system_default = pk == str(settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID)
    context = {"canonical_transcript_collection": canonical_transcript_collection,
               "summary": summary,
               "is_system_default": is_system_default}
    return render(request, 'genes/view_canonical_transcript_collection.html', context)


def gene_grid(request, columns_from_url=None):
    # All the work is done in genes.templatetags.gene_grid_tags
    context = {"columns_from_url": columns_from_url}
    return render(request, 'genes/gene_grid.html', context)


def canonical_transcripts(request):

    default_ctc = None
    default_ctc_id = settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID
    if default_ctc_id:
        try:
            default_ctc = CanonicalTranscriptCollection.objects.get(pk=default_ctc_id)
        except CanonicalTranscriptCollection.DoesNotExist:
            message = "Could not load CanonicalTranscriptCollection for " \
                "settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID={default_ctc_id}"
            messages.add_message(request, messages.ERROR, message)
    else:
        message = "Setting 'GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID' not set."
        messages.add_message(request, messages.WARNING, message)

    context = {"default_ctc": default_ctc}
    return render(request, 'genes/canonical_transcripts.html', context)



def sample_gene_lists_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    has_write_permission = sample.can_write(request.user)

    create_gene_list_form = NamedCustomGeneListForm(request.POST or None, username=request.user,
                                                    initial={"name": "Sample Gene List"})
    if request.method == "POST":
        if create_gene_list_form.is_valid():
            custom_text_gene_list = create_gene_list_form.save()
            category = GeneListCategory.get_or_create_category(GeneListCategory.SAMPLE_GENE_LIST, hidden=True)
            gene_list = custom_text_gene_list.gene_list
            gene_list.category = category
            gene_list.save()
            SampleGeneList.objects.create(sample=sample, gene_list=gene_list)

    sample_gene_lists_data = []
    gene_grid_arg_list = []
    for sgl in sample.samplegenelist_set.all():
        sgl_serializer = SampleGeneListSerializer(sgl)
        sgl_serializer.context["request"] = request
        sample_gene_lists_data.append(sgl_serializer.data)
        gene_grid_arg_list.append(f"gene-list-{sgl.gene_list_id}")

    if gene_grid_arg_list:
        gene_grid_url = reverse("passed_gene_grid", kwargs={"columns_from_url": "/".join(gene_grid_arg_list)})
    else:
        gene_grid_url = None

    context = {"sample": sample,
               "has_write_permission": has_write_permission,
               "gene_grid_url": gene_grid_url,
               "create_gene_list_form": create_gene_list_form,
               "sample_gene_lists_data": sample_gene_lists_data}
    return render(request, 'genes/sample_gene_lists_tab.html', context)


def gene_wiki(request):
    return render(request, "genes/gene_wiki.html")


def gene_list_graphs_tab(request, gene_list_id):
    gene_list = GeneList.get_for_user(request.user, gene_list_id, success_only=False)
    context = {'gene_list': gene_list}
    return render(request, 'genes/gene_list_graphs_tab.html', context)


def gene_list_chromosome_graph(request, gene_list_id):
    GeneList.get_for_user(request.user, gene_list_id)  # permission check
    graph_class_name = full_class_name(GeneListChromosomeGraph)
    cached_graph = graphcache.async_graph(graph_class_name, gene_list_id)
    return redirect(cached_graph)
