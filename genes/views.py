import logging
import uuid
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Tuple, List, Optional, Dict, Set, Iterable, Union

from django.conf import settings
from django.contrib import messages
from django.db.models import QuerySet
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.http.response import JsonResponse
from django.shortcuts import render, get_object_or_404
from django.urls import reverse
from django.utils.datastructures import OrderedSet
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST
from django.views.generic import TemplateView
from global_login_required import login_not_required
from lazy import lazy

from analysis.models import VariantTag
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models.models import AnnotationVersion, Citation, VariantAnnotation
from annotation.models.molecular_consequence_enums import MolecularConsequenceColors
from classification.enums import ShareLevel
from classification.models import ClassificationModification, Classification
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.views.classification_datatables import ClassificationColumns
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.forms import GeneListForm, NamedCustomGeneListForm, UserGeneListForm, CustomGeneListForm, \
    GeneSymbolForm, GeneAnnotationReleaseGenomeBuildForm
from genes.hgvs import HGVSMatcher
from genes.models import GeneInfo, CanonicalTranscriptCollection, GeneListCategory, \
    GeneList, GeneCoverageCollection, GeneCoverageCanonicalTranscript, \
    CustomTextGeneList, Transcript, Gene, TranscriptVersion, GeneSymbol, GeneCoverage, \
    PfamSequenceIdentifier, PanelAppServer, SampleGeneList, HGNC, GeneVersion, TranscriptVersionSequenceInfo, \
    NoTranscript
from genes.models_enums import AnnotationConsortium
from genes.serializers import SampleGeneListSerializer
from library.constants import MINUTE_SECS
from library.django_utils import get_field_counts, add_save_message
from library.utils import defaultdict_to_dict, LazyAttribute
from ontology.models import OntologySnake, OntologyService, OntologyTerm
from seqauto.models import EnrichmentKit
from snpdb.models import CohortGenotypeCollection, Cohort, VariantZygosityCountCollection, Sample
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.models.models_variant import Variant
from snpdb.variant_queries import get_has_classifications_q, get_variant_queryset_for_gene_symbol


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
    for tv in TranscriptVersion.objects.filter(gene_version__gene=gene).order_by("transcript_id", "version"):
        transcript_genome_build_ids.add(tv.genome_build_id)
        transcript_versions_by_id_and_build[tv.transcript_id][tv.genome_build_id].add(tv)
        transcript_biotype[tv.transcript_id].add(tv.biotype)

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


def _get_omim_and_hpo_for_gene_symbol(gene_symbol: GeneSymbol) -> List[Tuple[OntologyTerm, List[OntologyTerm]]]:
    omim_and_hpo_for_gene = []
    try:
        for omim in OntologySnake.terms_for_gene_symbol(gene_symbol, OntologyService.OMIM, max_depth=0).leafs():  # direct links only
            hpo_list = OntologySnake.snake_from(omim, OntologyService.HPO, max_depth=0).leafs()
            omim_and_hpo_for_gene.append((omim, hpo_list))
    except ValueError:  # in case we don't have this gene symbol available
        pass

    return omim_and_hpo_for_gene


@dataclass(frozen=True)
class HasVariants:
    has_tagged_variants: bool
    has_observed_variants: bool
    has_classified_variants: bool

    @property
    def has_variants(self) -> bool:
        return self.has_tagged_variants or self.has_observed_variants or self.has_classified_variants

    def __bool__(self):
        return self.has_variants


class GeneSymbolViewInfo:

    def __init__(self, gene_symbol: GeneSymbol, desired_genome_build: GenomeBuild, user):
        self.gene_symbol = gene_symbol
        self.desired_genome_build = desired_genome_build
        self.user = user

    @lazy
    def omim_and_hpo_for_gene(self) -> List[Tuple[OntologyTerm, List[OntologyTerm]]]:
        return _get_omim_and_hpo_for_gene_symbol(self.gene_symbol)

    @lazy
    def hgnc(self) -> Optional[HGNC]:
        return self.gene_symbol.hgnc_set.order_by("status").first()

    @lazy
    def citations(self) -> Union[QuerySet, Iterable[Citation]]:
        return Citation.objects.filter(genesymbolcitation__gene_symbol=self.gene_symbol).distinct()

    @lazy
    def gene_summary(self) -> Optional[str]:
        refseq_gene: Gene
        if refseq_gene := Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ,
                                              geneversion__gene_symbol=self.gene_symbol).first():
            return refseq_gene.summary
        return None

    @lazy
    def gene_version(self) -> Optional[GeneVersion]:
        gene_versions = self.gene_symbol.geneversion_set.all()
        if gene_versions.exists():
            # This page is shown using the users default genome build
            # However - it's possible the gene doesn't exist for a particular genome build.
            # If gene has a version for a build in user settings, use that and show a warning message
            gene_version = gene_versions.filter(genome_build=self.desired_genome_build).first()
            if not gene_version:
                gene_version = gene_versions.first()  # Try another build
            return gene_version
        return None

    @lazy
    def genome_build(self) -> GenomeBuild:
        if gene_version := self.gene_version:
            return gene_version.genome_build
        return self.desired_genome_build

    def warnings(self) -> List[str]:
        warnings = list()
        if self.gene_version:
            # This page is shown using the users default genome build
            # However - it's possible the gene doesn't exist for a particular genome build.
            # If gene has a version for a build in user settings, use that and show a warning message
            if self.genome_build != self.desired_genome_build:
                warnings.append(f"This symbol is not associated with any genes in build {self.desired_genome_build}, viewing in build {self.genome_build}")
        else:
            warnings.append("There are no genes linked against this symbol!")
        return warnings

    @lazy
    def has_variants(self) -> HasVariants:

        has_tagged_variants = False
        has_observed_variants = False
        has_classified_variants = False

        if self.gene_version:
            gene_variant_qs = get_variant_queryset_for_gene_symbol(gene_symbol=self.gene_symbol, genome_build=self.genome_build,
                                                                   traverse_aliases=True)
            gene_variant_qs, count_column = VariantZygosityCountCollection.annotate_global_germline_counts(
                gene_variant_qs)
            has_observed_variants = gene_variant_qs.filter(**{f"{count_column}__gt": 0}).exists()

            has_tagged_variants = VariantTag.get_for_build(self.genome_build, variant_qs=gene_variant_qs).exists()

            # has classifications isn't 100% in sync with the classification table: this code looks at VariantAlleles
            # wheras the classification table will filter on gene symbol and transcript evidence keys
            q = get_has_classifications_q(self.genome_build)
            has_classified_variants = gene_variant_qs.filter(q).exists()

        return HasVariants(
            has_tagged_variants=has_tagged_variants,
            has_observed_variants=has_observed_variants,
            has_classified_variants=has_classified_variants)

    @property
    def has_classified_variants(self):
        return self.has_variants.has_classified_variants

    @lazy
    def consortium_genes_and_aliases(self) -> Dict[str, Set[str]]:
        consortium_genes_and_aliases = defaultdict(lambda: defaultdict(set))
        gene: Gene
        for gene in self.gene_symbol.genes:
            aliases = consortium_genes_and_aliases[gene.get_annotation_consortium_display()][gene.identifier]
            aliases.update(gene.get_symbols().exclude(symbol=self.gene_symbol))
        return defaultdict_to_dict(consortium_genes_and_aliases)

    @lazy
    def gene_external_urls(self) -> Dict[str, str]:
        gene_external_urls: Dict[str, str] = dict()
        for gene in self.gene_symbol.genes:
            gene_external_urls[gene.identifier] = gene.get_external_url()
        return gene_external_urls

    @lazy
    def has_gene_coverage(self) -> bool:
        has_gene_coverage = GeneCoverage.get_for_symbol(self.genome_build, self.gene_symbol).exists()
        if has_gene_coverage:
            return True
        has_canonical_gene_coverage = GeneCoverageCanonicalTranscript.get_for_symbol(self.genome_build, self.gene_symbol).exists()
        return has_canonical_gene_coverage

    @lazy
    def gene_in_gene_lists(self) -> bool:
        gene_lists_qs = GeneList.filter_for_user(self.user)
        gene_in_gene_lists = GeneList.visible_gene_lists_containing_gene_symbol(gene_lists_qs, self.gene_symbol).exists()
        return gene_in_gene_lists

    @lazy
    def gene_infos(self):
        return GeneInfo.get_for_gene_symbol(self.gene_symbol)

    @lazy
    def classifications(self) -> List[ClassificationModification]:
        # Note this is loaded in Ajax
        classifications = list()
        if filters := classification_gene_symbol_filter(self.gene_symbol):
            classifications = ClassificationModification.objects.filter(filters).filter(
                is_last_published=True).exclude(classification__withdrawn=True)
            classifications = ClassificationModification.filter_for_user(user=self.user, queryset=classifications)
        return classifications

    def panel_app_servers(self) -> Union[QuerySet, Iterable[PanelAppServer]]:
        return PanelAppServer.objects.order_by("pk")

    def show_classifications_hotspot_graph(self) -> bool:
        return settings.VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH and self.has_variants.has_classified_variants

    def show_hotspot_graph(self) -> bool:
        return settings.VIEW_GENE_SHOW_HOTSPOT_GRAPH and self.has_variants.has_observed_variants


def view_gene_symbol(request, gene_symbol: str, genome_build_name: Optional[str] = None):

    view_info = GeneSymbolViewInfo(
        gene_symbol=get_object_or_404(GeneSymbol, pk=gene_symbol),
        desired_genome_build=UserSettings.get_genome_build_or_default(request.user, genome_build_name),
        user=request.user)

    for warning in view_info.warnings():
        messages.add_message(request, messages.WARNING, warning)

    context = LazyAttribute.lazy_context(
        view_info,
        [
            "consortium_genes_and_aliases",
            "citations",
            "gene_symbol",
            "gene_external_urls",
            "gene_in_gene_lists",
            "gene_infos",
            "gene_summary",
            "genome_build",
            "has_classified_variants",
            "has_gene_coverage",
            "hgnc",
            "panel_app_servers",
            "omim_and_hpo_for_gene",
            "has_variants",
            "show_classifications_hotspot_graph",
            "show_hotspot_graph"
        ]
    )
    context["show_wiki"] = settings.VIEW_GENE_SHOW_WIKI
    context["show_annotation"] = settings.VARIANT_DETAILS_SHOW_ANNOTATION
    context["datatable_config"] = ClassificationColumns(request)

    return render(request, "genes/view_gene_symbol.html", context)


def view_classifications(request, gene_symbol: str, genome_build_name: str):

    genome_build = GenomeBuild.get_from_fuzzy_string(genome_build_name)
    view_info = GeneSymbolViewInfo(
        gene_symbol=get_object_or_404(GeneSymbol, pk=gene_symbol),
        desired_genome_build=genome_build,
        user=request.user)

    return render(request, "genes/view_gene_symbol_classifications.html", {
        "classifications": view_info.classifications,
        "genome_build": genome_build
    })


def view_transcript(request, transcript_id):
    transcript = get_object_or_404(Transcript, pk=transcript_id)

    gene_by_build = defaultdict(set)
    transcripts_versions_by_build = defaultdict(dict)

    versions = set()
    for tv in transcript.transcriptversion_set.order_by("version"):
        gene_by_build[tv.genome_build].add(tv.gene)
        version = tv.version or 0  # 0 = unknown
        transcripts_versions_by_build[tv.genome_build][version] = tv
        versions.add(version)

    genome_builds = sorted(gene_by_build.keys())
    build_genes = [gene_by_build.get(genome_build) for genome_build in genome_builds]
    build_matcher = {genome_build: HGVSMatcher(genome_build) for genome_build in genome_builds}
    transcript_version_lengths = []
    transcript_version_hgvs_methods = []
    for version in sorted(versions):
        transcript_accession = f"{transcript}.{version}"
        length_row = [version]
        hgvs_row = [version]
        for genome_build in genome_builds:
            tv = transcripts_versions_by_build.get(genome_build, {}).get(version)
            length_row.append(tv)
            matcher = build_matcher[genome_build]
            hgvs_row.append(matcher.filter_best_transcripts_and_method_by_accession(transcript_accession))

        transcript_version_lengths.append(length_row)
        transcript_version_hgvs_methods.append(hgvs_row)

    context = {
        "transcript": transcript,
        "genome_builds": genome_builds,
        "build_genes": build_genes,
        "transcript_version_lengths": transcript_version_lengths,
        "transcript_version_hgvs_methods": transcript_version_hgvs_methods,
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
        other_chroms = False

        for tv in tv_set.order_by("genome_build__name"):
            if tv_sequence_info:
                poly_a_tail = max(poly_a_tail, tv.sequence_poly_a_tail)
            genome_build_id = tv.genome_build.pk
            alignment_gap |= tv.alignment_gap
            other_chroms |= bool(tv.data.get("other_chroms"))
            transcript_versions_by_build[genome_build_id] = tv
            if not tv.has_valid_data:
                builds_missing_data.add(tv.genome_build)

        if other_chroms:
            other_chroms_msg = f"This transcript maps to multiple coordinates in the genome."
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


def gene_lists_tab(request):
    context = {"categories": GeneListCategory.get_gene_list_categories()}
    return render(request, 'genes/gene_lists_tab.html', context=context)


def gene_lists(request):
    context = {"create_gene_list_form": NamedCustomGeneListForm(username=request.user),
               "gene_list_form": UserGeneListForm(),
               "categories": GeneListCategory.get_gene_list_categories()}
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
        valid = gl_form.is_valid()
        if valid:
            gene_list = gl_form.save()

        add_save_message(request, valid, "GeneList")

    add_gene_list_unmatched_genes_message(request, gene_list)

    context = {'gene_list': gene_list,
               'gene_list_form': gl_form,
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
        summary = get_field_counts(qs, "gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__enrichment_kit__name")

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


def gene_coverage_graphs(request, genome_build, gene_symbols):
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
            filter_q = Q(gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__enrichment_kit=enrichment_kit)
            enrichment_kit_data = get_coverage_stats(base_gene_coverage_qs, filter_q, fields)
            enrichment_kit_name = str(enrichment_kit)
            for field_name in fields:
                field_enrichment_kit_gene_json[field_name][enrichment_kit_name][gene_symbol.symbol] = enrichment_kit_data.get(field_name, [])

        if has_non_kit_coverage:
            filter_q = Q(gene_coverage_collection__qcgenecoverage__qc__isnull=True)
            other_data = get_coverage_stats(base_gene_coverage_qs, filter_q, fields)
            for field_name in fields:
                field_enrichment_kit_gene_json[field_name][NON_KIT_NAME][gene_symbol.symbol] = other_data.get(field_name, [])

    context = {'has_coverage': has_coverage,
               'fields': fields,
               'enrichment_kits_list': enrichment_kit_names,
               'gene_symbols': [gs.symbol for gs in gene_symbols],
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
    return gene_coverage_graphs(request, genome_build, [gene_symbol])


def qc_gene_list_coverage_graphs(request, genome_build_name, gene_list_id):
    gene_list = GeneList.get_for_user(request.user, gene_list_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    return gene_coverage_graphs(request, genome_build, gene_list.get_gene_names())


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
        sample_gene_lists_data.append(SampleGeneListSerializer(sgl).data)
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


class HotspotGraphView(TemplateView):
    template_name = "genes/hotspot_graph.html"
    has_graph_filter_toolbar = True

    def _get_title(self, transcript_version) -> str:
        return f"Het/Hom variants for {transcript_version}"

    def _get_y_label(self) -> str:
        return "# samples"

    def _get_variant_queryset(self, transcript_version):
        annotation_version = AnnotationVersion.latest(transcript_version.genome_build)
        qs = get_variant_queryset_for_annotation_version(annotation_version)
        return qs.filter(Variant.get_no_reference_q(),
                         varianttranscriptannotation__transcript_version=transcript_version,
                         varianttranscriptannotation__protein_position__isnull=False)

    def _get_values(self, transcript_version) -> Tuple[int, str, str, str, int]:
        """ :returns hgvs_p, pp, consequence, count, gnomad_af """
        qs = self._get_variant_queryset(transcript_version)
        qs, count_column = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        qs = qs.filter(**{count_column + "__gt": 0})
        return qs.values_list("varianttranscriptannotation__hgvs_p",
                              "varianttranscriptannotation__protein_position",
                              "varianttranscriptannotation__consequence",
                              "variantannotation__gnomad_af",
                              count_column)

    @lazy
    def genome_build(self):
        genome_build_name = self.kwargs["genome_build_name"]
        return GenomeBuild.get_name_or_alias(genome_build_name)

    @lazy
    def transcript(self) -> Optional[Transcript]:
        transcript_id = self.kwargs.get("transcript_id")
        gene_id = self.kwargs.get("gene_id")
        gene_symbol_id = self.kwargs.get("gene_symbol")
        transcript = None

        if transcript_id:
            transcript = Transcript.objects.get(pk=transcript_id)
        elif gene_id or gene_symbol_id:
            vav = self.genome_build.latest_variant_annotation_version
            if gene_id:
                gene = get_object_or_404(Gene, identifier=gene_id)
                transcript = PfamSequenceIdentifier.get_transcript_for_gene(gene, vav)
            elif gene_symbol_id:
                gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol_id)
                transcript = PfamSequenceIdentifier.get_transcript_for_gene_symbol(gene_symbol, vav)

            if transcript is None:
                logging.warning("No PfamSequenceIdentifier for %s", gene_id or gene_symbol_id)
        else:
            raise ValueError("At least one of 'gene_symbol', 'gene_id' or 'transcript_id' must be in url kwargs")
        return transcript

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if self.transcript:
            transcript_version = self.transcript.latest_version(self.genome_build)

            variant_data = []
            for hgvs_p, protein_position, consequence, gnomad_af, count in self._get_values(transcript_version):
                consequence = MolecularConsequenceColors.CONSEQUENCE_LOOKUPS.get(consequence,
                                                                                 MolecularConsequenceColors.OTHER)
                protein_aa1 = ""
                if hgvs_p:
                    protein_aa3 = hgvs_p.split(":p.")[1]
                    protein_aa1 = VariantAnnotation.amino_acid_3_to_1(protein_aa3)

                pp = VariantAnnotation.protein_position_to_int(protein_position)
                variant_data.append((protein_aa1, pp, consequence, gnomad_af, count))

            domains = transcript_version.get_domains()
            molecular_consequence_colors = dict(MolecularConsequenceColors.HOTSPOT_COLORS)
            context.update({"transcript_version": transcript_version,
                            "molecular_consequence_colors": molecular_consequence_colors,
                            "domains": list(domains),
                            "variant_data": variant_data,
                            "uuid": uuid.uuid4(),
                            "title": self._get_title(transcript_version),
                            "y_title": self._get_y_label(),
                            "has_graph_filter_toolbar": self.has_graph_filter_toolbar})
        return context


class ClassificationsHotspotGraphView(HotspotGraphView):
    has_graph_filter_toolbar = False

    def _get_title(self, transcript_version) -> str:
        return f"Classifications for {transcript_version}"

    def _get_y_label(self) -> str:
        return "# classifications"

    def _get_classifications(self):
        vcm_qs = ClassificationModification.latest_for_user(self.request.user, published=True)
        return Classification.objects.filter(classificationmodification__in=vcm_qs)

    def _get_values(self, transcript_version) -> Tuple[int, str, str, str, int]:
        """ :returns hgvs_p, pp, consequence, count, gnomad_af """
        qs = self._get_variant_queryset(transcript_version)
        vc_qs = self._get_classifications()
        # Need to join through Allele so we get classifications through all genome builds
        vc_column = "variantallele__allele__variantallele__variant__classification"
        qs = qs.filter(**{vc_column + "__in": vc_qs})
        count_column = "classifications_count"
        qs = qs.values("variantallele__allele").annotate(**{count_column: Count(vc_column)})
        qs = qs.filter(**{count_column + "__gt": 0})
        return qs.values_list("varianttranscriptannotation__hgvs_p",
                              "varianttranscriptannotation__protein_position",
                              "varianttranscriptannotation__consequence",
                              "variantannotation__gnomad_af",
                              count_column)


class CohortHotspotGraphView(HotspotGraphView):

    @lazy
    def genome_build(self):
        return self.cohort.genome_build

    @lazy
    def cohort(self):
        cohort_id = self.kwargs["cohort_id"]
        return Cohort.get_for_user(self.request.user, cohort_id)

    def _get_title(self, transcript_version) -> str:
        return f"{self.cohort} - {transcript_version}"

    def _get_values(self, transcript_version):
        """ :returns hgvs_p, pp, consequence, count, gnomad_af """
        qs = self._get_variant_queryset(transcript_version)
        qs = qs.filter(cohortgenotype__collection=self.cohort.cohort_genotype_collection)
        qs = qs.filter(variantannotation__protein_position__isnull=False)
        qs, count_column = CohortGenotypeCollection.annotate_all_counts(qs)
        qs = qs.filter(**{count_column + "__gt": 0})
        return qs.values_list("varianttranscriptannotation__hgvs_p",
                              "varianttranscriptannotation__protein_position",
                              "varianttranscriptannotation__consequence",
                              "variantannotation__gnomad_af",
                              count_column)


@method_decorator([cache_page(MINUTE_SECS), login_not_required], name='dispatch')
class PublicRUNX1HotspotGraphView(ClassificationsHotspotGraphView):
    """ RUNX1 would like a hotspot graph on the front page - but we don't want to expose
        every hotspot graph obviously, so make a special case one """

    @lazy
    def genome_build(self):
        return GenomeBuild.get_name_or_alias("GRCh37")  # Hardcoded for server

    @lazy
    def transcript(self) -> Optional[Transcript]:
        return Transcript.objects.get(pk="ENST00000300305")

    def _get_title(self, transcript_version) -> str:
        title = super()._get_title(transcript_version)
        return title + " (Germline)"

    def _get_classifications(self):
        # Modified from RUNX1_classified_damage - GeneCountType.get_classification_qs()
        kwargs = {"share_level__in": ShareLevel.same_and_higher(ShareLevel.ALL_USERS),
                  "published_evidence__allele_origin__value__in": ['germline', 'likely_germline'],
                  "is_last_published": True}
        vcm_qs = ClassificationModification.objects.filter(**kwargs)
        return Classification.objects.filter(pk__in=vcm_qs.values('classification'))
