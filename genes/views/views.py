from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations
from typing import Tuple, List, Optional, Dict, Set, Iterable, Union, Any

from django.conf import settings
from django.contrib import messages
from django.db.models import QuerySet
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.http.response import JsonResponse
from django.shortcuts import render, get_object_or_404
from django.urls import reverse
from django.utils.datastructures import OrderedSet
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST

from analysis.models import VariantTag
from annotation.models import Citation
from annotation.models.models import AnnotationVersion, DBNSFPGeneAnnotationVersion, DBNSFPGeneAnnotation
from classification.models import ClassificationModification
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.views.classification_datatables import ClassificationColumns
from classification.views.exports import ClassificationExportFormatterCSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_csv import FormatDetailsCSV
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.forms import GeneListForm, NamedCustomGeneListForm, UserGeneListForm, CustomGeneListForm, \
    GeneSymbolForm, GeneAnnotationReleaseGenomeBuildForm
from genes.hgvs import HGVSMatcher
from genes.models import GeneInfo, CanonicalTranscriptCollection, GeneListCategory, \
    GeneList, GeneCoverageCollection, GeneCoverageCanonicalTranscript, \
    CustomTextGeneList, Transcript, Gene, TranscriptVersion, GeneSymbol, GeneCoverage, \
    PanelAppServer, SampleGeneList, HGNC, GeneVersion, TranscriptVersionSequenceInfo, NoTranscript
from genes.models_enums import AnnotationConsortium
from genes.serializers import SampleGeneListSerializer
from library.constants import WEEK_SECS
from library.django_utils import get_field_counts, add_save_message
from library.utils import defaultdict_to_dict, LazyAttribute
from ontology.models import OntologySnake, OntologyService, OntologyTerm
from seqauto.models import EnrichmentKit
from snpdb.models import VariantZygosityCountCollection, Sample, VariantGridColumn
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
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


def _get_omim_and_hpo_for_gene_symbol(gene_symbol: GeneSymbol) -> List[Tuple[OntologyTerm, List[OntologyTerm]]]:
    omim_and_hpo_for_gene = []
    try:
        # max_depth = 0 for direct links only
        for omim in OntologySnake.terms_for_gene_symbol(gene_symbol, OntologyService.OMIM, max_depth=0).leafs():
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

    def __init__(self, gene_symbol: GeneSymbol, desired_genome_build: Optional[GenomeBuild], user,
                 tool_tips=None):
        self.gene_symbol = gene_symbol
        self.desired_genome_build = desired_genome_build
        self.tool_tips = tool_tips
        if user is not None:
            self.user = user
            if self.tool_tips is None:
                user_settings = UserSettings.get_for_user(user)
                self.tool_tips = user_settings.tool_tips

    @cached_property
    def omim_and_hpo_for_gene(self) -> List[Tuple[OntologyTerm, List[OntologyTerm]]]:
        return _get_omim_and_hpo_for_gene_symbol(self.gene_symbol)

    @cached_property
    def hgnc(self) -> Optional[HGNC]:
        return self.gene_symbol.hgnc_set.order_by("status").first()

    @cached_property
    def citations_ids(self) -> List[str]:
        return sorted(set(Citation.objects.filter(genesymbolcitation__gene_symbol=self.gene_symbol).values_list('id', flat=True)))

    @cached_property
    def dbnsfp_gene_annotation(self) -> Optional[DBNSFPGeneAnnotation]:
        dga = None
        if dbnsfp_gene_version := DBNSFPGeneAnnotationVersion.latest():
            dga = self.gene_symbol.dbnsfpgeneannotation_set.filter(version=dbnsfp_gene_version).first()
        return dga

    @cached_property
    def gene_summary(self) -> Optional[str]:
        refseq_gene: Gene
        if refseq_gene := Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ,
                                              geneversion__gene_symbol=self.gene_symbol).first():
            return refseq_gene.summary
        return None

    @cached_property
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

    @cached_property
    def genome_build(self) -> GenomeBuild:
        if gene_version := self.gene_version:
            return gene_version.genome_build
        return self.desired_genome_build

    def warnings(self) -> List[str]:
        warnings = []
        if self.gene_version:
            # This page is shown using the users default genome build
            # However - it's possible the gene doesn't exist for a particular genome build.
            # If gene has a version for a build in user settings, use that and show a warning message
            if self.genome_build != self.desired_genome_build:
                warnings.append(f"This symbol is not associated with any genes in build {self.desired_genome_build}, viewing in build {self.genome_build}")
        else:
            warnings.append("There are no genes linked against this symbol!")
        return warnings

    @cached_property
    def has_variants(self) -> HasVariants:

        has_tagged_variants = False
        has_observed_variants = False
        has_classified_variants = False

        if self.gene_version:
            annotation_version = AnnotationVersion.latest(self.genome_build)
            gene_variant_qs = get_variant_queryset_for_gene_symbol(self.gene_symbol, annotation_version,
                                                                   traverse_aliases=True)
            gene_variant_qs, vzcc = VariantZygosityCountCollection.annotate_global_germline_counts(gene_variant_qs)
            has_observed_variants = gene_variant_qs.filter(**{f"{vzcc.germline_counts_alias}__gt": 0}).exists()

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

    @cached_property
    def consortium_genes_and_aliases(self) -> Dict[str, Set[str]]:
        consortium_genes_and_aliases = defaultdict(lambda: defaultdict(set))
        gene: Gene
        for gene in self.gene_symbol.genes:
            aliases = consortium_genes_and_aliases[gene.get_annotation_consortium_display()][gene.identifier]
            aliases.update(gene.get_symbols().exclude(symbol=self.gene_symbol))
        return defaultdict_to_dict(consortium_genes_and_aliases)

    @cached_property
    def gene_external_urls(self) -> Dict[str, str]:
        gene_external_urls: Dict[str, str] = {}
        for gene in self.gene_symbol.genes:
            gene_external_urls[gene.identifier] = gene.get_external_url()
        return gene_external_urls

    @cached_property
    def annotation_description(self):
        descriptions = {}
        if self.tool_tips:
            descriptions = VariantGridColumn.get_column_descriptions()
            descriptions["gnomad_gene_constraint"] = """
            constraint score shown in gnomAD is the ratio of the observed / expected (oe) number of loss-of-function
            variants in that gene. The expected counts are based on a mutational model that takes sequence context,
            coverage and methylation into account. Low oe values are indicative of strong intolerance. Range is 90%
            confidence interval. <a href='http://gnomad-sg.org/help/constraint'>Details at gnomAD</a> """
            descriptions["essential_gene"] = f"""
                <p><b>CRISPR:</b> {descriptions['essential_gene_crispr']}</p>
                <p><b>CRISPR2:</b>{descriptions['essential_gene_crispr2']}</p>
                <p><b>Gene Trap:</b>{descriptions['essential_gene_gene_trap']}</p>
            """
        return descriptions

    @cached_property
    def has_gene_coverage(self) -> bool:
        if not settings.VIEW_GENE_SYMBOL_SHOW_GENE_COVERAGE:
            return False
        has_gene_coverage = GeneCoverage.get_for_symbol(self.genome_build, self.gene_symbol).exists()
        if has_gene_coverage:
            return True
        has_canonical_gene_coverage = GeneCoverageCanonicalTranscript.get_for_symbol(self.genome_build, self.gene_symbol).exists()
        return has_canonical_gene_coverage

    @property
    def has_samples_in_other_builds(self) -> bool:
        return Sample.objects.exclude(vcf__genome_build=self.genome_build).exists()

    @cached_property
    def gene_in_gene_lists(self) -> bool:
        gene_lists_qs = GeneList.filter_for_user(self.user)
        gene_in_gene_lists = GeneList.visible_gene_lists_containing_gene_symbol(gene_lists_qs, self.gene_symbol).exists()
        return gene_in_gene_lists

    @cached_property
    def gene_infos(self):
        return GeneInfo.get_for_gene_symbol(self.gene_symbol)

    @cached_property
    def classifications(self) -> QuerySet[ClassificationModification]:
        # Note this is loaded in Ajax
        classifications_qs = ClassificationModification.objects.none()
        if filters := classification_gene_symbol_filter(self.gene_symbol):
            classifications = ClassificationModification.objects.filter(filters).filter(
                is_last_published=True).exclude(classification__withdrawn=True)
            classifications_qs = ClassificationModification.filter_for_user(user=self.user, queryset=classifications)

        classifications_qs = classifications_qs.select_related('classification', 'classification__lab', 'classification__allele', 'classification__allele__clingen_allele')
        return classifications_qs

    def panel_app_servers(self) -> Union[QuerySet, Iterable[PanelAppServer]]:
        return PanelAppServer.objects.order_by("pk")

    def show_classifications_hotspot_graph(self) -> bool:
        return settings.VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS and self.has_variants.has_classified_variants

    def show_hotspot_graph(self) -> bool:
        return settings.VIEW_GENE_HOTSPOT_GRAPH and self.has_variants.has_observed_variants


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
        desired_genome_build=UserSettings.get_genome_build_or_default(request.user, genome_build_name),
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
            "gene_in_gene_lists",
            "gene_infos",
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
            "show_hotspot_graph"
        ]
    )
    context["show_wiki"] = settings.VIEW_GENE_WIKI
    context["show_annotation"] = settings.VARIANT_DETAILS_SHOW_ANNOTATION
    context["datatable_config"] = ClassificationColumns(request)

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
    genes: List[Gene]


@dataclass(frozen=True)
class TranscriptVersionDetails:
    tv: Optional[TranscriptVersion]
    version: int  # redundant to tv if provided
    genome_build: GenomeBuild  # redundant to tv if provided
    hgvs_method: Any


def view_transcript(request, transcript_id):
    transcript = get_object_or_404(Transcript, pk=transcript_id)

    gene_by_build: Dict[GenomeBuild, Set[Gene]] = defaultdict(set)
    transcripts_versions_by_build = defaultdict(dict)

    versions: Set[int] = set()
    transcript_chroms = set()
    for tv in transcript.transcriptversion_set.order_by("version"):
        gene_by_build[tv.genome_build].add(tv.gene)
        version = tv.version or 0  # 0 = unknown
        transcripts_versions_by_build[tv.genome_build][version] = tv
        versions.add(version)
        transcript_chroms.update(tv.get_chromosomes())

    genome_builds = sorted(gene_by_build.keys())
    genome_build_genes = [GenomeBuildGenes(genome_build, sorted(gene_by_build.get(genome_build))) for genome_build in genome_builds]
    transcript_version_details: List[TranscriptVersionDetails] = []

    build_matcher = {genome_build: HGVSMatcher(genome_build) for genome_build in genome_builds}
    for version in sorted(versions):
        transcript_accession = f"{transcript}.{version}"
        for genome_build in genome_builds:
            tv = transcripts_versions_by_build.get(genome_build, {}).get(version)
            matcher = build_matcher[genome_build]
            hgvs_method = matcher.filter_best_transcripts_and_method_by_accession(transcript_accession)

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
            filter_q = Q(gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__enrichment_kit=enrichment_kit)
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


def gene_wiki(request):
    return render(request, "genes/gene_wiki.html")
