import uuid
from collections import defaultdict
from typing import Tuple

from django.conf import settings
from django.contrib import messages
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.forms.models import model_to_dict
from django.http.response import JsonResponse, Http404
from django.shortcuts import render, get_object_or_404
from django.utils.datastructures import OrderedSet
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST
from itertools import combinations

from django.views.generic import TemplateView
from global_login_required import login_not_required
from lazy import lazy

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.forms import EnsemblGeneAnnotationVersionForm
from annotation.models.models import AnnotationVersion, EnsemblGeneAnnotation, Citation, VariantAnnotation, MIMGene
from annotation.models.molecular_consequence_enums import MolecularConsequenceColors
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.forms import GeneListForm, NamedCustomGeneListForm, GeneForm, UserGeneListForm, CustomGeneListForm, \
    GeneSymbolForm
from genes.models import GeneInfo, CanonicalTranscriptCollection, GeneListCategory, \
    GeneList, GeneCoverageCollection, GeneCoverageCanonicalTranscript, \
    CustomTextGeneList, Transcript, Gene, TranscriptVersion, GeneSymbol, GeneCoverage, GeneVersion, \
    PfamSequenceIdentifier, GeneSymbolAliasSummary
from library.constants import MINUTE_SECS
from library.django_utils import get_field_counts, add_save_message
from library.utils import defaultdict_to_dict
from seqauto.models import EnrichmentKit
from snpdb.models import CohortGenotypeCollection, Cohort, VariantZygosityCountCollection
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.models.models_variant import Variant
from snpdb.variant_queries import get_has_classifications_q, get_has_variant_tags, get_variant_queryset_for_gene_symbol
from classification.enums import ShareLevel
from classification.models import ClassificationModification, Classification
from classification.views.classification_datatables import ClassificationDatatableConfig
from variantopedia.variant_column_utils import get_gene_annotation_column_data


def genes(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    av = AnnotationVersion.latest(genome_build)
    ensembl_gene_annotation_version = av.ensembl_gene_annotation_version
    version_form = EnsemblGeneAnnotationVersionForm(initial={'version': ensembl_gene_annotation_version})

    context = {"genome_build": genome_build,
               "gene_form": GeneForm(),
               "version_form": version_form}
    return render(request, 'genes/genes.html', context)


def _get_latest_gene_version_for_user_settings(source, gene_versions, user_settings: UserSettings) -> GeneVersion:
    """ returns latest GeneVersion using default_genome_build, falling back on other user settings builds """

    if not gene_versions.exists():
        raise ValueError(f"{source} does not have any GeneVersions!")

    genome_builds = user_settings.get_genome_builds()
    user_build_gene_versions = gene_versions.filter(genome_build__in=genome_builds).order_by("-version")
    if not user_build_gene_versions.exists():
        user_builds = ", ".join([str(gb) for gb in genome_builds])
        supported_builds = ", ".join(gene_versions.values_list("genome_build", flat=True).distinct())
        msg = f"{source} has no GeneVersions for user genome builds: {user_builds} - supported: {supported_builds}"
        raise ValueError(msg)

    gene_version = user_build_gene_versions.filter(genome_build=user_settings.default_genome_build).first()
    if not gene_version:
        gene_version = user_build_gene_versions.first()  # Try another build
    return gene_version


def view_gene(request, gene_id):
    gene = get_object_or_404(Gene, pk=gene_id)

    gene_versions_by_build = defaultdict(dict)

    versions = set()
    for tv in gene.geneversion_set.order_by("version"):
        genome_build_id = tv.genome_build.pk
        version = tv.version or 0  # 0 = unknown
        gene_versions_by_build[genome_build_id][version] = tv
        versions.add(version)

    genome_build_ids = sorted(gene_versions_by_build.keys())
    gene_versions = []
    for version in sorted(versions):
        data = [version]
        for genome_build_id in genome_build_ids:
            gv = gene_versions_by_build.get(genome_build_id, {}).get(version)
            data.append(gv)
        gene_versions.append(data)

    transcript_versions = defaultdict(OrderedSet)
    transcript_biotype = defaultdict(set)
    for tv in TranscriptVersion.objects.filter(gene_version__gene=gene).order_by("transcript_id", "version"):
        transcript_versions[tv.transcript_id].add(tv.version)
        if tv.biotype:
            transcript_biotype[tv.transcript_id].add(tv.biotype)
    transcripts = {}
    for transcript_id, versions in transcript_versions.items():
        biotype = ", ".join(transcript_biotype.get(transcript_id, []))
        transcripts[transcript_id] = (versions, biotype)

    context = {
        "gene": gene,
        "genome_build_ids": genome_build_ids,
        "gene_versions": gene_versions,
        "transcripts": transcripts,
    }
    return render(request, "genes/view_gene.html", context)


def _get_mim_and_hpo_for_gene_symbol(gene_symbol: GeneSymbol):
    genes = gene_symbol.get_genes()
    mim_set = set((mim_gene.mim_morbid for mim_gene in MIMGene.objects.filter(gene__in=genes)))
    mim_and_hpo_for_gene = []
    for mim in sorted(mim_set, key=lambda m: m.accession):
        hpo_set = set((pm.hpo for pm in mim.phenotypemim_set.order_by("pk")))
        hpo_list = sorted(hpo_set, key=lambda h: h.pk)
        mim_and_hpo_for_gene.append((mim, hpo_list))
    return mim_and_hpo_for_gene


def view_gene_symbol(request, gene_symbol):
    gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
    consortium_genes_and_aliases = defaultdict(lambda: defaultdict(set))
    for gene in gene_symbol.get_genes():
        aliases = consortium_genes_and_aliases[gene.get_annotation_consortium_display()][gene.identifier]
        aliases.update(gene.get_symbols().exclude(symbol=gene_symbol))

    user_settings = UserSettings.get_for_user(request.user)
    genome_build = user_settings.default_genome_build

    has_classified_variants = False
    has_observed_variants = False
    has_tagged_variants = False

    gene_versions = gene_symbol.geneversion_set.all()
    if gene_versions.exists():
        # This page is shown using the users default genome build
        # However - it's possible the gene doesn't exist for a particular genome build.
        # If gene has a version for a build in user settings, use that and show a warning message
        gene_version = _get_latest_gene_version_for_user_settings(f"GeneSymbol: {gene_symbol}", gene_versions, user_settings)
        genome_build = gene_version.genome_build
        if genome_build != user_settings.default_genome_build:
            messages.add_message(request, messages.WARNING,
                                 f"This symbol is not associated with any genes in your default build, viewing in build {genome_build}")

        gene_variant_qs = get_variant_queryset_for_gene_symbol(gene_symbol=gene_symbol, genome_build=genome_build, traverse_aliases=True)
        gene_variant_qs, count_column = VariantZygosityCountCollection.annotate_global_germline_counts(gene_variant_qs)
        has_observed_variants = gene_variant_qs.filter(**{f"{count_column}__gt": 0}).exists()

        q = get_has_variant_tags()
        has_tagged_variants = gene_variant_qs.filter(q).exists()

        q = get_has_classifications_q(genome_build)
        has_classified_variants = gene_variant_qs.filter(q).exists()
    else:
        messages.add_message(request, messages.WARNING, "There are no genes linked against this symbol!")

    has_variants = has_observed_variants or has_classified_variants or has_tagged_variants

    gene_annotation = EnsemblGeneAnnotation.get_for_symbol(genome_build, gene_symbol)
    if gene_annotation:
        gene_level_columns = get_gene_annotation_column_data(gene_annotation)
        ega_qs = gene_annotation.gene.ensemblgeneannotation_set.filter(version__genome_build=genome_build)
        num_gene_annotation_versions = ega_qs.count()
    else:
        gene_level_columns = None
        num_gene_annotation_versions = 0

    mim_and_hpo_for_gene = _get_mim_and_hpo_for_gene_symbol(gene_symbol)
    gene_lists_qs = GeneList.filter_for_user(request.user)
    gene_in_gene_lists = GeneList.visible_gene_lists_containing_gene_symbol(gene_lists_qs, gene_symbol).exists()

    citations = Citation.objects.filter(genesymbolcitation__gene_symbol=gene_symbol).distinct()

    has_gene_coverage = GeneCoverage.get_for_symbol(genome_build, gene_symbol).exists()
    has_canonical_gene_coverage = GeneCoverageCanonicalTranscript.get_for_symbol(genome_build, gene_symbol).exists()

    aliases_out = [alias for alias in gene_symbol.aliases_to_for if not alias.my_symbol_is_main]
    aliases_out.sort()

    aliases_in = [alias for alias in gene_symbol.aliases_to_for if alias.my_symbol_is_main]
    aliases_in.sort()

    filter_aliases = [gs.symbol for gs in gene_symbol.traverse_aliases()]
    filter_aliases.sort()

    context = {
        "consortium_genes_and_aliases": defaultdict_to_dict(consortium_genes_and_aliases),
        "citations": citations,
        "gene_symbol": gene_symbol,
        "gene_annotation": gene_annotation,
        "gene_in_gene_lists": gene_in_gene_lists,
        "gene_infos": GeneInfo.get_for_gene_symbol(gene_symbol),
        "gene_level_columns": gene_level_columns,
        "genome_build": genome_build,
        "has_classified_variants": has_classified_variants,
        "show_classifications_hotspot_graph": settings.VIEW_GENE_SHOW_CLASSIFICATIONS_HOTSPOT_GRAPH and has_classified_variants,
        "show_hotspot_graph": settings.VIEW_GENE_SHOW_HOTSPOT_GRAPH and has_observed_variants,
        "has_gene_coverage": has_gene_coverage or has_canonical_gene_coverage,
        "has_tagged_variants": has_tagged_variants,
        "has_variants": has_variants,
        "mim_and_hpo_for_gene": mim_and_hpo_for_gene,
        "num_gene_annotation_versions": num_gene_annotation_versions,
        "show_wiki": settings.VIEW_GENE_SHOW_WIKI,
        "show_annotation": settings.VARIANT_DETAILS_SHOW_ANNOTATION,
        "datatable_config": ClassificationDatatableConfig(request),
        "aliases_out": aliases_out,
        "aliases_in": aliases_in,
        "filter_aliases": ', '.join(filter_aliases)
    }
    return render(request, "genes/view_gene_symbol.html", context)


def view_gene_annotation_history(request, genome_build_name, gene_symbol):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)

    ega = EnsemblGeneAnnotation.get_for_symbol(genome_build, gene_symbol)
    if ega is None:
        raise Http404(f"No EnsemblGeneAnnotation for {genome_build}/{gene_symbol}")

    gene_annotation_dicts_by_version = {}
    versions = []
    for ega in ega.gene.ensemblgeneannotation_set.filter(version__genome_build=genome_build):
        ega_dict = model_to_dict(ega)
        del ega_dict["id"]
        del ega_dict["version"]
        gene_annotation_dicts_by_version[ega.version.pk] = ega_dict
        versions.append((ega.version.pk, str(ega.version)))

    context = {"gene": ega.gene,
               "gene_annotation_dicts_by_version": gene_annotation_dicts_by_version,
               "versions": versions}
    return render(request, "genes/view_gene_annotation_history.html", context)


def view_transcript(request, transcript_id):
    transcript = get_object_or_404(Transcript, pk=transcript_id)

    gene_by_build = defaultdict(set)
    transcripts_versions_by_build = defaultdict(dict)

    versions = set()
    for tv in transcript.transcriptversion_set.order_by("version"):
        genome_build_id = tv.genome_build.pk
        gene_by_build[genome_build_id].add(tv.gene)
        version = tv.version or 0  # 0 = unknown
        transcripts_versions_by_build[genome_build_id][version] = tv
        versions.add(version)

    genome_build_ids = sorted(gene_by_build.keys())
    build_genes = [gene_by_build.get(genome_build_id) for genome_build_id in genome_build_ids]
    transcript_versions = []
    for version in sorted(versions):
        data = [version]
        for genome_build_id in genome_build_ids:
            tv = transcripts_versions_by_build.get(genome_build_id, {}).get(version)
            data.append(tv)

        transcript_versions.append(data)

    context = {"transcript": transcript,
               "genome_build_ids": genome_build_ids,
               "build_genes": build_genes,
               "transcript_versions": transcript_versions}
    return render(request, "genes/view_transcript.html", context)


def view_transcript_version(request, transcript_id, version):
    transcript = Transcript.objects.filter(pk=transcript_id).first()
    if not transcript:
        return render(request, "genes/view_transcript.html", {'transcript_id': transcript_id})
    tv_set = transcript.transcriptversion_set.filter(version=version)
    tv: TranscriptVersion = tv_set.first()
    data = transcript.transcriptversion_set.aggregate(Count("version", distinct=True))
    version_count = data["version__count"]

    context = {"accession": TranscriptVersion.get_accession(transcript_id, version),
               "transcript": transcript,
               "version_count": version_count}

    if not tv:
        # https://www.theonion.com/area-man-constantly-mentioning-he-doesnt-own-a-televisi-1819565469
        return render(request, "genes/view_transcript_version.html", context)

    accession = tv.accession
    transcripts_by_build = {}
    for t in tv_set.order_by("genome_build__name"):
        genome_build_id = t.genome_build.pk
        transcripts_by_build[genome_build_id] = t

    differences = []
    for (a, b) in combinations(transcripts_by_build.keys(), 2):
        t_a = transcripts_by_build[a]
        t_b = transcripts_by_build[b]
        diff = t_a.get_differences(t_b)
        if diff:
            differences.append(((a, b), diff))

    context = {**context, **{"accession": accession,
                             "transcripts_by_build": transcripts_by_build,
                             "differences": differences}}
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
    unmatched_genes = list(gene_list.unmatched_genes)
    if unmatched_genes:
        if instructions:
            messages.add_message(request, messages.WARNING, instructions)

        for unmatched_gene in unmatched_genes:
            msg = f"Unmatched gene symbol: {unmatched_gene.original_name}"
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
        for (k, v) in data.items():
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


class HotspotGraphView(TemplateView):
    template_name = "genes/hotspot_graph.html"
    has_graph_filter_toolbar = True

    def _get_title(self, transcript_version) -> str:
        return f"Variants for {transcript_version}"

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

    @property
    def transcript(self):
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
                raise PfamSequenceIdentifier.DoesNotExist(f"No PfamSequenceIdentifier for {gene_id or gene_symbol_id}")
        else:
            raise ValueError("At least one of 'gene_symbol', 'gene_id' or 'transcript_id' must be in url kwargs")
        return transcript

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
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
        return Classification.objects.filter(variantclassificationmodification__in=vcm_qs)

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
    def transcript(self):
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
