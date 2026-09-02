import csv
import json
from collections import defaultdict
from datetime import timedelta

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.forms import model_to_dict
from django.http import StreamingHttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils.text import slugify
from django.utils.timezone import localtime

from analysis.models import VariantTag
from annotation.models import (
    AnnotationVersion,
    Classification,
    ClassificationModification,
    ClinVarRecordCollection,
    VariantAnnotation,
    VariantTranscriptAnnotation,
)
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from eventlog.models import create_event
from genes.models import CanonicalTranscriptCollection, GeneSymbol
from library.django_utils import get_field_counts
from library.django_utils.grid_export import EXPORT_ROWS_PER_CHUNK
from library.git import Git
from library.log_utils import log_traceback
from library.utils import StashFile
from library.utils.date_utils import local_date_string
from pathtests.models import cases_for_user
from seqauto.models import VCFFromSequencingRun, get_20x_gene_coverage
from seqauto.seqauto_stats import get_sample_enrichment_kits_df
from snpdb.forms import TagForm, UserSelectForm, get_settings_form_features
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Sample, Tag, Variant, VariantGridColumn, get_igv_data
from snpdb.models.models_user_settings import UserSettings
from snpdb.search import search_data
from snpdb.serializers import VariantAlleleSerializer
from snpdb.utils import get_genome_build_or_404
from snpdb.variant_filters import (
    get_all_variant_types,
    get_all_variants_filters,
    get_gene_symbol_alias_strs,
    get_variant_type_label,
    resolve_gene_symbols,
)
from variantopedia import forms
from variantopedia.grids import VariantTagsColumns
from variantopedia.interesting_nearby import (
    get_method_summaries,
    get_nearby_qs,
    get_nearby_summaries,
)


def variants(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    initial_filters = get_all_variants_filters(request.user, genome_build)

    selected_gene_symbols = resolve_gene_symbols(initial_filters.get("gene_symbols"))
    # Aliases are traversed when querying, so tell the user which extra symbols that brought in
    gene_symbol_aliases = {}
    for gene_symbol in selected_gene_symbols:
        alias_symbol_strs = get_gene_symbol_alias_strs(gene_symbol)
        if extra_symbols := [s for s in alias_symbol_strs if s != gene_symbol.symbol]:
            gene_symbol_aliases[gene_symbol.symbol] = extra_symbols

    gene_symbol_form = forms.AllVariantsGeneSymbolForm(initial={"gene_symbols": selected_gene_symbols})
    context = {
        "genome_build": genome_build,
        "standard_contigs": genome_build.standard_contigs,
        "variant_types": [(vt, get_variant_type_label(vt)) for vt in get_all_variant_types()],
        "initial_filters": initial_filters,
        "initial_filters_json": json.dumps(initial_filters),
        "gene_symbol_form": gene_symbol_form,
        "gene_symbol_aliases": gene_symbol_aliases,
    }
    return render(request, "variantopedia/variants.html", context)


def dashboard(request):
    sample_enrichment_kits_df = None
    latest_sequencing_vcfs = []
    if settings.SEQAUTO_ENABLED:
        NUM_SEQUENCING_PROJECTS = 5
        qs = VCFFromSequencingRun.objects.all().order_by("-sequencing_run__name")
        latest_sequencing_vcfs = qs[:NUM_SEQUENCING_PROJECTS]
        sample_enrichment_kits_df = get_sample_enrichment_kits_df()

    user_has_cases = cases_for_user(request.user).exists()

    context = {'git': Git(settings.BASE_DIR),
               'user_has_cases': user_has_cases,
               'sample_enrichment_kits_df': sample_enrichment_kits_df,
               "latest_sequencing_vcfs": latest_sequencing_vcfs}
    return render(request, "variantopedia/dashboard.html", context)


def variant_tag_detail(request, variant_id, tag):
    """ Loaded via tags grid on variant page """

    variant = get_object_or_404(Variant, pk=variant_id)
    tag = get_object_or_404(Tag, pk=tag)
    if not VariantTag.filter_for_user(request.user).filter(variant=variant, tag=tag).exists():
        raise PermissionDenied
    context = {
        "variant": variant,
        "tag": tag,
    }
    return render(request, "variantopedia/variant_tag_detail.html", context)


VARIANT_GRID_ROW_DETAIL_MAX_TRANSCRIPTS = 20
VARIANT_GRID_ROW_DETAIL_MAX_CLASSIFICATIONS = 10


def variant_grid_row_detail(request, variant_id: int, annotation_version_id: int):
    """ The expanded row under a variant grid row - identifiers for the grid's annotation version.
        @see variantGridRowDetail in grid.js """
    variant = get_object_or_404(Variant, pk=variant_id)
    annotation_version = get_object_or_404(AnnotationVersion, pk=annotation_version_id)
    vav = annotation_version.variant_annotation_version
    variant_annotation = VariantAnnotation.objects.filter(variant=variant, version=vav).first()
    transcript_annotations = VariantTranscriptAnnotation.objects.filter(variant=variant, version=vav) \
        .exclude(hgvs_c__isnull=True).order_by("hgvs_c").select_related("transcript_version")
    # The representative transcript first, then the rest by accession
    representative_transcript_version_id = variant_annotation.transcript_version_id if variant_annotation else None
    transcript_annotations = sorted(transcript_annotations[:VARIANT_GRID_ROW_DETAIL_MAX_TRANSCRIPTS],
                                    key=lambda ta: (ta.transcript_version_id != representative_transcript_version_id,
                                                    ta.hgvs_c))
    # The records behind the Classifications column's internal chips
    classifications = []
    if allele := variant.allele:
        cm_qs = ClassificationModification.latest_for_user(request.user) \
            .filter(classification__allele=allele) \
            .select_related("classification", "classification__lab").order_by("classification__pk")
        classifications = list(cm_qs[:VARIANT_GRID_ROW_DETAIL_MAX_CLASSIFICATIONS])

    # Symbolic variants span genes rather than sitting in one, so name what they hit
    overlapping_symbols = None
    if variant.is_symbolic and variant_annotation:
        overlapping_symbols = variant_annotation.overlapping_symbols

    context = {
        "variant": variant,
        "variant_annotation": variant_annotation,
        "transcript_annotations": transcript_annotations,
        "build_variants": variant.all_build_variants,
        "clingen_allele": variant.allele.clingen_allele if variant.allele else None,
        "classifications": classifications,
        "overlapping_symbols": overlapping_symbols,
    }
    return render(request, "variantopedia/variant_grid_row_detail.html", context)


def view_variant(request, variant_id, genome_build_name=None):
    """ This is to open it with the normal menu around it (ie via search etc) """
    template = 'variantopedia/view_variant.html'

    variant = get_object_or_404(Variant, pk=variant_id)
    in_multiple_genome_builds = len(variant.genome_builds) > 1
    genome_build = None
    if genome_build_name:
        genome_build = get_genome_build_or_404(genome_build_name)
        if genome_build not in variant.genome_builds:
            return redirect(reverse('view_variant', kwargs={'variant_id': variant_id}))
    else:
        if in_multiple_genome_builds:
            user_settings = UserSettings.get_for_user(request.user)
            if user_settings.default_genome_build in variant.genome_builds:
                genome_build = user_settings.default_genome_build

    if genome_build is None:
        genome_build = variant.any_genome_build

    GenomeBuildManager.set_current_genome_build(genome_build)

    igv_data = get_igv_data(request.user, genome_build=genome_build)
    extra_context = {"igv_data": igv_data,
                     "in_multiple_genome_builds": in_multiple_genome_builds,
                     "view_variant": True,
                     "edit_clinical_groupings": request.GET.get('edit_clinical_groupings') == 'True'}

    annotation_version = AnnotationVersion.latest(genome_build)
    return variant_details_annotation_version(request, variant_id, annotation_version.pk,
                                              template=template, extra_context=extra_context)


def view_variant_annotation_history(request, variant_id):
    variant = get_object_or_404(Variant, pk=variant_id)

    annotation_versions = []
    variant_annotation_by_version = {}
    for va in variant.variantannotation_set.order_by("version"):
        annotation_versions.append((va.version.pk, str(va.version)))

        va_dict = model_to_dict(va, exclude=["id", "variant_id"])
        variant_annotation_by_version[va.version.pk] = va_dict

    context = {"variant": variant,
               "annotation_versions": annotation_versions,
               "annotation_by_version": variant_annotation_by_version}
    return render(request, "variantopedia/view_variant_annotation_history.html", context)


def variant_tags(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    variant_tags_qs = VariantTag.get_for_build(genome_build)
    tag_counts = sorted(get_field_counts(variant_tags_qs, "tag").items())
    month_ago = localtime() - timedelta(days=30)

    filter_user = None
    if (user_id := request.GET.get("user")) and user_id.isdigit():
        filter_user = User.objects.filter(pk=user_id).first()

    # The grids below show coordinates, so a tag needs a variant in this build to appear in them.
    # get_for_build is what they both start from, so ask it rather than re-deriving the rule - a tag
    # made in this build has its own variant and shows up before liftover assigns it an allele
    without_coordinate_qs = VariantTag.objects.exclude(
        pk__in=variant_tags_qs.values_list("pk", flat=True))
    if filter_user:
        without_coordinate_qs = without_coordinate_qs.filter(user=filter_user)

    context = {"genome_build": genome_build,
               "tag_counts": tag_counts,
               "tag_events": sum(count for _, count in tag_counts),
               "tag_events_last_month": variant_tags_qs.filter(created__gte=month_ago).count(),
               "tags_without_build_coordinate": without_coordinate_qs.count(),
               # Tags to start the grids filtered on - variants must carry all of them
               "initial_tags": request.GET.getlist("tag"),
               # User to start the grids filtered on - @see user page tags card
               "filter_user": filter_user,
               "user_form": UserSelectForm()}
    return render(request, 'variantopedia/variant_tags.html', context)


def search(request):
    search_string = None

    if request.GET.get('classify') is not None:
        form = forms.SearchAndClassifyForm(request.GET)
    else:
        form = forms.SearchForm(request.GET)

    user_settings = UserSettings.get_for_user(request.user)

    preview_mode = False
    classify = False
    if form.is_valid() and form.cleaned_data['search']:
        search_string = form.cleaned_data['search']
        classify = form.cleaned_data.get('classify')
        if mode := form.cleaned_data.get('mode'):
            preview_mode = mode == "preview"

    # always perform a "search" so we can get told what kind of searches are enabled
    # note that searching on "" doesn't actually invoke any of the other search logic

    search_results = search_data(user=request.user, search_string=search_string, classify=classify)

    #results, _search_types, _search_errors = search_results.results, search_results.search_types, search_results.search_errors
    details = search_results.summary
    create_event(request.user, 'search', details=details)

    single_preferred_result = search_results.single_preferred_result()
    if not preview_mode and single_preferred_result:
        return redirect(single_preferred_result.preview.internal_url)

    # Attempt to give hints on why nothing was found
    # for search_error, genome_builds in search_results.search_errors.items():
    #     text = f"{search_error.search_type}: {search_error.error}"
    #     if genome_builds:
    #         genome_builds_str = ", ".join(gb.name for gb in sorted(genome_builds))
    #         text += f" ({genome_builds_str})"
    #     messages.add_message(request, search_error.log_level, text)
    #
    # epk_qs = ExternalPK.objects.values_list("external_type", flat=True)
    # external_codes = list(sorted(epk_qs.distinct()))

    context = {
        "user_settings": user_settings,
        "form": form,
        "classify": classify,
        "search": search_string,
        "search_results": search_results,
        "single_preferred_result": single_preferred_result
        # "external_codes": external_codes,
    }
    return render(request, "variantopedia/search.html", context)


def variant_wiki(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    context = {
        "genome_build": genome_build,
    }
    return render(request, "variantopedia/variant_wiki.html", context)


def view_allele_from_variant(request, variant_id):
    variant = get_object_or_404(Variant, pk=variant_id)
    allele = variant.allele
    if allele and settings.PREFER_ALLELE_LINKS:
        return redirect(reverse('view_allele', kwargs={"allele_id": allele.id}))
    return redirect(reverse('view_variant', kwargs={"variant_id": variant_id}))


def get_genes_canonical_transcripts(variant, annotation_version):
    """ returns dict of list of (enrichment kit description, original_transcript) """

    vav = annotation_version.variant_annotation_version
    vst_anno = variant.varianttranscriptannotation_set.filter(version=vav)
    default_canonical_transcript_collection = CanonicalTranscriptCollection.get_default()

    genes_canonical_transcripts = defaultdict(list)
    for gene_symbol in GeneSymbol.objects.filter(pk__in=vst_anno.values_list("symbol")):
        for ct in gene_symbol.canonicaltranscript_set.filter(collection__enrichmentkit__isnull=False):
            for enrichment_kit in ct.collection.enrichmentkit_set.order_by("name", "version"):
                description = str(enrichment_kit)
                if ct.collection == default_canonical_transcript_collection:
                    description += " (default)"
                genes_canonical_transcripts[gene_symbol].append((description, ct.original_transcript))

    return dict(genes_canonical_transcripts)


def variant_details_annotation_version(request, variant_id, annotation_version_id,
                                       template='variantopedia/variant_details.html',
                                       extra_context: dict = None):
    """ Main Variant Details page """
    variant = get_object_or_404(Variant, pk=variant_id)
    annotation_version = get_object_or_404(AnnotationVersion, pk=annotation_version_id)
    genome_build = annotation_version.genome_build
    latest_annotation_version = AnnotationVersion.latest(genome_build)
    variant_annotation = None
    vts = None
    genes_canonical_transcripts = None
    num_variant_annotation_versions = variant.variantannotation_set.count()
    user_settings = UserSettings.get_for_user(request.user)

    if variant.can_have_annotation:
        try:
            vts = VariantTranscriptSelections(variant, genome_build, annotation_version)
            variant_annotation = vts.variant_annotation
            for w_msg in vts.warning_messages:
                messages.add_message(request, messages.WARNING, w_msg)

            for e_msg in vts.error_messages:
                messages.add_message(request, messages.ERROR, e_msg)

            genes_canonical_transcripts = get_genes_canonical_transcripts(variant, annotation_version)

        except Exception:  # May not have been annotated?
            log_traceback()

    modified_normalised_variants = variant.modifiedimportedvariant_set.filter(old_variant__isnull=False)
    modified_normalised_variants = modified_normalised_variants.values_list("old_variant", flat=True).distinct()

    variant_allele_data = None
    # We don't really handle very rare case with having multiple VariantAlleles for a variant
    if variant_allele := variant.variantallele_set.filter(genome_build=genome_build).first():
        # If we require a ClinGen call (eg have Allele but no ClinGen, and settings say we can get it then don't
        # provide the data so we will do an async call)
        if not variant_allele.needs_clingen_call():
            variant_allele_data = VariantAlleleSerializer.data_with_link_data(variant_allele)
        ClinVarRecordCollection.set_allele_for_variants(variant_allele.allele)

    annotation_description = {}
    if user_settings.tool_tips:
        annotation_description = VariantGridColumn.get_column_descriptions()
        annotation_description["allele"] = "An Allele is genome build independent - ie GRCh37 and GRCh38 variants for" \
                                           " the same change are linked by an allele"
        annotation_description["gnomad"] = "<a href='https://gnomad.broadinstitute.org/'>gnomAD</a> is the world's largest population frequency database"
        annotation_description["mastermind"] = "Mastermind indexes medical literature for variants and genes, see <a href='https://www.genomenon.com/data'>Mastermind</a>"  # Used for no results
        annotation_description["mave_db"] = "<a href='https://www.mavedb.org/'>MaveDB</a> is a public repository for datasets from Multiplexed Assays of Variant Effect (MAVEs), such as those generated by deep mutational scanning or massively parallel reporter assay experiments."
        annotation_description["maxentscan"] = "<a href='http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html'>MaxEntScan</a> scores for human 5 prime splice sites."
        annotation_description["spliceai"] = "Deep Learning splicing predictor - see <a href='https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub'>SpliceAI</a>"

    has_tags = VariantTag.get_for_build(genome_build, variant_qs=variant.equivalent_variants).exists()

    if variant_annotation and variant_annotation.hgvs_g:
        hgvs_g = variant_annotation.hgvs_g
    else:
        hgvs_g = VariantAnnotation.get_hgvs_g(variant)  # Calculate for reference variants

    context = {
        "annotation_description": annotation_description,
        "annotation_version": annotation_version,
        "can_create_classification": Classification.can_create_via_web_form(request.user),
        "genes_canonical_transcripts": genes_canonical_transcripts,
        "genome_build": genome_build,
        "has_tags": has_tags,
        "hgvs_g": hgvs_g,
        "latest_annotation_version": latest_annotation_version,
        "modified_normalised_variants": modified_normalised_variants,
        "num_variant_annotation_versions": num_variant_annotation_versions,
        "tag_form": TagForm(),
        "tool_tips": user_settings.tool_tips,
        "igv_links_enabled": get_settings_form_features().igv_links_enabled,
        "variant": variant,
        "variant_allele": variant_allele_data,
        "variant_annotation": variant_annotation,
        # Names the page and the analysis' variant details tab, which have no room for a transcript
        "variant_short_label": variant_annotation.get_short_label() if variant_annotation else hgvs_g or str(variant),
        "variant_tag_stale_days": user_settings.variant_tag_stale_days,
        "visible_fields": variant_annotation.visible_columns if variant_annotation else frozenset(),
        "vts": vts,
    }
    if extra_context:
        context.update(extra_context)
    return render(request, template, context)


def variant_sample_information(request, variant_id, genome_build_name):
    """ Shell only - the samples grid, locus counts and multi-allelic list are drawn client side
        from the variant_sample_genotypes API, one request per variant """
    variant = get_object_or_404(Variant, pk=variant_id)
    get_genome_build_or_404(genome_build_name)  # Validate, builds we search come from the variants' contigs

    context = {
        "variant": variant,
        "variant_ids": [v.pk for v in variant.all_build_variants],
        "has_samples_in_other_builds":
            Sample.objects.exclude(vcf__genome_build__in=variant.all_genome_builds).exists(),
    }
    return render(request, "variantopedia/variant_sample_information.html", context)


def gene_coverage(request, gene_symbol_id):
    gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol_id)

    twenty_x_coverage_percent_counts = None
    if settings.SEQAUTO_ENABLED:
        twenty_x_coverage_percent_counts = []
        for min_coverage in [1, 50, 99, 100]:
            count = get_20x_gene_coverage(gene_symbol, min_coverage)
            twenty_x_coverage_percent_counts.append((min_coverage, count))

    context = {"gene_symbol": gene_symbol,
               "twenty_x_coverage_percent_counts": twenty_x_coverage_percent_counts}
    return render(request, "variantopedia/gene_coverage.html", context)


def nearby_variants_tab(request, variant_id, annotation_version_id):
    variant = get_object_or_404(Variant, pk=variant_id)
    annotation_version = get_object_or_404(AnnotationVersion, pk=annotation_version_id)
    distance: int = settings.VARIANT_DETAILS_NEARBY_RANGE
    context = {
        'variant': variant,
        "annotation_version": annotation_version,
        "distance": distance,
        "distance_str": str(distance)
    }
    context.update(get_nearby_summaries(request.user, variant, annotation_version,
                                        distance=distance, clinical_significance=True))
    return render(request, "variantopedia/nearby_variants_tab.html", context)


def nearby_variants(request, variant_id, annotation_version_id):
    variant = get_object_or_404(Variant, pk=variant_id)
    annotation_version = get_object_or_404(AnnotationVersion, pk=annotation_version_id)

    variant_annotation_version = annotation_version.variant_annotation_version
    variant_annotation = variant.variantannotation_set.filter(version=variant_annotation_version).first()
    context = {
        "method_summaries": get_method_summaries(variant, annotation_version,
                                                 distance=settings.VARIANT_DETAILS_NEARBY_RANGE),
        "genome_build": annotation_version.genome_build,
        "variant": variant,
        "variant_annotation": variant_annotation
    }
    context.update(get_nearby_qs(variant, annotation_version))
    return render(request, "variantopedia/nearby_variants.html", context)


def _get_grid_name(request, name) -> str:
    name_parts = [name]
    tag_id = request.GET.get("tag")
    if extra_filters := request.GET.get("extra_filters"):
        tag_id = json.loads(extra_filters).get("tag")
    if tag_id:
        name_parts.extend(["tag", tag_id])
    return "_".join(name_parts)


# The DataTable pages at 100 rows, so the variant tags CSV comes from the queryset directly
VARIANT_TAGS_EXPORT_COLUMNS = {
    "variant_string": "Variant",
    "gene_symbol": "Gene",
    "tag__id": "Tag",
    "analysis__id": "Analysis ID",
    "analysis__name": "Analysis",
    "user__username": "Username",
    "created": "Created",
    "variant__id": "Variant ID",
    "id": "Tag ID",
}


def variant_tags_export(request, genome_build_name):
    config = VariantTagsColumns(request)
    qs = config.filter_queryset(config.get_initial_queryset()).order_by("variant_string")
    basename = f"{_get_grid_name(request, 'variant_tags_export')}_{local_date_string()}"

    pseudo_buffer = StashFile()
    writer = csv.writer(pseudo_buffer, dialect='excel', quoting=csv.QUOTE_MINIMAL)

    def iter_rows():
        writer.writerow(["Genome Build"] + list(VARIANT_TAGS_EXPORT_COLUMNS.values()))
        yield pseudo_buffer.value
        for i, row in enumerate(qs.values(*VARIANT_TAGS_EXPORT_COLUMNS), start=1):
            writer.writerow([genome_build_name] + [row[k] for k in VARIANT_TAGS_EXPORT_COLUMNS])
            if i % EXPORT_ROWS_PER_CHUNK == 0:
                yield pseudo_buffer.value
        if remaining := pseudo_buffer.value:
            yield remaining

    response = StreamingHttpResponse(iter_rows(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{slugify(basename)}.csv"'
    return response


