import json
import operator
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta
from functools import reduce
from typing import Any

from django.conf import settings
from django.contrib import messages
from django.db import connection
from django.forms import model_to_dict
from django.shortcuts import get_object_or_404, render, redirect
from django.urls import reverse
from django.utils import timezone
from django.utils.timezone import localtime
from django.views.decorators.http import require_POST

from analysis.models import VariantTag
from annotation.models import AnnotationRun, AnnotationVersion, ClassificationModification, Classification, \
    VariantAnnotationVersion, VariantAnnotation, AnnotationStatus
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from classification.models.classification_import_run import ClassificationImportRun
from classification.variant_card import AlleleCard
from classification.views.exports import ClassificationExportFormatterCSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_csv import FormatDetailsCSV
from eventlog.models import create_event
from genes.hgvs import HGVSMatcher
from genes.models import CanonicalTranscriptCollection, GeneSymbol
from library.django_utils import require_superuser, highest_pk, get_field_counts
from library.django_utils.jqgrid_view import JQGridView
from library.git import Git
from library.guardian_utils import admin_bot
from library.health_check import HealthCheckRequest, health_check_overall_stats_signal
from library.log_utils import log_traceback, report_message, slack_bot_username, AdminNotificationBuilder
from library.utils import flatten_nested_lists
from pathtests.models import cases_for_user
from patients.models import Clinician
from seqauto.models import VCFFromSequencingRun, get_20x_gene_coverage
from seqauto.seqauto_stats import get_sample_enrichment_kits_df
from snpdb.clingen_allele import link_allele_to_existing_variants
from snpdb.forms import TagForm, get_settings_form_features
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import Variant, Sample, VCF, get_igv_data, Allele, AlleleConversionTool, ImportSource, AlleleOrigin, \
    VariantGridColumn, Tag
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.search import search_data
from snpdb.serializers import VariantAlleleSerializer
from snpdb.variant_sample_information import VariantSampleInformation
from upload.models import ModifiedImportedVariant
from upload.upload_stats import get_vcf_variant_upload_stats
from variantgrid.celery import app
from variantgrid.tasks.server_monitoring_tasks import get_disk_messages
from variantopedia import forms
from variantopedia.grids import VariantTagsGrid, TaggedVariantGrid
from variantopedia.interesting_nearby import get_nearby_qs, get_method_summaries, get_nearby_summaries
from variantopedia.server_status import get_dashboard_notices
from variantopedia.tasks.server_status_tasks import notify_server_status_now


def variants(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    context = {"genome_build": genome_build}
    return render(request, "variantopedia/variants.html", context)


def strip_celery_from_keys(celery_state):
    worker_status = {}
    if celery_state:
        for worker_string, data in celery_state.items():
            m = re.match(r".*@(.*?)$", worker_string)
            if m:
                worker = m.group(1)
                worker_status[worker] = data
    return worker_status


def dashboard(request):
    if Clinician.user_is_clinician(request.user):
        return redirect('clinician_login')

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


@require_superuser
def server_status(request):
    if request.method == "POST":
        action = request.POST.get('action')
        if action == 'Test Slack':
            nb = AdminNotificationBuilder(message="Slack Check")
            nb.add_markdown("This is a Slack Test :ladybug:")
            nb.send()
            messages.add_message(request, level=messages.INFO, message="Slack should have been sent a test message.")
        if action == 'Health Check':
            notify_server_status_now()
            messages.add_message(request, level=messages.INFO, message="Slack should have been sent the health check.")
        elif action == 'Test Rollbar':
            report_message("Testing Rollbar", level='error')
            messages.add_message(request, level=messages.INFO, message="Rollbar should have been sent an error.")
        elif action == 'Test Message Branding':
            messages.success(request, "Success message")
            messages.info(request, "Info message")
            messages.warning(request, "Warning message")
            messages.error(request, "Error message")

        elif action == 'kill-pid':
            pid = int(request.POST.get('pid'))
            with connection.cursor() as cursor:
                cursor.execute("SELECT pg_terminate_backend(%s)", [pid])
                terminated = cursor.fetchone()[0]
                messages.add_message(request, level=messages.INFO, message=f"Query {pid} Terminated = {terminated}")
        else:
            print(f"Unrecognised action {action}")

        # return redirect(reverse('server_status'))

        # TODO should redirect to read-only version of the page

    celery_workers = {}
    if settings.CELERY_ENABLED:
        # This relies on the services being started with the "-n worker_name" with a separate one for each service
        worker_names = settings.CELERY_WORKER_NAMES.copy()
        if settings.URLS_APP_REGISTER["analysis"]:
            worker_names.extend(settings.CELERY_ANALYSIS_WORKER_NAMES)
        if settings.SEQAUTO_ENABLED:
            worker_names.extend(settings.CELERY_SEQAUTO_WORKER_NAMES)

        i = app.control.inspect()
        ping = strip_celery_from_keys(i.ping())
        stats = strip_celery_from_keys(i.stats())
        active = strip_celery_from_keys(i.active())
        scheduled = strip_celery_from_keys(i.scheduled())

        for worker in worker_names:
            num_workers = "?"
            status = 'ERROR - no workers found'
            ok = False

            # Sometimes stats fails - just use ping
            pong = ping.get(worker, {})
            if pong.get("ok") == "pong":
                status = "OK"
                ok = True

            data = stats.get(worker)
            if data:
                processes = data.get("pool", {}).get("processes")
                if processes:
                    num_workers = len(processes)
                    status = "OK"
                    ok = True

            celery_workers[worker] = {
                "status": status,
                "ok": ok,
                "active": len(active.get(worker, [])),
                "scheduled": len(scheduled.get(worker, [])),
                "num_workers": num_workers,
            }

    can_access_reference = True
    try:
        for genome_build in GenomeBuild.builds_with_annotation():
            _ = genome_build.reference_fasta  # Throws exception on error
    except (KeyError, FileNotFoundError):
        can_access_reference = False

    # Variant Annotation - incredibly quick check
    highest_variant_annotated = {}
    try:
        q = reduce(operator.and_, VariantAnnotation.VARIANT_ANNOTATION_Q)
        highest_variant = Variant.objects.filter(q).order_by("pk").last()
        genome_build = next(iter(highest_variant.genome_builds))  # Just pick one if spans multiple
        vav = VariantAnnotationVersion.latest(genome_build)
        annotated = highest_variant.variantannotation_set.filter(version=vav).exists()
        if annotated:
            highest_variant_annotated["status"] = "info"
            highest_variant_annotated["message"] = "OK"
        else:
            try:
                ar_qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav,
                                                     annotation_range_lock__min_variant__lte=highest_variant.pk,
                                                     annotation_range_lock__max_variant__gte=highest_variant.pk)
                ar_qs = ar_qs.exclude(status=AnnotationStatus.ERROR)
                annotation_run_message = f"AnnotationRuns: {', '.join([str(ar) for ar in ar_qs])}"

                highest_variant_annotated["status"] = "warning"
                highest_variant_annotated["message"] = annotation_run_message
            except AnnotationRun.DoesNotExist:
                highest_variant_annotated["status"] = "danger"
                highest_variant_annotated["message"] = "Not annotated, no AnnotationRun!"
    except Exception as e:
        highest_variant_annotated["status"] = "danger"
        highest_variant_annotated["message"] = str(e)

    sample_enrichment_kits_df = None
    if settings.SEQAUTO_ENABLED:
        sample_enrichment_kits_df = get_sample_enrichment_kits_df()
    disk_messages = get_disk_messages(info_messages=True)
    disk_free = {"status": "info", "messages": []}
    for status, message in disk_messages:
        if status == "warning":
            disk_free["status"] = "warning"
        disk_free["messages"].append(message)

    context = {
        "celery_workers": celery_workers,
        "queries": long_running_sql(0),
        "can_access_reference": can_access_reference,
        "highest_variant_annotated": highest_variant_annotated,
        "sample_enrichment_kits_df": sample_enrichment_kits_df
    }
    return render(request, "variantopedia/server_status.html", context)


@require_superuser
def server_status_activity(request, days_ago: int):
    dashboard_notices = get_dashboard_notices(request.user, days_ago)
    return render(request, "variantopedia/server_status_activity_detail.html", {"dashboard_notices": dashboard_notices})


@require_superuser
def server_status_settings(request):
    slack_emoji = (settings.SLACK or {}).get('emoji') or ':dna:'
    slack_username = f"{slack_emoji} {slack_bot_username()}"

    hgvs_matcher_id = settings.HGVS_DEFAULT_METHOD
    hgvs_matcher = HGVSMatcher(GenomeBuild.grch38(), hgvs_converter_type=hgvs_matcher_id)

    return render(request, "variantopedia/server_status_settings_detail.html", {
        "settings": settings,
        "slack_bot_username": slack_username,
        "ongoing_imports": ClassificationImportRun.ongoing_imports(),
        "hgvs_matcher": hgvs_matcher
    })


@require_superuser
def health_check_details(request):
    now = localtime()
    since = now - timedelta(days=1)
    health_request = HealthCheckRequest(since=since, now=now)

    results = []
    for _, result in health_check_overall_stats_signal.send_robust(sender=None, health_request=health_request):
        if not isinstance(result, Exception):
            results.append(result)

    checks = flatten_nested_lists(results)
    checks = sorted(checks, key=lambda hc: (hc.sort_order(), hc.name if hasattr(hc, "name") else "z"))
    overall_lines = []
    for check in checks:
        line_content = check.as_html()
        overall_lines.append(line_content)

    context = {
        'overall_lines': overall_lines,
    }
    return render(request, "variantopedia/health_check_details.html", context)


@dataclass
class RunningQuery:
    pid: int
    duration: Any  # is actually a Duration
    query: str
    state: str


def long_running_sql(min_age_in_seconds: int = 30):
    with connection.cursor() as cursor:
        db_name = connection.settings_dict['NAME']
        cursor.execute(
            """
            SELECT
              pid,
              now() - pg_stat_activity.query_start AS duration,
              query,
              state
            FROM pg_stat_activity
            WHERE (now() - pg_stat_activity.query_start) > interval %s
            AND datname = %s
            ORDER BY now() - pg_stat_activity.query_start desc;
            """,
            [f"{min_age_in_seconds} seconds", db_name]
        )

        def to_obj(row) -> RunningQuery:
            return RunningQuery(
                pid=row[0],
                duration=row[1],
                query=row[2],
                state=row[3]
            )

        return [to_obj(result) for result in cursor.fetchall()]


@require_superuser
def database_statistics(request):
    max_variant_id = highest_pk(Variant)
    num_vcfs = VCF.objects.count()
    num_samples = Sample.objects.count()

    vcf_variant_stats_df = get_vcf_variant_upload_stats()
    variant_stats = {}
    for col in ["cumulative_samples", "cumulative_variants", "percent_known"]:
        variant_stats[col] = vcf_variant_stats_df[col].tolist()

    context = {"max_variant_id": max_variant_id,
               "num_vcfs": num_vcfs,
               "num_samples": num_samples,
               "variant_stats": variant_stats}
    return render(request, "variantopedia/database_statistics_detail.html", context)


def variant_tag_detail(request, variant_id, tag):
    """ Loaded via tags grid on variant page """

    variant = Variant.objects.get(pk=variant_id)
    tag = get_object_or_404(Tag, pk=tag)
    context = {
        "variant": variant,
        "tag": tag,
    }
    return render(request, "variantopedia/variant_tag_detail.html", context)


def view_variant(request, variant_id, genome_build_name=None):
    """ This is to open it with the normal menu around it (ie via search etc) """
    template = 'variantopedia/view_variant.html'

    variant = get_object_or_404(Variant, pk=variant_id)
    in_multiple_genome_builds = len(variant.genome_builds) > 1
    genome_build = None
    if genome_build_name:
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        if genome_build not in variant.genome_builds:
            raise ValueError(f"Variant {variant} not in '{genome_build}'")
    else:
        if in_multiple_genome_builds:
            user_settings = UserSettings.get_for_user(request.user)
            if user_settings.default_genome_build in variant.genome_builds:
                genome_build = user_settings.default_genome_build

    if genome_build is None:
        genome_build = next(iter(variant.genome_builds))

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
    context = {"genome_build": genome_build,
               "tag_counts": tag_counts}
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


def view_allele(request, allele_id: int):
    allele: Allele = get_object_or_404(Allele, pk=allele_id)
    link_allele_to_existing_variants(allele, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)

    latest_classifications = ClassificationModification.latest_for_user(
        user=request.user,
        allele=allele,
        published=True
    ).select_related(
        'classification__allele_info__allele',
        'classification__allele_info__imported_genome_build_patch_version',
        'classification__allele_info__latest_validation'
    )

    context = {
        "allele_card": AlleleCard(user=request.user, allele=allele),
        "allele": allele,
        "edit_clinical_groupings": request.GET.get('edit_clinical_groupings') == 'True',
        "classifications": latest_classifications
    }
    return render(request, "variantopedia/view_allele.html", context)


def export_classifications_allele(request, allele_id: int):
    """
    CSV export of what is currently filtered into the classification grid
    """
    allele = Allele.objects.get(pk=allele_id)
    return ClassificationExportFormatterCSV(
        ClassificationFilter(
            user=request.user,
            genome_build=GenomeBuildManager.get_current_genome_build(),
            allele=allele_id,
            file_prefix=f"classifications_allele_{allele:CA}"
        ),
        FormatDetailsCSV()
    ).serve()


@require_POST
def create_variant_for_allele(request, allele_id, genome_build_name):
    """ Shortcut to create manual variant, but as a POST """
    allele = get_object_or_404(Allele, pk=allele_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    non_liftover_origin = [AlleleOrigin.IMPORTED_TO_DATABASE, AlleleOrigin.IMPORTED_NORMALIZED]
    if variant_allele := allele.variantallele_set.filter(origin__in=non_liftover_origin).first():
        create_liftover_pipelines(admin_bot(), [allele], ImportSource.WEB, variant_allele.genome_build, [genome_build])
    return redirect(allele)


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
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
    genome_build = annotation_version.genome_build
    latest_annotation_version = AnnotationVersion.latest(genome_build)
    variant_annotation = None
    vts = None
    genes_canonical_transcripts = None
    num_variant_annotation_versions = variant.variantannotation_set.count()
    user_settings = UserSettings.get_for_user(request.user)

    latest_classifications = ClassificationModification.latest_for_user(
        user=request.user,
        variant=variant.equivalent_variants,
        published=True
    )

    if variant.can_have_annotation:
        try:
            vts = VariantTranscriptSelections(variant, genome_build, annotation_version)
            variant_annotation = vts.variant_annotation
            for w_msg in vts.warning_messages:
                messages.add_message(request, messages.WARNING, w_msg)

            for e_msg in vts.error_messages:
                messages.add_message(request, messages.ERROR, e_msg)

            genes_canonical_transcripts = get_genes_canonical_transcripts(variant, annotation_version)

        except:  # May not have been annotated?
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

    annotation_description = {}
    if user_settings.tool_tips:
        annotation_description = VariantGridColumn.get_column_descriptions()
        annotation_description["allele"] = "An Allele is genome build independent - ie GRCh37 and GRCh38 variants for" \
                                           " the same change are linked by an allele"
        annotation_description["maxentscan"] = "<a href='http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html'>MaxEntScan</a> scores for human 5 prime splice sites."
        annotation_description["spliceai"] = "Deep Learning splicing predictor - see <a href='https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub'>SpliceAI</a>"

    has_tags = VariantTag.get_for_build(genome_build, variant_qs=variant.equivalent_variants).exists()

    context = {
        "annotation_description": annotation_description,
        "annotation_version": annotation_version,
        "can_create_classification": Classification.can_create_via_web_form(request.user),
        "classifications": latest_classifications,
        "genes_canonical_transcripts": genes_canonical_transcripts,
        "genome_build": genome_build,
        "has_tags": has_tags,
        "latest_annotation_version": latest_annotation_version,
        "modified_normalised_variants": modified_normalised_variants,
        "num_variant_annotation_versions": num_variant_annotation_versions,
        "tag_form": TagForm(),
        "tool_tips": user_settings.tool_tips,
        "igv_links_enabled": get_settings_form_features().igv_links_enabled,
        "variant": variant,
        "variant_allele": variant_allele_data,
        "variant_annotation": variant_annotation,
        "vts": vts,
    }
    if extra_context:
        context.update(extra_context)
    return render(request, template, context)


def variant_sample_information(request, variant_id, genome_build_name):
    variant = get_object_or_404(Variant, pk=variant_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    vsi = VariantSampleInformation(request.user, variant, genome_build)
    other_loci_variants_by_multiallelic = ModifiedImportedVariant.get_other_loci_variants_by_multiallelic(variant)

    context = {
        "variant": variant,
        "vsi": vsi,
        "visible_rows": vsi.visible_rows,
        "other_loci_variants_by_multiallelic": other_loci_variants_by_multiallelic,
        "has_samples_in_other_builds": Sample.objects.exclude(vcf__genome_build=genome_build).exists(),
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
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)

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
    name_parts = [
        name,
        timezone.now().strftime('%Y-%m-%d'),
    ]
    if extra_filters := request.GET.get("extra_filters"):
        extra_filters = json.loads(extra_filters)
        if tag_id := extra_filters.get("tag"):
            name_parts.extend(["tag", tag_id])
    return "_".join(name_parts)


def variant_tags_export(request, genome_build_name):
    class SortVariantTagsGrid(VariantTagsGrid):
        def _get_sidx_and_sord(self, request) -> tuple:
            return "variant_string", "asc"

    basename = _get_grid_name(request, "variant_tags_export")
    return JQGridView.export_grid_as_csv(request, grid_klass=SortVariantTagsGrid,
                                         basename=basename, genome_build_name=genome_build_name)


def tagged_variant_export(request, genome_build_name):
    class SortTaggedVariantGrid(TaggedVariantGrid):
        def _get_sidx_and_sord(self, request) -> tuple:
            return "locus__position", "asc"

    basename = _get_grid_name(request, "tagged_variant_export")
    return JQGridView.export_grid_as_csv(request, grid_klass=SortTaggedVariantGrid,
                                         basename=basename, genome_build_name=genome_build_name)
