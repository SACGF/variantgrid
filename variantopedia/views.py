import json
import logging
import operator
import re
from collections import defaultdict
from datetime import timedelta
from functools import reduce
from random import randint
from typing import Optional

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User
from django.db.models import Q
from django.forms import model_to_dict
from django.shortcuts import get_object_or_404, render, redirect
from django.urls import reverse
from django.utils import timezone
from django.utils.timesince import timesince
from django.views.decorators.http import require_POST

from analysis.models import VariantTag
from analysis.models.nodes.analysis_node import Analysis
from analysis.serializers import VariantTagSerializer
from annotation.models import AnnotationRun, AnnotationVersion, ClassificationModification, Classification, ClinVar, \
    VariantAnnotationVersion, VariantAnnotation
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from eventlog.models import Event, create_event
from genes.hgvs import HGVSMatcher
from genes.models import CanonicalTranscriptCollection, GeneSymbol
from library.django_utils import require_superuser, highest_pk, get_url_from_view_path, get_field_counts
from library.enums.log_level import LogLevel
from library.git import Git
from library.guardian_utils import admin_bot
from library.log_utils import report_exc_info, log_traceback, NotificationBuilder
from library.utils import count, pretty_label
from pathtests.models import cases_for_user
from patients.models import ExternalPK, Clinician
from seqauto.models import VCFFromSequencingRun, get_20x_gene_coverage
from seqauto.seqauto_stats import get_sample_enrichment_kits_df
from snpdb.clingen_allele import link_allele_to_existing_variants
from snpdb.forms import TagForm
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import Variant, Sample, VCF, get_igv_data, Allele, AlleleMergeLog, VariantAllele, \
    AlleleConversionTool, ImportSource, AlleleOrigin, VariantAlleleSource
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings, UserPageAck
from snpdb.serializers import VariantAlleleSerializer
from snpdb.variant_pk_lookup import VariantPKLookup
from snpdb.variant_sample_information import VariantSampleInformation
from upload.upload_stats import get_vcf_variant_upload_stats
from variantgrid.celery import app
from variantgrid.perm_path import get_visible_url_names
from variantgrid.tasks.server_monitoring_tasks import get_disk_messages
from variantopedia import forms
from variantopedia.interesting_nearby import get_nearby_qs, get_method_summaries
from variantopedia.search import search_data, SearchResults


def variants(request):
    context = {}
    return render(request, "variantopedia/variants.html", context)


def get_dashboard_notices(user: User, days_ago: Optional[int]) -> dict:
    """ returns {} if nothing to show """
    MAX_PAST_DAYS = 30

    start_time = 0
    notice_header = ""

    if days_ago:
        # if days ago was passed in, don't update upa
        days_ago = min(days_ago, MAX_PAST_DAYS)
        start_time = timezone.now() - timedelta(days=days_ago)
        notice_header = f"Since the last {days_ago} day{'s' if days_ago > 1 else ''}"
    else:
        upa, created = UserPageAck.objects.get_or_create(user=user, page_id="server_status")
        if created:
            if user.last_login:
                start_time = user.last_login
                notice_header = f"Since last login ({timesince(start_time)} ago)"
        else:
            start_time = upa.modified
            notice_header = f"Since last visit to this page ({timesince(start_time)} ago)"
            upa.save()  # Update last modified timestamp

        max_days_ago = timezone.now() - timedelta(days=MAX_PAST_DAYS)
        if max_days_ago > start_time:
            start_time = max_days_ago
            notice_header = f"Since the last {MAX_PAST_DAYS} days"

    if user.is_superuser:
        events = Event.objects.filter(date__gte=start_time, severity=LogLevel.ERROR)
    else:
        events = Event.objects.none()

    vcfs = VCF.filter_for_user(user, True).filter(date__gte=start_time)
    analyses = Analysis.filter_for_user(user)
    analyses_created = analyses.filter(created__gte=start_time)
    analyses_modified = analyses.filter(created__lt=start_time, modified__gte=start_time)
    classifications_of_interest = Classification.dashboard_report_classifications_of_interest(since=start_time)
    new_classification_count = Classification.dashboard_report_new_classifications(since=start_time)

    return {
        "notice_header": notice_header,
        "events": events,
        "classifications_of_interest": classifications_of_interest,
        "classifications_created": new_classification_count,
        "vcfs": vcfs,
        "analyses_created": analyses_created,
        "analyses_modified": analyses_modified,
    }


def strip_celery_from_keys(celery_state):
    worker_status = {}
    if celery_state:
        for worker_string, data in celery_state.items():
            m = re.match(".*@(.*?)$", worker_string)
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


def notify_server_status():
    dashboard_notices = get_dashboard_notices(admin_bot(), days_ago=1)
    url = get_url_from_view_path(reverse('server_status')) + '?days=1'

    emoji = ":male-doctor:" if randint(0, 1) else ":female-doctor:"
    nb = NotificationBuilder(message="Health Check", emoji=emoji)
    nb.add_header("Health Check")
    nb.add_markdown(f"URL : <{url}|{url}>")
    nb.add_markdown("*Disk usage*")
    disk_usage = list()
    for _, message in get_disk_messages(info_messages=True):
        disk_usage.append(f":floppy_disk: {message}")
    nb.add_markdown("\n".join(disk_usage), indented=True)
    nb.add_markdown("*In the last 24 hours*")

    keys = set(dashboard_notices.keys())
    keys.discard('events')
    keys.discard('notice_header')
    visible_urls = get_visible_url_names()
    if not visible_urls.get('analyses'):
        for exclude_key in ['vcfs', 'analyses_created', 'analyses_modified']:
            keys.discard(exclude_key)
    sorted_keys = sorted(list(keys))

    lines = list()
    for key in sorted_keys:
        emoji = ":blue_book:"
        if 'analyses' in key:
            emoji = ":orange_book:"
        elif 'vcf' in key:
            emoji = ":green_book:"
        count_display = count(dashboard_notices.get(key))
        if count_display:
            count_display = f"*{count_display}*"

        lines.append(f"{emoji} {count_display} : {pretty_label(key)}")
    nb.add_markdown("\n".join(lines), indented=True)
    nb.send()


@require_superuser
def server_status(request):
    if request.method == "POST":
        notify_server_status()

    celery_workers = {}
    if settings.CELERY_ENABLED:
        # This relies on the services being started with the "-n worker_name" with a separate one for each service
        worker_names = settings.CELERY_WORKER_NAMES.copy()
        if settings.URLS_APP_REGISTER["analysis"]:
            worker_names.extend(settings.CELERY_ANALYSIS_WORKER_NAMES)
        if settings.SEQAUTO_ENABLED:
            worker_names.extend(settings.CELERY_SEQAUTO_WORKER_NAMES)

        i = app.control.inspect()
        pong = strip_celery_from_keys(i.ping())
        active = strip_celery_from_keys(i.active())
        scheduled = strip_celery_from_keys(i.scheduled())

        for worker in worker_names:
            data = pong.get(worker)
            ok = False
            if data:
                if data.get('ok') == 'pong':
                    status = 'ok'
                    ok = True
                else:
                    status = data
            else:
                status = 'ERROR - no workers found'

            celery_workers[worker] = {"status": status,
                                      "ok": ok,
                                      "active": len(active.get(worker, [])),
                                      "scheduled": len(scheduled.get(worker, []))}

    can_access_reference = True
    try:
        for genome_build in GenomeBuild.builds_with_annotation():
            _ = genome_build.reference_fasta  # Throws exception on error
    except (KeyError, FileNotFoundError):
        can_access_reference = False

    try:
        variant_pk_lookup = VariantPKLookup()
        if variant_pk_lookup.redis_check:
            redis_status = "info"
            redis_message = "OK"
        else:
            redis_status = "warning"
            redis_message = "Inserts running - unable to check"
    except ValueError as ve:
        redis_status = "danger"
        redis_message = str(ve)

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
                ar = AnnotationRun.objects.get(annotation_range_lock__version=vav,
                                               annotation_range_lock__min_variant__lte=highest_variant.pk,
                                               annotation_range_lock__max_variant__gte=highest_variant.pk)
                highest_variant_annotated["status"] = "warning"
                highest_variant_annotated["message"] = f"AnnotationRun: {ar}"
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

    days_ago: Optional[int] = None
    try:
        if days_ago_str := request.GET.get('days'):
            days_ago = int(days_ago_str)
    except ValueError:
        pass

    dashboard_notices = get_dashboard_notices(request.user, days_ago)

    context = {"celery_workers": celery_workers,
               "can_access_reference": can_access_reference,
               "redis_status": redis_status,
               "redis_message": redis_message,
               "disk_free": disk_free,
               "highest_variant_annotated": highest_variant_annotated,
               "sample_enrichment_kits_df": sample_enrichment_kits_df,
               "dashboard_notices": dashboard_notices}
    return render(request, "variantopedia/server_status.html", context)


def database_statistics(request):
    max_variant_id = highest_pk(Variant)
    num_vcfs = VCF.objects.all().count()
    num_samples = Sample.objects.all().count()

    vcf_variant_stats_df = get_vcf_variant_upload_stats()
    variant_stats = {}
    for col in ["cumulative_samples", "total_variants", "percent_known"]:
        variant_stats[col] = vcf_variant_stats_df[col].tolist()

    context = {"max_variant_id": max_variant_id,
               "num_vcfs": num_vcfs,
               "num_samples": num_samples,
               "variant_stats": variant_stats}
    return render(request, "variantopedia/database_statistics.html", context)


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

    igv_data = get_igv_data(request.user, genome_build=genome_build)
    extra_context = {"igv_data": igv_data,
                     "in_multiple_genome_builds": in_multiple_genome_builds,
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
    variant_tags_qs = VariantTag.objects.filter(analysis__genome_build=genome_build)
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

    search_results: Optional[SearchResults] = None
    if form.is_valid() and form.cleaned_data['search']:
        search_string = form.cleaned_data['search']
        classify = form.cleaned_data.get('classify')
        try:
            search_results = search_data(request.user, search_string, classify)
            results, _search_types, _search_errors = search_results.non_debug_results, search_results.search_types, search_results.search_errors
            details = f"'{search_string}' calculated {len(results)} results."
            create_event(request.user, 'search', details=details)

            # don't auto load
            if preferred_result := search_results.single_preferred_result():
                return redirect(preferred_result.record)

            # Attempt to give hints on why nothing was found
            if search_results.search_errors:
                for search_type, e, genome_build in search_results.search_errors:
                    messages.add_message(request, messages.ERROR, f"{search_type}: {e} ({genome_build})")

        except Exception as e:
            report_exc_info(extra_data={
                'search_string': search_string,
                'classify': classify
            })
            msg = f"An error occurred during search: {e}"
            messages.add_message(request, messages.ERROR, msg, extra_tags='save-message')

    epk_qs = ExternalPK.objects.values_list("external_type", flat=True)
    external_codes = list(sorted(epk_qs.distinct()))

    context = {"user_settings": user_settings,
               "form": form,
               "search": search_string,
               "search_results": search_results,
               "external_codes": external_codes,
               "variant_vcf_db_prefix": settings.VARIANT_VCF_DB_PREFIX,
               "variant_summary": settings.SEARCH_VARIANT_SHOW_SUMMARY}
    return render(request, "variantopedia/search.html", context)


def wiki(request):
    context = {}
    return render(request, "variantopedia/wiki.html", context)


def view_allele_from_variant(request, variant_id):
    variant = get_object_or_404(Variant, pk=variant_id)
    allele = variant.allele
    if allele and settings.PREFER_ALLELE_LINKS:
        return redirect(reverse('view_allele', kwargs={"pk": allele.id}))
    return redirect(reverse('view_variant', kwargs={"variant_id": variant_id}))


def view_allele(request, pk):
    allele: Allele = get_object_or_404(Allele, pk=pk)
    link_allele_to_existing_variants(allele, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)

    latest_classifications = ClassificationModification.latest_for_user(
        user=request.user,
        allele=allele,
        published=True
    )

    allele_merge_log_qs = AlleleMergeLog.objects.filter(Q(old_allele=allele) | Q(new_allele=allele)).order_by("pk")
    context = {"allele": allele,
               "edit_clinical_groupings": request.GET.get('edit_clinical_groupings') == 'True',
               "allele_merge_log_qs": allele_merge_log_qs,
               "clingen_url": settings.CLINGEN_ALLELE_REGISTRY_DOMAIN,
               "classifications": latest_classifications,
               "annotated_builds": GenomeBuild.builds_with_annotation()}
    return render(request, "variantopedia/view_allele.html", context)


@require_POST
def create_variant_for_allele(request, allele_id, genome_build_name):
    """ Shortcut to create manual variant, but as a POST """
    allele = get_object_or_404(Allele, pk=allele_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    non_liftover_origin = [AlleleOrigin.IMPORTED_TO_DATABASE, AlleleOrigin.IMPORTED_NORMALIZED]
    if variant_allele := allele.variantallele_set.filter(origin__in=non_liftover_origin).first():
        allele_source = VariantAlleleSource.objects.create(variant_allele=variant_allele)
        create_liftover_pipelines(admin_bot(), allele_source, ImportSource.WEB, variant_allele.genome_build, [genome_build])
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
    g_hgvs = HGVSMatcher(genome_build).variant_to_g_hgvs(variant)
    latest_annotation_version = AnnotationVersion.latest(genome_build)
    variant_annotation = None
    vts = None
    clinvar = None
    genes_canonical_transcripts = None
    num_clinvar_citations = 0
    clinvar_citations = None
    num_variant_annotation_versions = variant.variantannotation_set.count()

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

            clinvar_qs = ClinVar.objects.filter(variant=variant, version=annotation_version.clinvar_version)
            try:
                clinvar = clinvar_qs.get()
            except ClinVar.MultipleObjectsReturned:
                # Report this - but carry on for the user
                details = f"Multiple ClinVar entries found for variant {variant.pk}"
                create_event(request.user, "duplicate_annotation", details, LogLevel.WARNING)
                clinvar = clinvar_qs.first()
            except ClinVar.DoesNotExist:
                pass

            if clinvar:
                clinvar_citations = [{'idx': c.citation_id, 'db': c.get_citation_source_display()} for c in clinvar.get_citations()]
                num_clinvar_citations = len(clinvar_citations)
        except:  # May not have been annotated?
            log_traceback()

    modified_normalised_variants = variant.modifiedimportedvariant_set.all().filter(old_variant__isnull=False)
    modified_normalised_variants = modified_normalised_variants.values_list("old_variant", flat=True).distinct()

    try:
        variant_allele = variant.variantallele_set.get(genome_build=genome_build)
        # If we require a ClinGen call (eg have Allele but no ClinGen, and settings say we can get it then don't
        # provide the data so we will do an async call)
        if not variant_allele.needs_clinvar_call():
            logging.info("Skipping as we need clinvar call")
            variant_allele_data = VariantAlleleSerializer.data_with_link_data(variant_allele)
        else:
            variant_allele_data = None
    except VariantAllele.DoesNotExist:
        variant_allele_data = None

    variant_tags = []
    for vt in VariantTag.get_for_build(genome_build, variant_qs=variant.equivalent_variants):
        variant_tags.append(VariantTagSerializer(vt, context={"request": request}).data)

    context = {
        "ANNOTATION_PUBMED_SEARCH_TERMS_ENABLED": settings.ANNOTATION_PUBMED_SEARCH_TERMS_ENABLED,
        "annotation_version": annotation_version,
        "can_create_classification": Classification.can_create_via_web_form(request.user),
        "classifications": latest_classifications,
        "clinvar": clinvar,
        "genes_canonical_transcripts": genes_canonical_transcripts,
        "genome_build": genome_build,
        "g_hgvs": g_hgvs,
        "latest_annotation_version": latest_annotation_version,
        "modified_normalised_variants": modified_normalised_variants,
        "num_variant_annotation_versions": num_variant_annotation_versions,
        "num_clinvar_citations": num_clinvar_citations,
        "clinvar_citations": clinvar_citations,
        "show_annotation": settings.VARIANT_DETAILS_SHOW_ANNOTATION,
        "show_samples": settings.VARIANT_DETAILS_SHOW_SAMPLES,
        "tag_form": TagForm(),
        "variant": variant,
        "variant_allele": variant_allele_data,
        "variant_annotation": variant_annotation,
        "variant_tags": json.dumps(variant_tags),
        "vts": vts,
    }
    if extra_context:
        context.update(extra_context)
    return render(request, template, context)


def variant_sample_information(request, variant_id, genome_build_name):
    variant = get_object_or_404(Variant, pk=variant_id)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    vsi = VariantSampleInformation(request.user, variant, genome_build)

    context = {"variant": variant,
               "vsi": vsi,
               "visible_rows": vsi.visible_rows}
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


def nearby_variants(request, variant_id, annotation_version_id):
    variant = get_object_or_404(Variant, pk=variant_id)
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)

    variant_annotation_version = annotation_version.variant_annotation_version
    variant_annotation = variant.variantannotation_set.filter(version=variant_annotation_version).first()
    context = {
        "method_summaries": get_method_summaries(variant, distance=settings.VARIANT_DETAILS_NEARBY_RANGE),
        "genome_build": annotation_version.genome_build,
        "variant": variant,
        "variant_annotation": variant_annotation
    }
    context.update(get_nearby_qs(variant, annotation_version))
    return render(request, "variantopedia/nearby_variants.html", context)
