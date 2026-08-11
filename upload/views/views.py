import json
import logging
import operator
import os
from datetime import date, timedelta
from functools import cached_property, reduce

from django.conf import settings
from django.contrib import messages
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.db.models import Q
from django.http.response import HttpResponse, HttpResponseRedirect, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse
from django.utils import timezone
from django.utils.timesince import timesince
from django.views.decorators.cache import never_cache
from django.views.decorators.http import require_http_methods, require_POST
from django_downloadview import PathDownloadView

from analysis.models import AnalysisTemplate
from annotation.models import AnnotationRun
from annotation.views import get_build_contigs
from eventlog.models import create_event
from library.django_utils.file_uploads import filepond_process_response, filepond_upload_receive
from library.enums.log_level import LogLevel
from library.log_utils import log_traceback
from library.utils.django_utils import render_ajax_view
from snpdb.models import VCF
from snpdb.models.models_enums import ImportStatus
from upload import forms, upload_processing, upload_stats
from upload.models import (
    FileUpload,
    ImportSource,
    ModifiedImportedVariant,
    ProcessingStatus,
    SimpleVCFImportInfo,
    TimeFilterMethod,
    UploadedFileTypes,
    UploadPipeline,
    UploadSettings,
    UploadStep,
    VCFImportInfo,
    VCFSkippedContigs,
)
from upload.uploaded_file_type import (
    get_import_tasks_by_extension,
    get_upload_data_for_uploaded_file,
    get_uploaded_file_type,
    get_url_and_data_for_uploaded_file_data,
    retry_upload_pipeline,
)
from upload.upload_metadata import (
    UploadMetadataError,
    get_metadata_keys_for_file_type,
    validate_upload_metadata,
)

UPLOADED_FILE_CONTEXT = {UploadedFileTypes.VCF: "uploaded_vcf",
                         UploadedFileTypes.GENE_LIST: "uploaded_gene_list"}


def get_status_icon(status):
    ICONS = {
        ProcessingStatus.CREATED: {'icon': 'fa-clock', 'title': 'Queued'},
        ProcessingStatus.PROCESSING: {'icon': 'fa-spinner fa-spin', 'title': 'Processing'},
        ProcessingStatus.ERROR: {'icon': 'fa-times-circle', 'css': 'text-danger', 'title': 'Error'},
        ProcessingStatus.SUCCESS: {'icon': 'fa-check-circle', 'css': 'text-success', 'title': 'Success'},
        ProcessingStatus.TERMINATED_EARLY: {'icon': 'fa-exclamation-triangle', 'css': 'text-warning', 'title': 'Terminated early'},
    }
    return ICONS.get(status, {})

def _get_basic_uploaded_file_context(file_upload) -> dict:
    data_url, upload_data = get_url_and_data_for_uploaded_file_data(file_upload)
    file_type = None
    if file_upload.file_type:
        file_type = UploadedFileTypes(file_upload.file_type).label

    data = {
        'file_type': file_type,
        'file_type_code': file_upload.file_type,
        'data_url': data_url,
    }
    if upload_data:
        data["upload_data"] = upload_data.get_upload_context()
    return data


def uploadedfile_dict(file_upload) -> dict:
    try:
        size = file_upload.file_field.size
    except Exception:
        size = None

    time_since = timesince(file_upload.created)

    file_upload_id = file_upload.pk
    data = {
        'file_upload_id': file_upload_id,
        'uploaded_file_id': file_upload_id,  # deprecated alias
        'name': file_upload.name,
        'size': size,
        'user': file_upload.user.get_full_name(),
        'time_since': f"{time_since} ago",
        'deleteUrl': reverse('upload_file_delete', kwargs={'pk': file_upload.pk}),
        'deleteType': 'POST',
    }
    data.update(_get_basic_uploaded_file_context(file_upload))

    if not file_upload.file_type:
        data['error'] = f'Could not determine how to read file: "{file_upload.name}"'

    try:
        upload_pipeline = UploadPipeline.objects.get(file_upload=file_upload)
        data["upload_pipeline_id"] = upload_pipeline.pk
        try:
            if upload_pipeline.genome_build:
                data["genome_build"] = str(upload_pipeline.genome_build)
                if upload_pipeline.status == ProcessingStatus.PROCESSING:
                    try:
                        uploaded_vcf = file_upload.uploadedvcf
                        data['remaining_annotation_runs'] = get_remaining_annotation_runs(uploaded_vcf, upload_pipeline.genome_build)
                    except ObjectDoesNotExist:
                        pass
        except Exception:
            pass  # Genome build is optional

        status = upload_pipeline.status
        url = reverse('view_upload_pipeline', kwargs={'upload_pipeline_id': upload_pipeline.pk})
    except Exception:
        status = ProcessingStatus.ERROR
        url = reverse('view_uploaded_file', kwargs={'file_upload_id': file_upload.pk})

    data['processing_status'] = status
    data['status_icon'] = get_status_icon(status)
    data["url"] = url
    return data


def get_remaining_annotation_runs(uploaded_vcf, genome_build) -> int:
    max_variant_id = uploaded_vcf.max_variant_id
    if max_variant_id is None:
        # VCF not fully imported yet, so highest known variant is unknown - no remaining runs to report
        return 0
    ar_qs = AnnotationRun.get_active_runs(genome_build)
    return ar_qs.filter(annotation_range_lock__max_variant_id__lte=max_variant_id).count()


def handle_file_upload(user, django_uploaded_file, path=None, metadata=None) -> FileUpload:
    original_filename = django_uploaded_file._name
    kwargs = {
        "name": original_filename,
        "file_field": django_uploaded_file,
        "import_source": ImportSource.WEB_UPLOAD,
        "user": user,
        "path": path,
    }
    file_upload = FileUpload.objects.create(**kwargs)
    # Save 1st to actually create file (need to open handling unicode)
    file_upload.file_type = get_uploaded_file_type(file_upload, original_filename)

    # Validate while the client is still connected - the file type is known by here, so a bad key or
    # an unresolvable build is theirs to fix now rather than a failed import several stages later
    try:
        file_upload.metadata = validate_upload_metadata(metadata,
                                                        get_metadata_keys_for_file_type(file_upload.file_type))
    except UploadMetadataError:
        file_upload.delete()
        raise
    file_upload.save()

    # File is on disk now - store hash so uploads can be de-duped / polled by content (API + web)
    file_upload.store_sha256_hash()

    if file_upload.file_type:
        upload_processing.process_uploaded_file(file_upload)
    return file_upload


def _cohort_export_templates_configured() -> bool:
    """ Whether this deployment can produce annotated cohort downloads (VCF/CSV export). """
    try:
        AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
        return True
    except ValueError:
        return False


def get_upload_status_dict(file_upload) -> dict:
    """ Token-API status payload for a FileUpload - import + annotation progress and, for VCFs,
        the resulting vcf/samples plus whether annotated downloads are ready. """
    file_type = None
    if file_upload.file_type:
        file_type = UploadedFileTypes(file_upload.file_type).label

    file_upload_id = file_upload.pk
    data = {
        "file_upload_id": file_upload_id,
        "uploaded_file_id": file_upload_id,  # deprecated alias
        "sha256_hash": file_upload.sha256_hash,
        "file_type": file_type,
        "pipeline_status": None,
        "progress_percent": None,
        "import_status": None,
        "remaining_annotation_runs": None,
        "annotation_complete": False,
        "vcf_id": None,
        "samples": [],
        "error": None,
        "downloads_available": False,
    }

    upload_pipeline = UploadPipeline.objects.filter(file_upload=file_upload).first()
    if upload_pipeline is None:
        data["error"] = f'Could not determine how to read file: "{file_upload.name}"'
        return data

    data["pipeline_status"] = ProcessingStatus(upload_pipeline.status).label
    data["progress_percent"] = upload_pipeline.progress_percent
    if upload_pipeline.status == ProcessingStatus.ERROR:
        data["error"] = upload_pipeline.progress_status

    try:
        uploaded_vcf = file_upload.uploadedvcf
    except ObjectDoesNotExist:
        uploaded_vcf = None

    remaining_annotation_runs = None
    if uploaded_vcf:
        # uploadedvcf can exist before its vcf is created (import still starting) - guard against
        # the resulting race, as UploadPipeline.genome_build dereferences uploadedvcf.vcf
        if vcf := uploaded_vcf.vcf:
            data["vcf_id"] = vcf.pk
            data["import_status"] = ImportStatus(vcf.import_status).label
            data["samples"] = [{"sample_id": s.pk, "name": s.name}
                               for s in vcf.sample_set.all()]
            if genome_build := upload_pipeline.genome_build:
                remaining_annotation_runs = get_remaining_annotation_runs(uploaded_vcf, genome_build)
            data["remaining_annotation_runs"] = remaining_annotation_runs

    annotation_complete = (upload_pipeline.status == ProcessingStatus.SUCCESS
                           and (remaining_annotation_runs or 0) == 0)
    data["annotation_complete"] = annotation_complete
    data["downloads_available"] = annotation_complete and _cohort_export_templates_configured()
    return data

@require_POST
def upload_file(request):
    """ FilePond ``process`` endpoint: receive one file and start the import pipeline. """
    if not settings.UPLOAD_ENABLED:
        raise PermissionDenied("Uploads are currently disabled (settings.UPLOAD_ENABLED=False)")

    try:
        try:
            django_uploaded_file = filepond_upload_receive(request)
        except ValueError as e:
            create_event(request.user, str(e), severity=LogLevel.ERROR)
            raise

        file_upload = handle_file_upload(request.user, django_uploaded_file)
    except Exception as e:
        logging.error(e)
        log_traceback()
        return HttpResponse("Upload failed. Please try again or contact support.", status=500)

    return filepond_process_response(file_upload.pk)


@require_http_methods(["DELETE", "POST"])
def upload_file_delete(request, pk):
    """ FilePond ``revert`` endpoint (also used by table-row delete on already-processed files). """
    try:
        instance = FileUpload.objects.get(pk=pk)
    except FileUpload.DoesNotExist:
        return HttpResponse(status=404)

    if not (request.user.is_superuser or request.user == instance.user):
        raise PermissionDenied(f"You don't own uploaded file {pk}")

    instance.delete()
    return HttpResponse(status=200)


def get_file_dicts_list(upload_settings):
    file_types = upload_settings.uploadsettingsfiletype_set.values_list("file_type", flat=True)
    filters = [Q(file_type__in=file_types)]
    if not upload_settings.user.is_superuser:
        filters.append(Q(user=upload_settings.user))

    if upload_settings.time_filter_method == TimeFilterMethod.DAYS:
        start_date = date.today() - timedelta(days=upload_settings.time_filter_value)
        filters.append(Q(created__gte=start_date))

    q = reduce(operator.and_, filters)
    qs = FileUpload.objects.filter(q).order_by("-created")
    if upload_settings.time_filter_method == TimeFilterMethod.RECORDS:
        qs = qs[:upload_settings.time_filter_value]

    file_dicts = []
    for file_upload in qs:
        file_dicts.append(uploadedfile_dict(file_upload))
    file_dicts = list(reversed(file_dicts))  # render newest-first
    return file_dicts


@never_cache
def upload_poll(request):
    upload_settings, _ = UploadSettings.objects.get_or_create(user=request.user)
    return JsonResponse(get_file_dicts_list(upload_settings), safe=False)


def upload(request):
    upload_settings, _ = UploadSettings.objects.get_or_create(user=request.user)

    form = forms.UploadSettingsForm(request.POST or None, instance=upload_settings, user=request.user)
    if request.method == "POST":
        if form.is_valid():
            upload_settings = form.save()

    file_dicts_list = get_file_dicts_list(upload_settings)
    extensions = get_import_tasks_by_extension().keys()
    accept_file_types = fr"/(\.|\/)({'|'.join(extensions)})$/i"
    context = {'existing_files': file_dicts_list,
               'form': form,
               "upload_enabled": settings.UPLOAD_ENABLED,
               "accept_file_types": accept_file_types}
    return render(request, 'upload/upload.html', context)


def view_uploaded_file(request, file_upload_id):
    file_upload = get_object_or_404(FileUpload, pk=file_upload_id)
    file_upload.check_can_view(request.user)
    context = {'file_upload': file_upload}
    return render(request, 'upload/view_uploaded_file.html', context)


def view_upload_stats(request):
    DAYS = 30

    end_date = timezone.now() - timedelta(days=DAYS)
    ups_qs = UploadStep.objects.filter(end_date__gte=end_date, name__isnull=False)
    us = list(map(json.dumps, upload_stats.get_upload_stats(ups_qs, max_step_names=10)))
    num_upload_pipelines, total_times, time_per_kilo_variant = us

    context = {"days": DAYS,
               "num_upload_pipelines": num_upload_pipelines,
               "total_times": total_times,
               "time_per_kilo_variant": time_per_kilo_variant}
    return render(request, 'upload/view_upload_stats_detail.html', context)


def view_upload_step_detail(request, upload_step_id: int):
    upload_step = get_object_or_404(UploadStep, pk=upload_step_id)
    upload_step.upload_pipeline.file_upload.check_can_view(request.user)
    return render_ajax_view(request, 'upload/upload_step.html', {
        "upload_step": upload_step
    })


def view_upload_pipeline(request, upload_pipeline_id):
    upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
    file_upload = upload_pipeline.file_upload
    file_upload.check_can_view(request.user)

    filename = file_upload.get_filename()
    file_exists = filename and os.path.exists(filename)
    allow_retry_import = (file_upload.user == request.user or request.user.is_superuser) and file_exists

    if not file_exists:
        status = messages.WARNING
        messages.add_message(request, status, "File does not exist on disk, cannot reload", extra_tags='import-message')

    step_order, step_start_end_lines = upload_stats.get_step_order_and_step_start_end_lines(upload_pipeline)
    more_warning_or_error_details = False
    for vii in upload_pipeline.get_errors(hide_accepted=False):
        msg = f"Error: {vii.message}"
        messages.add_message(request, messages.ERROR, msg, extra_tags='import-message')
        more_warning_or_error_details |= vii.has_more_details

    for vii in upload_pipeline.get_warnings(hide_accepted=False):
        msg = f"Warning: {vii.message}"
        messages.add_message(request, messages.WARNING, msg, extra_tags='import-message')
        more_warning_or_error_details |= vii.has_more_details

    context = {
        'upload_pipeline': upload_pipeline,
        'has_upload_steps': upload_pipeline.uploadstep_set.exists(),
        "allow_retry_import": allow_retry_import and settings.UPLOAD_ENABLED,
        "more_warning_or_error_details": more_warning_or_error_details,
        # this data is redundant if there isn't multiple runs for the same step
        "step_total_stats": upload_stats.get_step_total_stats(upload_pipeline, only_if_multiple_runs=True),
        "step_order": list(step_order),
        "step_start_end_lines": step_start_end_lines,
    }
    context.update(_get_basic_uploaded_file_context(file_upload))

    context_name = UPLOADED_FILE_CONTEXT.get(upload_pipeline.file_type)
    if context_name:
        try:
            context[context_name] = get_upload_data_for_uploaded_file(upload_pipeline.file_upload)
        except Exception:
            pass

    return render(request, 'upload/view_upload_pipeline.html', context)


def view_upload_pipeline_warnings_and_errors(request, upload_pipeline_id):
    upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
    upload_pipeline.file_upload.check_can_view(request.user)

    skipped_contigs = VCFSkippedContigs.objects.filter(upload_step__upload_pipeline=upload_pipeline).first()

    contigs_import = get_build_contigs()
    has_miv = ModifiedImportedVariant.objects.filter(import_info__upload_step__upload_pipeline=upload_pipeline).exists()
    skipped_annotation = SimpleVCFImportInfo.objects.filter(upload_step__upload_pipeline=upload_pipeline,
                                                            type=SimpleVCFImportInfo.ANNOTATION_SKIPPED).first()

    for vii in upload_pipeline.get_errors(hide_accepted=False):
        msg = f"Error: {vii.message}"
        messages.add_message(request, messages.ERROR, msg, extra_tags='import-message')

    for vii in upload_pipeline.get_warnings(hide_accepted=False):
        msg = f"Warning: {vii.message}"
        messages.add_message(request, messages.WARNING, msg, extra_tags='import-message')

    context = {
        'upload_pipeline': upload_pipeline,
        'skipped_contigs': skipped_contigs,
        "has_modified_imported_variants": has_miv,
        "skipped_annotation": skipped_annotation,
        'contigs_import': contigs_import,
    }
    return render(request, 'upload/view_upload_pipeline_warnings_and_errors.html', context)


@require_POST
def upload_retry_import(request, upload_pipeline_id):
    upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
    upload_pipeline = retry_upload_pipeline(upload_pipeline)

    msg = 'Attempting to re-importing file'
    status = messages.INFO
    messages.add_message(request, status, msg, extra_tags='import-message')

    return HttpResponseRedirect(reverse("view_upload_pipeline", kwargs={"upload_pipeline_id": upload_pipeline.pk}))


@require_POST
def accept_vcf_import_info_tag(request, vcf_import_info_id):
    vii = VCFImportInfo.objects.get_subclass(pk=vcf_import_info_id)
    vcf_id = vii.upload_step.upload_pipeline.file_upload.uploadedvcf.vcf.pk
    VCF.get_for_user(request.user, vcf_id)  # Permission check
    vii.accepted_date = timezone.now()
    vii.save()

    return JsonResponse({})


class DownloadUploadedFile(PathDownloadView):
    @cached_property
    def file_upload(self):
        file_upload_id = self.kwargs["pk"]
        file_upload = get_object_or_404(FileUpload, pk=file_upload_id)
        upload_data = get_upload_data_for_uploaded_file(file_upload)
        data = upload_data.get_data()
        # TODO: use check_can_view once everything implements GuardianPermissionsMixin
        if not data.can_view(self.request.user):
            raise PermissionDenied(f"You do not have permission to access: {data}")
        return file_upload

    def get_mimetype(self):
        """ Firefox 86 downloads XX.vcf.gz as XX.vcf.vcf - so provide mimetype to force .gz extension
            @see https://stackoverflow.com/a/65596550/295724 """
        mimetype = super().get_mimetype()
        filename = self.get_path()
        if filename.endswith(".gz"):
            mimetype = "application/gzip"
        return mimetype

    def get_path(self):
        return self.file_upload.get_filename()
