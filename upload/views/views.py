from datetime import date, timedelta
from django.conf import settings
from django.contrib import messages
from django.core.exceptions import PermissionDenied
from django.http.response import HttpResponseRedirect, JsonResponse
from django.shortcuts import render, get_object_or_404
from django.urls.base import reverse
from django.utils import timezone
from django.utils.timesince import timesince
from django.views.decorators.http import require_POST
from django_downloadview import PathDownloadView
from jfu.http import upload_receive, UploadResponse, JFUResponse
import json
import logging
import os

from lazy import lazy

from annotation.views import get_build_contigs
from eventlog.models import create_event
from library.enums.log_level import LogLevel
from library.log_utils import log_traceback
from snpdb.models import VCF
from upload import forms, upload_processing, upload_stats
from upload.models import UploadPipeline, UploadedFile, ProcessingStatus, UploadedFileTypes, \
    UploadSettings, ImportSource, UploadStep, VCFSkippedContigs, \
    VCFImportInfo, SimpleVCFImportInfo, ModifiedImportedVariant, TimeFilterMethod
from upload.uploaded_file_type import get_upload_data_for_uploaded_file, \
    get_uploaded_file_type, get_url_and_data_for_uploaded_file_data, \
    retry_upload_pipeline


UPLOADED_FILE_CONTEXT = {UploadedFileTypes.VCF: "uploaded_vcf",
                         UploadedFileTypes.GENE_LIST: "uploaded_gene_list",
                         UploadedFileTypes.CUFFDIFF: 'uploaded_cuffdiff'}


def get_icon_for_uploaded_file_status(status):
    THUMBNAILS = {ProcessingStatus.CREATED: 'queued.png',
                  ProcessingStatus.PROCESSING: 'loading.gif',
                  ProcessingStatus.ERROR: 'cross.png',
                  ProcessingStatus.SUCCESS: 'tick.png',
                  ProcessingStatus.TERMINATED_EARLY: 'warning.png'}
    icon = THUMBNAILS.get(status)
    if icon:
        return os.path.join(settings.STATIC_URL, 'icons', icon)
    return None


def uploadedfile_dict(uploaded_file) -> dict:
    try:
        size = uploaded_file.uploaded_file.size
    except:
        size = None

    data_url, upload_data = get_url_and_data_for_uploaded_file_data(uploaded_file)
    time_since = timesince(uploaded_file.created)
    file_type = None
    if uploaded_file.file_type:
        file_type = UploadedFileTypes(uploaded_file.file_type).label

    data = {
        'uploaded_file_id': uploaded_file.pk,
        'name': uploaded_file.name,
        'size': size,
        'user': uploaded_file.user.get_full_name(),
        'time_since': f"{time_since} ago",
        'file_type': file_type,
        'file_type_code': uploaded_file.file_type,
        'data_url': data_url,
        'deleteUrl': reverse('jfu_delete', kwargs={'pk': uploaded_file.pk}),
        'deleteType': 'POST',
    }

    if upload_data:
        try:
            data["upload_data"] = upload_data.get_upload_context()
        except:
            pass

    if not uploaded_file.file_type:
        data['error'] = f'Could not determine how to read file: "{uploaded_file.name}"'

    try:
        upload_pipeline = UploadPipeline.objects.get(uploaded_file=uploaded_file)
        data["upload_pipeline_id"] = upload_pipeline.pk
        status = upload_pipeline.status
        url = reverse('view_upload_pipeline', kwargs={'upload_pipeline_id': upload_pipeline.pk})
    except:
        status = ProcessingStatus.ERROR
        url = reverse('view_uploaded_file', kwargs={'uploaded_file_id': uploaded_file.pk})

    data['processing_status'] = status
    data['status_image'] = get_icon_for_uploaded_file_status(status)
    data["url"] = url
    return data


def handle_file_upload(user, django_uploaded_file, import_source=ImportSource.WEB_UPLOAD, path=None) -> UploadedFile:
    original_filename = django_uploaded_file._name
    kwargs = {
        "name": original_filename,
        "uploaded_file": django_uploaded_file,
        "import_source": import_source,
        "user": user,
        "path": path,
    }
    uploaded_file = UploadedFile.objects.create(**kwargs)
    # Save 1st to actually create file (need to open handling unicode)
    uploaded_file.file_type = get_uploaded_file_type(uploaded_file, original_filename)
    uploaded_file.save()

    if uploaded_file.file_type:
        upload_processing.process_uploaded_file(uploaded_file)
    return uploaded_file

@require_POST
def jfu_upload(request):
    # The assumption here is that jQuery File Upload
    # has been configured to send files one at a time.
    # If multiple files can be uploaded simulatenously,
    # 'file' may be a list of files.

    if not settings.UPLOAD_ENABLED:
        raise PermissionDenied("Uploads are currently disabled (settings.UPLOAD_ENABLED=False)")

    try:
        message = "jfu_upload POST had empty FILES - this is often due to running out of disk space for Nginx temp file storage."
        if not request.FILES:
            create_event(request.user, message, severity=LogLevel.ERROR)
            raise ValueError(message)

        django_uploaded_file = upload_receive(request)
        uploaded_file = handle_file_upload(request.user, django_uploaded_file)

        file_dict = uploadedfile_dict(uploaded_file)
    except Exception as e:
        logging.error(e)
        log_traceback()
        file_dict = {"error": str(e)}

    logging.debug("views.upload() END")
    return UploadResponse(request, file_dict)


@require_POST
def jfu_upload_delete(request, pk):
    try:
        instance = UploadedFile.objects.get(pk=pk)
        if not (request.user.is_superuser or request.user == instance.user):
            msg = f"You don't own uploaded file {pk}"
            raise PermissionDenied(msg)

        data = {'uploaded_file_id': instance.pk}
        try:
            if os.path.exists(instance.uploaded_file.path):
                os.unlink(instance.uploaded_file.path)
        except:  # Someone may have deleted MEDIA_ROOT file already - causing uploaded_file to error here
            pass
        instance.delete()
    except UploadedFile.DoesNotExist:
        data = False

    return JFUResponse(request, data)


def get_file_dicts_list(upload_settings):
    qs = UploadedFile.objects.order_by("-created")
    if not upload_settings.user.is_superuser:
        qs = qs.filter(user=upload_settings.user)
    if not upload_settings.show_all:
        qs = qs.filter(visible=True)
    if upload_settings.time_filter_method == TimeFilterMethod.DAYS:
        start_date = date.today() - timedelta(days=upload_settings.time_filter_value)
        qs = qs.filter(created__gte=start_date)
    elif upload_settings.time_filter_method == TimeFilterMethod.RECORDS:
        qs = qs[:upload_settings.time_filter_value]

    file_dicts = []
    for uploaded_file in qs:
        file_dicts.append(uploadedfile_dict(uploaded_file))
    file_dicts = list(reversed(file_dicts))  # JFU adds most recent at the end
    return file_dicts


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
    context = {'existing_files': file_dicts_list,
               'form': form,
               "upload_enabled": settings.UPLOAD_ENABLED}
    return render(request, 'upload/upload.html', context)


def view_uploaded_file(request, uploaded_file_id):
    uploaded_file = get_object_or_404(UploadedFile, pk=uploaded_file_id)
    uploaded_file.check_can_view(request.user)
    context = {'uploaded_file': uploaded_file}
    return render(request, 'upload/view_uploaded_file.html', context)


def view_upload_stats(request):
    DAYS = 30

    end_date = date.today() - timedelta(days=DAYS)
    ups_qs = UploadStep.objects.filter(end_date__gte=end_date, name__isnull=False)
    us = list(map(json.dumps, upload_stats.get_upload_stats(ups_qs, max_step_names=10)))
    num_upload_pipelines, total_times, time_per_kilo_variant = us

    context = {"days": DAYS,
               "num_upload_pipelines": num_upload_pipelines,
               "total_times": total_times,
               "time_per_kilo_variant": time_per_kilo_variant}
    return render(request, 'upload/view_upload_stats.html', context)


def view_upload_pipeline(request, upload_pipeline_id):
    upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
    uploaded_file = upload_pipeline.uploaded_file
    uploaded_file.check_can_view(request.user)

    filename = uploaded_file.get_filename()
    file_exists = os.path.exists(filename)
    allow_retry_import = (uploaded_file.user == request.user or request.user.is_superuser) and file_exists

    if not file_exists:
        status = messages.WARNING
        messages.add_message(request, status, "File does not exist on disk, cannot reload", extra_tags='import-message')

    if upload_pipeline.status == ProcessingStatus.ERROR:
        msg = f"Import Error: This file failed to import due to: {upload_pipeline.progress_status}. "
        status = messages.ERROR
        messages.add_message(request, status, msg, extra_tags='import-message')

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

    context = {'upload_pipeline': upload_pipeline,
               'has_upload_steps': upload_pipeline.uploadstep_set.exists(),
               "allow_retry_import": allow_retry_import and settings.UPLOAD_ENABLED,
               "more_warning_or_error_details": more_warning_or_error_details,
               "step_total_stats": upload_stats.get_step_total_stats(upload_pipeline),
               "step_order": list(step_order),
               "step_start_end_lines": step_start_end_lines}

    context_name = UPLOADED_FILE_CONTEXT.get(upload_pipeline.file_type)
    if context_name:
        try:
            context[context_name] = get_upload_data_for_uploaded_file(upload_pipeline.uploaded_file)
        except Exception:
            pass

    return render(request, 'upload/view_upload_pipeline.html', context)


def view_upload_pipeline_warnings_and_errors(request, upload_pipeline_id):
    upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
    upload_pipeline.uploaded_file.check_can_view(request.user)

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
    vcf_id = vii.upload_step.upload_pipeline.uploaded_file.uploadedvcf.vcf.pk
    VCF.get_for_user(request.user, vcf_id)  # Permission check
    vii.accepted_date = timezone.now()
    vii.save()

    return JsonResponse({})


class DownloadUploadedFile(PathDownloadView):
    @lazy
    def uploaded_file(self):
        uploaded_file_id = self.kwargs["pk"]
        uploaded_file = get_object_or_404(UploadedFile, pk=uploaded_file_id)
        upload_data = get_upload_data_for_uploaded_file(uploaded_file)
        data = upload_data.get_data()
        # TODO: use check_can_view once everything implements GuardianPermissionsMixin
        if not data.can_view(self.request.user):
            raise PermissionDenied(f"You do not have permission to access: {data}")
        return uploaded_file

    def get_mimetype(self):
        """ Firefox 86 downloads XX.vcf.gz as XX.vcf.vcf - so provide mimetype to force .gz extension
            @see https://stackoverflow.com/a/65596550/295724 """
        mimetype = super().get_mimetype()
        filename = self.get_path()
        if filename.endswith(".gz"):
            mimetype = "application/gzip"
        return mimetype

    def get_path(self):
        return self.uploaded_file.get_filename()
