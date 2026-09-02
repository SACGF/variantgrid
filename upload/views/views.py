import json
import os
from datetime import timedelta
from functools import cached_property

from django.conf import settings
from django.contrib import messages
from django.core.exceptions import PermissionDenied
from django.http.response import HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse
from django.utils import timezone
from django.views.decorators.http import require_POST
from django_downloadview import PathDownloadView

from annotation.views import get_build_contigs
from library.utils.django_utils import render_ajax_view
from upload import forms, upload_stats
from upload.models import (
    FileUpload,
    ModifiedImportedVariant,
    SimpleVCFImportInfo,
    UploadedFileTypes,
    UploadPipeline,
    UploadSettings,
    UploadStep,
    VCFSkippedContigs,
)
from upload.uploaded_file_type import (
    get_import_tasks_by_extension,
    get_upload_data_for_uploaded_file,
    retry_upload_pipeline,
)
from upload.views.views_json import _get_basic_uploaded_file_context, get_file_dicts_list

UPLOADED_FILE_CONTEXT = {UploadedFileTypes.VCF: "uploaded_vcf",
                         UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS: "uploaded_vcf",
                         UploadedFileTypes.GENE_LIST: "uploaded_gene_list"}


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
