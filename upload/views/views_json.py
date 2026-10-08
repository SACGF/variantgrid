import logging

from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.http.response import HttpResponse, JsonResponse
from django.utils import timezone
from django.views.decorators.http import require_http_methods, require_POST

from analysis.models import AnalysisTemplate
from annotation.models import AnnotationRun, VariantAnnotationVersion
from eventlog.models import create_event
from library.django_utils.file_uploads import filepond_process_response, filepond_upload_receive
from library.enums.log_level import LogLevel
from library.log_utils import log_traceback
from snpdb.models import VCF
from snpdb.models.models_enums import ImportStatus
from upload import upload_processing
from upload.models import (
    FileUpload,
    ImportSource,
    ProcessingStatus,
    UploadedFileTypes,
    UploadPipeline,
    VCFImportInfo,
)
from upload.models.models_enums import VCFImportInfoSeverity
from upload.upload_metadata import (
    UploadMetadataError,
    get_metadata_keys_for_file_type,
    validate_upload_metadata,
)
from upload.uploaded_file_type import get_uploaded_file_type, get_url_and_data_for_uploaded_file_data


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
        data["requires_user_input"] = upload_data.requires_user_input
    return data


def get_remaining_annotation_runs(uploaded_vcf, genome_build) -> int:
    """ Unfinished runs of the build's active annotation version that cover this VCF's variants. Only the
        active version counts: a run left behind on a retired version never executes, and would otherwise
        report every later upload as still annotating (and hide its downloads) forever """
    max_variant_id = uploaded_vcf.max_variant_id
    if max_variant_id is None:
        # VCF not fully imported yet, so highest known variant is unknown - no remaining runs to report
        return 0
    variant_annotation_version = VariantAnnotationVersion.latest(genome_build)
    if variant_annotation_version is None:
        return 0
    ar_qs = AnnotationRun.get_active_runs(genome_build).filter(annotation_range_lock__version=variant_annotation_version)
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
        "warnings": [],
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
    # What the VCF page asks the user to accept, eg REF bases that mostly mismatch the build (#2030)
    data["warnings"] = [{"severity": VCFImportInfoSeverity(vii.severity).label, "message": vii.message}
                        for vii in upload_pipeline.get_vcf_import_info()]

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


@require_POST
def accept_vcf_import_info_tag(request, vcf_import_info_id):
    vii = VCFImportInfo.objects.get_subclass(pk=vcf_import_info_id)
    vcf_id = vii.upload_step.upload_pipeline.file_upload.uploadedvcf.vcf.pk
    VCF.get_for_user(request.user, vcf_id)  # Permission check
    vii.accepted_date = timezone.now()
    vii.save()

    return JsonResponse({})
