import os
from hashlib import sha256
from http import HTTPStatus
from typing import Optional

from django.core.exceptions import ObjectDoesNotExist
from django.http import JsonResponse, FileResponse, Http404
from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema
from rest_framework.views import APIView

from analysis.models import AnalysisTemplate
from analysis.tasks.analysis_grid_export_tasks import export_cohort_to_downloadable_file, \
    get_grid_downloadable_file_params_hash
from library.log_utils import report_exc_info
from snpdb.models import Cohort, CachedGeneratedFile
from snpdb.models.models_enums import ImportStatus
from upload.models import UploadedFile
from upload.uploaded_file_type import get_upload_data_for_uploaded_file
from upload.views.views import handle_file_upload, get_upload_status_dict


def _get_uploaded_file_for_user(user, *, uploaded_file_id=None, sha256_hash=None) -> UploadedFile:
    """ Resolve an UploadedFile by pk or by content hash, then run the ownership check. """
    if uploaded_file_id is not None:
        uploaded_file = get_object_or_404(UploadedFile, pk=uploaded_file_id)
    elif sha256_hash:
        # Most-recent matching upload the user can see
        uploaded_file = get_object_or_404(UploadedFile.objects.order_by("-pk"), sha256_hash=sha256_hash)
    else:
        raise Http404("Provide uploaded_file_id or sha256_hash")
    uploaded_file.check_can_view(user)
    return uploaded_file


@method_decorator(csrf_exempt, name='dispatch')
class APIFileUploadView(APIView):
    """ Re-implemented uploads in DRF so we can use API tokens for all client work """

    @extend_schema(
        summary="Upload a file (multipart form field 'file'), returning the uploaded file ID",
        request=OpenApiTypes.BINARY,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, *args, **kwargs):
        try:
            django_uploaded_file = request.FILES['file']
            path = request.query_params.get("path")
            force = request.query_params.get("force", False)

            response_data = {}
            existing_uploaded_file: Optional[UploadedFile] = None
            if not force:
                sha256_hash = self._django_file_sha256(django_uploaded_file)
                existing_uploaded_file = self._get_existing_uploaded_file(request.user, sha256_hash, path)

            if existing_uploaded_file:
                response_data["message"] = "Existing uploaded file with matching hash found"
                uploaded_file = existing_uploaded_file
            else:
                # This is still a Django UploadedFile - even though it's API we save as import_source.WEB_UPLOAD
                # So we know to handle it as UploadedFile rather than 'path' (which is the path on the API client)
                uploaded_file = handle_file_upload(request.user, django_uploaded_file, path=path)

            response_data["uploaded_file_id"] = uploaded_file.pk
            response_data["sha256_hash"] = uploaded_file.sha256_hash
            return JsonResponse(response_data)
        except Exception:
            report_exc_info(request=request)
            return JsonResponse({"error": "Upload failed"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    @classmethod
    def _get_existing_uploaded_file(cls, user, sha256_hash, path=None) -> Optional[UploadedFile]:
        """ Most-recent successfully-processed upload with the same content the user can view.
            Prefers a matching 'path' when supplied. """
        matches = []
        for existing_uf in UploadedFile.objects.filter(sha256_hash=sha256_hash).order_by("-pk"):
            if not existing_uf.can_view(user):
                continue
            # Only re-use uploaded files that have been successfully processed
            upload_data = get_upload_data_for_uploaded_file(existing_uf)
            if upload_data is not None and existing_uf.uploaded_file:
                matches.append(existing_uf)

        if not matches:
            return None
        if path:
            for existing_uf in matches:
                if existing_uf.path == path:
                    return existing_uf
        return matches[0]

    @staticmethod
    def _django_file_sha256(django_uploaded_file) -> str:
        hasher = sha256()
        for chunk in django_uploaded_file.chunks():
            hasher.update(chunk)
        django_uploaded_file.seek(0)  # rewind so handle_file_upload can re-read it
        return hasher.hexdigest()


class APIUploadStatusView(APIView):
    """ Poll import + annotation status for an uploaded file (by id or by sha256 hash) """

    @extend_schema(
        summary="Get import/annotation status for an uploaded file (by id or sha256 hash)",
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, *args, **kwargs):
        uploaded_file = _get_uploaded_file_for_user(
            request.user,
            uploaded_file_id=kwargs.get("uploaded_file_id"),
            sha256_hash=kwargs.get("sha256_hash"))
        return JsonResponse(get_upload_status_dict(uploaded_file))


class APIAnnotatedDownloadView(APIView):
    """ Download the annotated VCF/CSV for an uploaded VCF (by id or sha256 hash).

        Returns 202 with progress while the export is still generating, then streams the file. """
    EXPORT_TYPES = {"vcf", "csv"}

    @extend_schema(
        summary="Download the annotated VCF/CSV for an uploaded VCF (cohort-level export)",
        responses=OpenApiTypes.BINARY,
    )
    def get(self, request, *args, **kwargs):
        export_type = kwargs["export_type"]
        if export_type not in self.EXPORT_TYPES:
            return JsonResponse({"error": f"export_type must be one of {sorted(self.EXPORT_TYPES)}"},
                                status=HTTPStatus.BAD_REQUEST)

        uploaded_file = _get_uploaded_file_for_user(
            request.user,
            uploaded_file_id=kwargs.get("uploaded_file_id"),
            sha256_hash=kwargs.get("sha256_hash"))

        try:
            uploaded_vcf = uploaded_file.uploadedvcf
        except ObjectDoesNotExist:
            uploaded_vcf = None
        vcf = uploaded_vcf.vcf if uploaded_vcf else None
        if vcf is None:
            return JsonResponse({"error": "Uploaded file is not a VCF, or import has not created a VCF yet"},
                                status=HTTPStatus.BAD_REQUEST)
        if vcf.import_status != ImportStatus.SUCCESS:
            import_status_label = ImportStatus(vcf.import_status).label
            return JsonResponse({"error": f"VCF import not complete (import_status={import_status_label})"},
                                status=HTTPStatus.BAD_REQUEST)

        cohort = vcf.cohort
        Cohort.get_for_user(request.user, cohort.pk)  # Permission check

        try:
            AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
        except ValueError as e:
            return JsonResponse({"error": str(e)}, status=HTTPStatus.NOT_IMPLEMENTED)

        params_hash = get_grid_downloadable_file_params_hash(cohort.pk, export_type)
        task = export_cohort_to_downloadable_file.si(cohort.pk, export_type)
        cgf = CachedGeneratedFile.get_or_create_and_launch("export_cohort_to_downloadable_file", params_hash, task)

        if cgf.exception:
            return JsonResponse({"error": cgf.exception}, status=HTTPStatus.INTERNAL_SERVER_ERROR)
        if not cgf.filename:
            return JsonResponse({"status": "generating", "progress": cgf.progress or 0.0},
                                status=HTTPStatus.ACCEPTED)

        basename = os.path.basename(cgf.filename)
        if basename.endswith(".zip"):
            content_type = "application/zip"
        elif basename.endswith(".gz"):
            content_type = "application/gzip"
        else:
            content_type = "application/octet-stream"
        return FileResponse(open(cgf.filename, "rb"), as_attachment=True, filename=basename,
                            content_type=content_type)
