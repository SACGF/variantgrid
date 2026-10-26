import os
from hashlib import sha256
from http import HTTPStatus
from typing import Optional

from django.core.exceptions import ObjectDoesNotExist
from django.http import FileResponse, Http404, JsonResponse
from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import OpenApiParameter, extend_schema
from rest_framework.views import APIView

from analysis.models import AnalysisTemplate
from analysis.tasks.analysis_grid_export_tasks import (
    export_cohort_to_downloadable_file,
    get_export_annotation_version,
    get_grid_downloadable_file_params_hash,
    get_zipped_csv_filename,
)
from library.log_utils import report_exc_info
from snpdb.enums import ImportStatus
from snpdb.models import CachedGeneratedFile, Cohort
from upload.models import FileUpload
from upload.upload_metadata import (
    UploadMetadataError,
    get_metadata_keys_for_file_type,
    validate_upload_metadata,
)
from upload.uploaded_file_type import get_upload_data_for_uploaded_file
from upload.views.views_json import EmptyUploadError, get_upload_status_dict, handle_file_upload


def _get_file_upload_for_user(user, *, file_upload_id=None, sha256_hash=None) -> FileUpload:
    """ Resolve a FileUpload by pk or by content hash, then run the ownership check. """
    if file_upload_id is not None:
        file_upload = get_object_or_404(FileUpload, pk=file_upload_id)
    elif sha256_hash:
        # Most-recent matching upload the user can see
        file_uploads = FileUpload.objects.filter(sha256_hash=sha256_hash).order_by("-pk")
        if not user.is_superuser:
            file_uploads = file_uploads.filter(user=user)
        file_upload = file_uploads.first()
        if file_upload is None:
            raise Http404(f"No upload with sha256_hash={sha256_hash}")
    else:
        raise Http404("Provide file_upload_id or sha256_hash")
    file_upload.check_can_view(user)
    return file_upload


@method_decorator(csrf_exempt, name='dispatch')
class APIFileUploadView(APIView):
    """ Re-implemented uploads in DRF so we can use API tokens for all client work.

        The optional 'path' query param is a SeqAuto-only hint: on SeqAuto deployments it links the
        upload to a pre-registered backend sequencing VCF at that server path. General upload / annotate
        / download clients should omit it and identify uploads by the returned file_upload_id and/or
        sha256_hash (content hash, which is stable across servers). 'uploaded_file_id' is a retained
        alias for 'file_upload_id'.

        Every other query param is upload metadata - facts about the file it doesn't carry itself, whose
        accepted keys depend on the detected file type (a VCF takes 'genome_build' and 'source'). An
        unknown key or an unresolvable build is a 400. @see upload.upload_metadata """

    RESERVED_QUERY_PARAMS = {"path", "force"}

    @extend_schema(
        summary="Upload a file (multipart form field 'file'), returning the uploaded file ID",
        description="Optional query params: 'path' (SeqAuto backend-link hint only - omit for general "
                    "uploads) and 'force' (skip sha256 de-duplication). Any other query param is treated "
                    "as upload metadata for the detected file type - a VCF accepts 'genome_build' and "
                    "'source'; an unknown key is a 400. Send a build's own name ('GRCh37') rather than "
                    "an alias ('hg19'). Returns file_upload_id and sha256_hash. "
                    "'uploaded_file_id' is a retained alias for 'file_upload_id'.",
        request=OpenApiTypes.BINARY,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, *args, **kwargs):
        try:
            django_uploaded_file = request.FILES['file']
            path = request.query_params.get("path")
            force = request.query_params.get("force", "").lower() in {"1", "true", "yes"}
            metadata = {k: v for k, v in request.query_params.items()
                        if k not in self.RESERVED_QUERY_PARAMS}

            response_data = {}
            existing_file_upload: Optional[FileUpload] = None
            if not force:
                sha256_hash = self._django_file_sha256(django_uploaded_file)
                existing_file_upload = self._get_existing_file_upload(request.user, sha256_hash, path, metadata)

            if existing_file_upload:
                response_data["message"] = "Existing uploaded file with matching hash found"
                file_upload = existing_file_upload
            else:
                # This is still a Django UploadedFile - even though it's API we save as import_source.WEB_UPLOAD
                # So we know to handle it as a web upload rather than 'path' (which is the path on the API client)
                file_upload = handle_file_upload(request.user, django_uploaded_file, path=path,
                                                 metadata=metadata)

            response_data["file_upload_id"] = file_upload.pk
            response_data["uploaded_file_id"] = file_upload.pk  # deprecated alias
            response_data["sha256_hash"] = file_upload.sha256_hash
            return JsonResponse(response_data)
        except (UploadMetadataError, EmptyUploadError) as e:
            return JsonResponse({"error": str(e)}, status=HTTPStatus.BAD_REQUEST)
        except Exception:
            report_exc_info(request=request)
            return JsonResponse({"error": "Upload failed"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    @classmethod
    def _get_existing_file_upload(cls, user, sha256_hash, path=None, metadata=None) -> Optional[FileUpload]:
        """ Most-recent successfully-processed upload with the same content and metadata the user can view.
            Prefers a matching 'path' when supplied. Different metadata (eg a corrected genome_build) is a
            new upload, not the old one. """
        matches = []
        for existing_uf in FileUpload.objects.filter(sha256_hash=sha256_hash).order_by("-pk"):
            if not existing_uf.can_view(user):
                continue
            metadata_keys = get_metadata_keys_for_file_type(existing_uf.file_type)
            if validate_upload_metadata(metadata, metadata_keys) != (existing_uf.metadata or {}):
                continue
            # Only re-use uploaded files that have been successfully processed
            upload_data = get_upload_data_for_uploaded_file(existing_uf)
            if upload_data is not None and existing_uf.file_field:
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
        description="Returns file_upload_id and sha256_hash along with import/annotation progress. "
                    "'uploaded_file_id' is a retained alias for 'file_upload_id'.",
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, *args, **kwargs):
        file_upload = _get_file_upload_for_user(
            request.user,
            file_upload_id=kwargs.get("file_upload_id"),
            sha256_hash=kwargs.get("sha256_hash"))
        return JsonResponse(get_upload_status_dict(file_upload))


class APIAnnotatedDownloadView(APIView):
    """ Download the annotated VCF/CSV for an uploaded VCF (by id or sha256 hash).

        Returns 202 with progress while the export is still generating, then streams the file.
        ?compression=zip|none (CSV only, default zip): the VCF is always gzipped """
    EXPORT_TYPES = {"vcf", "csv"}
    CSV_COMPRESSIONS = {"zip", "none"}

    @extend_schema(
        summary="Download the annotated VCF/CSV for an uploaded VCF (cohort-level export)",
        parameters=[OpenApiParameter("compression", str, enum=["zip", "none"], default="zip",
                                     description="CSV only: 'none' for a plain .csv")],
        responses=OpenApiTypes.BINARY,
    )
    def get(self, request, *args, **kwargs):
        export_type = kwargs["export_type"]
        if export_type not in self.EXPORT_TYPES:
            return JsonResponse({"error": f"export_type must be one of {sorted(self.EXPORT_TYPES)}"},
                                status=HTTPStatus.BAD_REQUEST)
        compression = request.GET.get("compression", "zip")
        if compression not in self.CSV_COMPRESSIONS:
            return JsonResponse({"error": f"compression must be one of {sorted(self.CSV_COMPRESSIONS)}"},
                                status=HTTPStatus.BAD_REQUEST)

        file_upload = _get_file_upload_for_user(
            request.user,
            file_upload_id=kwargs.get("file_upload_id"),
            sha256_hash=kwargs.get("sha256_hash"))

        try:
            uploaded_vcf = file_upload.uploadedvcf
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

        annotation_version = get_export_annotation_version(cohort)
        params_hash = get_grid_downloadable_file_params_hash(cohort.pk, export_type, annotation_version.pk)
        task = export_cohort_to_downloadable_file.si(cohort.pk, export_type, annotation_version.pk)
        cgf = CachedGeneratedFile.get_or_create_and_launch("export_cohort_to_downloadable_file", params_hash, task)

        if cgf.exception:
            return JsonResponse({"error": cgf.exception}, status=HTTPStatus.INTERNAL_SERVER_ERROR)
        if not cgf.filename:
            return JsonResponse({"status": "generating", "progress": cgf.progress or 0.0},
                                status=HTTPStatus.ACCEPTED)

        filename = cgf.filename
        if export_type == "csv" and compression == "zip" and not filename.endswith(".zip"):
            filename = get_zipped_csv_filename(cgf)
        basename = os.path.basename(filename)
        if basename.endswith(".zip"):
            content_type = "application/zip"
        elif basename.endswith(".gz"):
            content_type = "application/gzip"
        elif basename.endswith(".csv"):
            content_type = "text/csv"
        else:
            content_type = "application/octet-stream"
        return FileResponse(open(filename, "rb"), as_attachment=True, filename=basename,
                            content_type=content_type)
