import os
from hashlib import md5
from http import HTTPStatus
from typing import IO, Optional

from django.conf import settings
from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema
from rest_framework.views import APIView

from library.log_utils import report_exc_info
from upload.models import UploadedFile
from upload.uploaded_file_type import get_upload_data_for_uploaded_file
from upload.views.views import handle_file_upload


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
            if path and not force:
                existing_uploaded_file = self._get_existing_uploaded_file(request.user,
                                                                          path, django_uploaded_file)

            if existing_uploaded_file:
                response_data["message"] = "Existing path/hash uploaded file found"
                uploaded_file = existing_uploaded_file
            else:
                # This is still a Django UploadedFile - even though it's API we save as import_source.WEB_UPLOAD
                # So we know to handle it as UploadedFile rather than 'path' (which is the path on the API client)
                uploaded_file = handle_file_upload(request.user, django_uploaded_file, path=path)

            response_data["uploaded_file_id"] = uploaded_file.pk
            return JsonResponse(response_data)
        except Exception:
            report_exc_info(request=request)
            return JsonResponse({"error": "Upload failed"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    @classmethod
    def _get_existing_uploaded_file(cls, user, path, django_uploaded_file) -> Optional[UploadedFile]:
        existing_ufs = list(UploadedFile.objects.filter(path=path).order_by("pk"))
        if not existing_ufs:
            return None

        new_hash = cls._md5sum(django_uploaded_file)
        for existing_uf in existing_ufs:
            existing_uf.check_can_view(user)
            # Only look at uploaded files that have been successfully processed
            upload_data = get_upload_data_for_uploaded_file(existing_uf)
            if upload_data is not None and existing_uf.uploaded_file:
                if new_hash == cls._uploaded_file_md5sum(existing_uf):
                    return existing_uf
        return None

    @classmethod
    def _uploaded_file_md5sum(cls, uploaded_file: UploadedFile) -> str:
        """ md5sum of the stored file, ensuring the resolved path stays within UPLOAD_DIR """
        upload_dir = os.path.realpath(settings.UPLOAD_DIR)
        file_path = os.path.realpath(uploaded_file.uploaded_file.path)
        if not file_path.startswith(upload_dir + os.sep):
            raise ValueError(f"Uploaded file path '{file_path}' is not within UPLOAD_DIR")
        with open(file_path, "rb") as f:
            return cls._md5sum(f)

    @staticmethod
    def _md5sum(f: IO) -> str:
        # MD5 is used here only as a fast file-equivalency checksum (detecting identical
        # file contents), NOT for any security/cryptographic purpose, so collision
        # resistance is not required.
        m = md5()
        for chunk in iter(lambda: f.read(8192), b''):
            m.update(chunk)
        return m.hexdigest()
