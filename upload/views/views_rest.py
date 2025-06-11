import logging
from http import HTTPStatus
from typing import Optional

from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt

from rest_framework.views import APIView

from library.file_utils import file_or_filename_md5sum
from upload.models import UploadedFile
from upload.views.views import handle_file_upload


@method_decorator(csrf_exempt, name='dispatch')
class APIFileUploadView(APIView):
    """ Re-implemented uploads in DRF so we can use API tokens for all client work """
    def post(self, request, *args, **kwargs):
        try:
            django_uploaded_file = request.FILES['file']
            path = request.query_params.get("path")
            force = request.query_params.get("force", False)

            response_data = {}
            existing_uploaded_file: Optional[UploadedFile] = None
            if path and not force:
                existing_uploaded_file = self._get_existing_uploaded_file(path, django_uploaded_file.path)

            if existing_uploaded_file:
                response_data["message"] = "Existing path/hash uploaded file found"
                uploaded_file = existing_uploaded_file
            else:
                # This is still a Django UploadedFile - even though it's API we save as import_source.WEB_UPLOAD
                # So we know to handle it as UploadedFile rather than 'path' (which is the path on the API client)
                uploaded_file = handle_file_upload(request.user, django_uploaded_file, path=path)

            response_data["uploaded_file_id"] = uploaded_file.pk
            return JsonResponse(response_data)
        except Exception as e:
            logging.error(e)
            return JsonResponse({"error": str(e)}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    @staticmethod
    def _get_existing_uploaded_file(path, django_uploaded_file_path) -> Optional[UploadedFile]:
        if existing_ufs := list(UploadedFile.objects.filter(path=path).order_by("pk")):
            new_hash = file_or_filename_md5sum(django_uploaded_file_path)
            for existing_uf in existing_ufs:
                existing_hash = file_or_filename_md5sum(existing_uf.uploaded_file.path)
                if new_hash == existing_hash:
                    return existing_uf
        return None

