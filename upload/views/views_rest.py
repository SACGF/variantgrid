from http import HTTPStatus

from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt

from rest_framework.views import APIView

from snpdb.models import ImportSource
from upload.views.views import handle_file_upload


@method_decorator(csrf_exempt, name='dispatch')
class APIFileUploadView(APIView):
    """ Re-implemented uploads in DRF so we can use API tokens for all client work """
    def post(self, request, *args, **kwargs):
        try:
            django_uploaded_file = request.FILES['file']
            path = request.query_params.get("path")
            uploaded_file = handle_file_upload(request.user, django_uploaded_file,
                                               import_source=ImportSource.API, path=path)
            return JsonResponse({"uploaded_file_id": uploaded_file.pk})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=HTTPStatus.INTERNAL_SERVER_ERROR)
