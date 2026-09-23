from django.conf import settings
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema
from rest_framework.response import Response
from rest_framework.views import APIView

from library.git import Git
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import UploadedFileTypes

# The client contract: each name is a fact about this codebase a client may rely on. Append one in the same
# change as a client-visible feature, and keep it while this endpoint exists.
API_FEATURES = (
    "patients",
    "link_extraction",
    "upload_status",
    "joint_called_vcf_cross_run",
    "upload_metadata",
)

# Import factories the server drives itself, rather than files a client uploads (UploadedFileTypes names, lower case)
INTERNAL_UPLOAD_FILE_TYPES = frozenset({
    "analysis",
    "clinvar",
    "liftover",
    "manual_variant_entry",
    "variant_tags",
    "wiki_gene",
    "wiki_variant",
})


class CapabilitiesView(APIView):
    """ What this deployment accepts, so one API client can work against servers of different ages """

    @extend_schema(
        summary="Features and upload file types this server supports",
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, *args, **kwargs):
        upload_file_types = {UploadedFileTypes(factory.get_uploaded_file_type()).name.lower()
                             for factory in get_import_task_factories()}
        return Response({
            "version": settings.VARIANTGRID_VERSION,
            "git_hash": Git(settings.BASE_DIR).hash,
            "features": list(API_FEATURES),
            "upload_file_types": sorted(upload_file_types - INTERNAL_UPLOAD_FILE_TYPES),
        })
