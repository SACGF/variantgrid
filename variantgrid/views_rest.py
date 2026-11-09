import re

from django.conf import settings
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework.views import APIView

from library.git import Git
from upload.enums import UploadedFileTypes
from upload.import_task_factories.import_task_factory import get_import_task_factories

# The client contract: each name is a fact about this codebase a client may rely on. Append one in the same
# change as a client-visible feature, and keep it while this endpoint exists.
API_FEATURES = (
    "patients",
    "link_extraction",
    "link_extractions_bulk",
    "upload_status",
    "joint_called_vcf_cross_run",
    "upload_metadata",
    "cram_alignment_files",  # A sequencing file's 'bam_file' path can be a CRAM
    "alignment_files",  # SequencingFile, SingleSampleVCF and QC take an 'alignment_files' list (BAM/CRAM)
    "annotated_download_compression",  # Annotated CSV download takes ?compression=zip|none
    "symbolic_alts",  # VCF import stores <DEL>/<DUP>/<INV>/<CNV> alts (with SVLEN/END) as structural variants
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


def release_version(version: str) -> str:
    """ The release tag a git-describe version counts from ('vg4.0-12-gc174556' -> 'vg4.0', 'vg4-gc174556' -> 'vg4'):
        the commit count would pin the exact commit as surely as the hash """
    return re.sub(r"(-\d+)?-g[0-9a-f]+$", "", version)


class CapabilitiesView(APIView):
    """ What this deployment accepts, so one API client can work against servers of different ages.
        Public, so a client can probe a server before it has credentials, but only an authenticated caller
        learns the exact commit the server runs (#2732) """
    permission_classes = [AllowAny]

    @extend_schema(
        summary="Features and upload file types this server supports",
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, *args, **kwargs):
        upload_file_types = {UploadedFileTypes(factory.get_uploaded_file_type()).name.lower()
                             for factory in get_import_task_factories()}
        data = {
            "version": settings.VARIANTGRID_VERSION,
            "features": list(API_FEATURES),
            "upload_file_types": sorted(upload_file_types - INTERNAL_UPLOAD_FILE_TYPES),
        }
        if request.user.is_authenticated:
            data["git_hash"] = Git(settings.BASE_DIR).hash
        else:
            data["version"] = release_version(data["version"])
        return Response(data)
