"""Thin server-side helpers for the FilePond client.

FilePond's server contract is described at https://pqina.nl/filepond/docs/api/server/.
We implement only the small subset of endpoints we use: ``process`` (uploads)
and ``revert`` (delete just-uploaded). Existing files are seeded into FilePond
through the ``files`` instance option using the ``"local"`` server-file type,
so the ``load`` endpoint is not required.

serve_attachment returns an uploaded attachment: images and PDFs display in the browser, anything else
downloads, so an uploaded page is never rendered on the site.
"""
import os

from django.http import FileResponse, HttpResponse

INLINE_ATTACHMENT_EXTENSIONS = {'.png', '.gif', '.jpg', '.jpeg', '.pdf'}


def filepond_upload_receive(request, field_name="filepond"):
    """Return the Django UploadedFile FilePond posted under ``field_name``.

    Raises ValueError when nothing was posted (typically a nginx temp-storage
    issue rather than a client bug).
    """
    if not request.FILES:
        raise ValueError(
            "Upload POST had empty FILES - this is often due to running out of "
            "disk space for Nginx temp file storage."
        )
    django_uploaded_file = request.FILES.get(field_name)
    if django_uploaded_file is None:
        # FilePond's default field name is ``filepond`` but a caller may
        # supply a different name; fall back to the first file present.
        django_uploaded_file = next(iter(request.FILES.values()))
    return django_uploaded_file


def filepond_process_response(file_id) -> HttpResponse:
    """Return the plain-text response FilePond expects from ``process``.

    FilePond reads the response body as the opaque server id used later for
    revert/load/restore. We use the model PK as the id.
    """
    return HttpResponse(str(file_id), content_type="text/plain")




def serve_attachment(filename: str) -> FileResponse:
    """ Images and PDFs display inline, every other type downloads """
    _, ext = os.path.splitext(filename)
    as_attachment = ext.lower() not in INLINE_ATTACHMENT_EXTENSIONS
    return FileResponse(open(filename, "rb"), as_attachment=as_attachment)
