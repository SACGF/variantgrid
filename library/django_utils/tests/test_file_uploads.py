import os
import tempfile

from django.test import SimpleTestCase

from library.django_utils.file_uploads import serve_attachment


class ServeAttachmentTest(SimpleTestCase):
    def test_only_images_and_pdfs_display_inline(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            for basename, inline in [("scan.PNG", True), ("report.pdf", True), ("page.html", False), ("pic.svg", False)]:
                filename = os.path.join(tmp_dir, basename)
                with open(filename, "w") as f:
                    f.write("x")
                response = serve_attachment(filename)
                response.close()
                disposition = response.get("Content-Disposition", "")
                self.assertEqual(inline, not disposition.startswith("attachment"), basename)
