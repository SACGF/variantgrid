from django.test import TestCase
from classification.models import AttachmentFileType


class TestVariantClassificationModels(TestCase):

    def setUp(self):
        pass

    def test_file_types(self):
        IMAGES = ["foo.jpg",
                  "foo.jpeg",
                  "FOO.JPG",
                  "foo.bar.baz.png"]
        TEXT = ["png.TXT"]
        OTHER = ["foo.", "foo"]
        PDF = ["foo.pdf", 'foo.PDF']

        for image in IMAGES:
            self.assertEqual(AttachmentFileType.IMAGE, AttachmentFileType.get_type_for_file(image))

        for text in TEXT:
            self.assertEqual(AttachmentFileType.TEXT, AttachmentFileType.get_type_for_file(text))

        for pdf in PDF:
            self.assertEqual(AttachmentFileType.PDF, AttachmentFileType.get_type_for_file(pdf))

        for other in OTHER:
            self.assertEqual(AttachmentFileType.OTHER, AttachmentFileType.get_type_for_file(other))
