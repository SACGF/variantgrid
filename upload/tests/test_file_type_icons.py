from django.test import SimpleTestCase

from upload.file_type_icons import file_type_icon_html
from upload.models.models_enums import UploadedFileTypes


class FileTypeIconsTest(SimpleTestCase):
    def test_every_file_type_has_an_icon(self):
        """ A new UploadedFileTypes value needs an entry in FILE_TYPE_ICONS, or the upload page shows a blank """
        for file_type in UploadedFileTypes:
            html = file_type_icon_html(file_type.value)
            self.assertIn(f'title="{file_type.label}"', html, file_type)
            self.assertTrue('<i class="fa-' in html or '<use href="#' in html, f"{file_type} has no glyph: {html}")

    def test_vcf_sub_types_are_badged_vcfs(self):
        for file_type in (UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY, UploadedFileTypes.MANUAL_VARIANT_ENTRY,
                          UploadedFileTypes.GENE_LEVEL_INSERT_VARIANTS_ONLY):
            html = file_type_icon_html(file_type.value)
            self.assertIn('href="#file-icon-vcf"', html)
            self.assertIn('class="file-type-badge ', html)
        self.assertNotIn("file-type-badge", file_type_icon_html(UploadedFileTypes.VCF.value))

    def test_unknown_type_is_an_empty_block(self):
        self.assertEqual(file_type_icon_html(None), '<div class="file-type-icon"></div>')
