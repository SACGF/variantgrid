from django.test import SimpleTestCase, override_settings

from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import UploadedFileTypes

GENE_LEVEL_FILE_TYPES = {
    UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS,
    UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT,
    UploadedFileTypes.GENE_LEVEL_CNV_VCF,
    UploadedFileTypes.GENE_LEVEL_INSERT_VARIANTS_ONLY,
}


class ImportTaskFactoryEnabledTest(SimpleTestCase):

    @staticmethod
    def _file_types() -> set[str]:
        return {itf.get_uploaded_file_type() for itf in get_import_task_factories()}

    def test_gene_level_file_types_withdrawn_when_disabled(self):
        self.assertTrue(GENE_LEVEL_FILE_TYPES.issubset(self._file_types()))
        with override_settings(VARIANT_GENE_LEVEL_ENABLED=False):
            file_types = self._file_types()
        self.assertFalse(GENE_LEVEL_FILE_TYPES & file_types)
        self.assertIn(UploadedFileTypes.VCF, file_types)
