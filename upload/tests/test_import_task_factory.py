import os

from django.conf import settings
from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase, override_settings

from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import UploadedFileTypes
from upload.upload_processing import get_vcf_file_type

TSO500_PAIR_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500", "ExampleSample_2600000001")

GENE_LEVEL_FILE_TYPES = {
    UploadedFileTypes.DRAGEN_TSO500_ALL_FUSIONS,
    UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT,
    UploadedFileTypes.GENE_LEVEL_CNV_VCF,
    UploadedFileTypes.GENE_LEVEL_SPLICE_VCF,
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


class ProcessVCFFileTypeTest(TestCase):
    """ process_vcf_file (import_vcf, seqauto) picks a VCF's file type off its contents, as an upload does """

    def test_each_vcf_takes_the_path_its_contents_call_for(self):
        user = User.objects.get_or_create(username="testuser")[0]
        expected = {
            "ExampleSample_RNA_2600000001B/ExampleSample_RNA_2600000001B_SpliceVariants.vcf":
                UploadedFileTypes.GENE_LEVEL_SPLICE_VCF,
            "ExampleSample_DNA_2600000001C/ExampleSample_DNA_2600000001C.cnv.vcf": UploadedFileTypes.GENE_LEVEL_CNV_VCF,
            "ExampleSample_DNA_2600000001C/ExampleSample_DNA_2600000001C.hard-filtered.vcf": UploadedFileTypes.VCF,
        }
        for filename, file_type in expected.items():
            self.assertEqual(file_type, get_vcf_file_type(user, os.path.join(TSO500_PAIR_DIR, filename)), filename)
