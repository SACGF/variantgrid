import os
import tempfile

from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models import ImportSource
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline, UploadStep
from upload.tasks.vcf.genotype_vcf_tasks import ImportCreateVCFModelForGenotypeVCFTask
from upload.vcf.vcf_import import VCFSamplesRemovedException

VCF_HEADER = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=249250621>
"""


class TestVCFReloadSamples(TestCase):
    """ A VCF reloaded in place from a file whose sample columns have changed (#1191) """

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        user = User.objects.get_or_create(username='testuser')[0]
        vcf_filename = self._write_vcf("original.vcf", ["A", "B"])
        file_upload = FileUpload.objects.create(path=vcf_filename, import_source=ImportSource.COMMAND_LINE,
                                                user=user, file_type=UploadedFileTypes.VCF,
                                                metadata={"genome_build": "GRCh37"})
        self.upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        self._create_vcf_model(vcf_filename)
        self.vcf = self.upload_pipeline.uploadedvcf.vcf

    def _write_vcf(self, name, sample_names) -> str:
        filename = os.path.join(self.temp_dir.name, name)
        columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *sample_names]
        with open(filename, "w") as f:
            f.write(VCF_HEADER + "\t".join(columns) + "\n")
        return filename

    def _create_vcf_model(self, vcf_filename):
        upload_step = UploadStep.objects.create(upload_pipeline_id=self.upload_pipeline.pk,  # fresh, as in celery
                                                input_filename=vcf_filename, sort_order=0)
        ImportCreateVCFModelForGenotypeVCFTask.process_items(upload_step)

    def _packed_and_sort_by_vcf_sample_name(self) -> dict:
        return {cs.sample.vcf_sample_name: (cs.cohort_genotype_packed_field_index, cs.sort_order)
                for cs in self.vcf.cohort.cohortsample_set.select_related("sample")}

    def test_reload_with_reordered_and_added_samples(self):
        self._create_vcf_model(self._write_vcf("reloaded.vcf", ["B", "A", "C"]))

        self.assertEqual({"B": (0, 0), "A": (1, 1), "C": (2, 2)}, self._packed_and_sort_by_vcf_sample_name())
        self.vcf.refresh_from_db()
        self.assertEqual(3, self.vcf.genotype_samples)
        self.assertEqual(3, self.vcf.cohort.cohort_genotype_collection.num_samples)

    def test_reload_missing_a_sample_fails(self):
        with self.assertRaises(VCFSamplesRemovedException):
            self._create_vcf_model(self._write_vcf("reloaded.vcf", ["B"]))
        self.assertEqual({"A": (0, 0), "B": (1, 1)}, self._packed_and_sort_by_vcf_sample_name())
