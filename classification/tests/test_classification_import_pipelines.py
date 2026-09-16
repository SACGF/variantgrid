"""An import holding both ordinary and gene-level coordinates runs a pipeline for each (#1835)."""
import os
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from classification.classification_import import (
    API_UPLOAD_NAME,
    GENE_LEVEL_API_UPLOAD_NAME,
    _classification_upload_pipeline,
)
from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from classification.models.classification import ClassificationImport
from classification.tasks.classification_import_process_variants_task import (
    ClassificationImportProcessVariantsTask,
)
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME, GENE_LEVEL_REF, GENE_LEVEL_SVLEN
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion, ImportSource, VariantCoordinate
from upload.models import FileUpload, UploadedClassificationImport, UploadPipeline, UploadStep
from upload.models.models_enums import UploadedFileTypes, UploadStepOrigin, UploadStepTaskType


class TestClassificationImportPipelines(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]

    def _run(self, variant_coordinates) -> list[UploadPipeline]:
        classification_import = ClassificationImport.objects.create(user=self.user,
                                                                    genome_build=GenomeBuild.grch37())
        with patch("classification.classification_import.process_upload_pipeline"), \
                patch("classification.classification_import.write_vcf_from_variant_coordinates"), \
                patch("classification.classification_import.get_contigs_header_lines", return_value=[]):
            _classification_upload_pipeline(classification_import, variant_coordinates, ImportSource.API)
        return list(UploadPipeline.objects.filter(
            file_upload__uploadedclassificationimport__classification_import=classification_import
        ).order_by("pk"))

    @staticmethod
    def _gene_level_coordinate() -> VariantCoordinate:
        return VariantCoordinate(chrom=GENE_LEVEL_CONTIG_NAME, position=644, ref=GENE_LEVEL_REF,
                                 alt="<SPLICE:HGNC:644:V7>", svlen=GENE_LEVEL_SVLEN)

    def test_ordinary_only_runs_one_pipeline(self):
        ordinary = VariantCoordinate.from_explicit_no_svlen("1", 169519049, "T", "C")
        pipelines = self._run([ordinary])
        self.assertEqual([p.file_upload.file_type for p in pipelines],
                         [UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY])

    def test_gene_level_runs_its_own_pipeline(self):
        """ Both pipelines link back to the one ClassificationImport - it used to be a OneToOne """
        ordinary = VariantCoordinate.from_explicit_no_svlen("1", 169519049, "T", "C")
        pipelines = self._run([ordinary, self._gene_level_coordinate()])

        self.assertEqual([p.file_upload.file_type for p in pipelines],
                         [UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY,
                          UploadedFileTypes.GENE_LEVEL_INSERT_VARIANTS_ONLY])
        self.assertEqual([p.file_upload.name for p in pipelines],
                         [API_UPLOAD_NAME, GENE_LEVEL_API_UPLOAD_NAME])
        self.assertEqual(UploadedClassificationImport.objects.count(), 2)

    def test_pipelines_do_not_share_an_input_dir(self):
        """ UploadPipeline.remove_generated_input_file deletes the dir, taking a shared sibling's VCF """
        ordinary = VariantCoordinate.from_explicit_no_svlen("1", 169519049, "T", "C")
        pipelines = self._run([ordinary, self._gene_level_coordinate()])
        input_dirs = {os.path.dirname(p.file_upload.path) for p in pipelines}
        self.assertEqual(len(input_dirs), 2, input_dirs)


class TestLinkInsertedVariants(TestCase):
    """ link_inserted_variants hashes what came back unmatched, and a hash needs a Sequence row (#1835) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]

    def test_uninserted_gene_level_coordinate_does_not_kill_the_whole_import(self):
        """ A gene-level alt is unique to the event, so one that never got inserted has no Sequence
            row and hashing it raised - taking every other record in the import down with it """
        genome_build = GenomeBuild.grch37()
        classification_import = ClassificationImport.objects.create(user=self.user, genome_build=genome_build)
        gbpv = GenomeBuildPatchVersion.get_unspecified_patch_version_for(genome_build)
        allele_info = ImportedAlleleInfo.objects.create(
            imported_genome_build_patch_version=gbpv,
            imported_c_hgvs="METex14skip",
            variant_coordinate=f"{GENE_LEVEL_CONTIG_NAME}:7029-7029 <SPLICE:HGNC:7029:ex14skip>",
            classification_import=classification_import)

        file_upload = UploadedClassificationImport.objects.create(
            classification_import=classification_import,
            file_upload=self._file_upload()).file_upload
        pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        upload_step = UploadStep.objects.create(upload_pipeline=pipeline, name="link", sort_order=0,
                                                task_type=UploadStepTaskType.CELERY,
                                                origin=UploadStepOrigin.IMPORT_TASK_FACTORY)

        ClassificationImportProcessVariantsTask.link_inserted_variants(
            genome_build, classification_import, upload_step)

        allele_info.refresh_from_db()
        self.assertEqual(allele_info.status, ImportedAlleleInfoStatus.FAILED)
        self.assertIn("not inserted", allele_info.message)

    def _file_upload(self) -> FileUpload:
        return FileUpload.objects.create(path="/tmp/gene_level.vcf", name=GENE_LEVEL_API_UPLOAD_NAME,
                                         file_type=UploadedFileTypes.GENE_LEVEL_INSERT_VARIANTS_ONLY,
                                         import_source=ImportSource.API, user=self.user)
