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
from classification.models.classification import ClassificationImport
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME, GENE_LEVEL_REF, GENE_LEVEL_SVLEN
from snpdb.models import GenomeBuild, ImportSource, VariantCoordinate
from upload.models import UploadedClassificationImport, UploadPipeline
from upload.models.models_enums import UploadedFileTypes


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
