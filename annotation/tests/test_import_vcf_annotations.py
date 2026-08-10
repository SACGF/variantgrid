"""
#1701: the import lane must fail a run whose annotated VCF is short, rather than importing it and
writing vep_skipped_reason=UNKNOWN rows over the records VEP lost.
"""
import os
from unittest import mock

from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.fake_annotation import get_fake_annotation_settings_dict
from annotation.models import VariantAnnotation
from annotation.models.models import AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import VEPSkippedReason
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import (
    VEPConfig,
    get_vep_version_from_vcf,
    vep_dict_to_variant_annotation_version_kwargs,
)
from annotation.vep_warnings import TruncatedVEPOutputError
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf

TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
TEST_ANNOTATION_VCF = os.path.join(TEST_DATA_DIR, "test_columns_version1_grch37.vep_annotated.vcf")

# Two records VEP said it deliberately dropped - see annotation/tests/test_vep_warnings.py for the format
VEP_WARNINGS_TWO_SKIPS = (
    "WARNING: line 5 skipped (1 200 . GG A . . variant_id=1...): "
    "Length of reference allele (GG length 2) does not match coordinates 200-200\n"
    "WARNING: line 6 skipped (1 300 . A <FOO> . . ...): a brand new reason\n"
)


@override_settings(**get_fake_annotation_settings_dict(columns_version=1))
class TruncatedVEPOutputTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        vep_config = VEPConfig(cls.genome_build)
        vep_dict = get_vep_version_from_vcf(TEST_ANNOTATION_VCF)
        kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_dict)
        kwargs["genome_build"] = cls.genome_build
        cls.vav, _ = VariantAnnotationVersion.objects.get_or_create(**kwargs)
        slowly_create_loci_and_variants_for_vcf(cls.genome_build, TEST_ANNOTATION_VCF,
                                                get_variant_id_from_info=True)
        with open(TEST_ANNOTATION_VCF) as f:
            cls.num_records = sum(1 for line in f if not line.startswith("#"))

    def _annotation_run(self, **kwargs) -> AnnotationRun:
        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(self.vav)
        annotation_range_lock.save()
        # task_id: the import lane always holds the execution lock, and the failure path logs against it
        return AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                            vcf_annotated_filename=TEST_ANNOTATION_VCF,
                                            task_id="test-task", **kwargs)

    def _import(self, annotation_run):
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

    def test_complete_file_imports(self):
        annotation_run = self._annotation_run(dump_count=self.num_records)
        self._import(annotation_run)
        self.assertEqual(annotation_run.annotated_count, self.num_records)
        self.assertEqual(VariantAnnotation.objects.filter(annotation_run=annotation_run).count(),
                         self.num_records)

    def test_short_file_with_no_warnings_raises(self):
        annotation_run = self._annotation_run(dump_count=self.num_records + 2)
        with self.assertRaises(TruncatedVEPOutputError):
            self._import(annotation_run)
        # The lost records must NOT have been written off as skipped
        self.assertFalse(VariantAnnotation.objects.filter(annotation_run=annotation_run,
                                                          vep_skipped_reason=VEPSkippedReason.UNKNOWN).exists())

    def test_short_file_with_matching_warnings_imports(self):
        annotation_run = self._annotation_run(dump_count=self.num_records + 2,
                                              vep_warnings=VEP_WARNINGS_TWO_SKIPS)
        # handle_vep_skipped writes its skip rows through the version partition, which under UNIT_TEST is
        # empty (bulk_create goes to the base table) - so it would re-offer the records just imported
        with mock.patch("annotation.vcf_files.import_vcf_annotations.handle_vep_skipped") as mock_skipped:
            self._import(annotation_run)
        mock_skipped.assert_called_once()
        self.assertEqual(annotation_run.annotated_count, self.num_records)

    def test_no_dump_count_skips_check(self):
        """ Legacy rows predating dump_count (and every existing test) are left alone """
        annotation_run = self._annotation_run()
        self._import(annotation_run)
        self.assertEqual(annotation_run.annotated_count, self.num_records)

    def test_skipped_variants_file_preferred_over_warnings(self):
        skipped_variants_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, "skipped_variants.tsv")
        os.makedirs(settings.ANNOTATION_VCF_DUMP_DIR, exist_ok=True)
        with open(skipped_variants_filename, "w") as f:
            f.write("[VEP skipped the following variants from dump.vcf]\n")
            f.write("Line 5    \t1 200 . GG A . . variant_id=1  \tsomething VEP dropped\n")
        self.addCleanup(os.remove, skipped_variants_filename)

        # One row in the file, but the warnings text mentions two - the file wins, so this is 1 short
        annotation_run = self._annotation_run(dump_count=self.num_records + 2,
                                              vep_warnings=VEP_WARNINGS_TWO_SKIPS,
                                              vep_skipped_variants_filename=skipped_variants_filename)
        with self.assertRaises(TruncatedVEPOutputError):
            self._import(annotation_run)
