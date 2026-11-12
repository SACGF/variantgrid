import gzip
import os
import shutil
import tempfile
from io import StringIO
from unittest.mock import patch

import cyvcf2
from django.conf import settings
from django.core.management import call_command
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.enums import (
    AnnotationStatus,
    VariantAnnotationPipelineType,
    VEPSkippedReason,
)
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
    retire_seeded_annotation_version,
)
from annotation.models import (
    AnnotationVersion,
    VariantAnnotation,
    VariantAnnotationVersion,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.models.models_annotation_run import (
    AnnotationPipelineVersion,
    AnnotationRangeLock,
    AnnotationRun,
)
from annotation.vep_annotation import (
    VEPConfig,
    get_vep_skipped_variants_filename,
    get_vep_version_from_vcf,
    vep_dict_to_variant_annotation_version_kwargs,
)
from annotation.vep_warnings import SKIPPED_VARIANTS_FILE_HEADER_PREFIX
from genes.enums import AnnotationConsortium
from genes.fake_data import create_fake_gene_version
from library.django_utils.django_partition import bulk_create_in_partition
from library.genomics.vcf_enums import VariantClass
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf

COMMAND = "one_off_reannotate_converted_symbolic_variants"
MODULE = f"annotation.management.commands.{COMMAND}"
TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
# STRUCTURAL_VARIANT pipeline output: <INV> 1.18 Mb (variant 202), <DEL> 14.5 kb (203), <DUP> 6 kb (205)
SV_VEP_OUTPUT = os.path.join(TEST_DATA_DIR, "test_columns_version3_grch37_sv.vep_annotated.vcf")
ANNOTSV_TSV = os.path.join(TEST_DATA_DIR, "annotsv", "test_grch37_sv.annotated.tsv")
INV, DEL, DUP = 202, 203, 205
STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL_VARIANT = VariantAnnotationPipelineType.STRUCTURAL_VARIANT
ANNOTSV = VariantAnnotationPipelineType.ANNOTSV


def fake_sv_vep(_vav, vcf_dump_filename, vcf_annotated_filename):
    """ VEP's output for the dumped variants, taken from the fixture """
    dumped_ids = {v.INFO["variant_id"] for v in cyvcf2.VCF(vcf_dump_filename)}
    with open(SV_VEP_OUTPUT) as f_in, gzip.open(vcf_annotated_filename, "wt") as f_out:
        for line in f_in:
            if line.startswith("#") or int(line.split("variant_id=")[1].split(";")[0]) in dumped_ids:
                f_out.write(line)
    with open(get_vep_skipped_variants_filename(vcf_annotated_filename), "w") as f:
        f.write(f"{SKIPPED_VARIANTS_FILE_HEADER_PREFIX} {vcf_dump_filename}]\n")


@override_settings(**{**get_fake_annotation_settings_dict(columns_version=3),
                      "ANNOTATION_VEP_SV_CONSERVATION_PYBIGWIG_ENABLED": False,
                      "ANNOTATION_ANNOTSV_ENABLED": False,
                      "ANNOTATION_VEP_SV_MAX_SIZE": 1_000_000})  # leaves the inversion out of the dump
@patch(f"{MODULE}.run_sv_vep", side_effect=fake_sv_vep)
class OneOffReannotateConvertedSymbolicVariantsTest(TestCase):
    """ One finished range lock over the fixture's SVs: the inversion and duplication were converted in place
        (STANDARD rows), the deletion was annotated by the STRUCTURAL_VARIANT run """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        seeded_vav = get_fake_annotation_version(cls.grch37).variant_annotation_version
        gene_annotation_release = seeded_vav.gene_annotation_release
        retire_seeded_annotation_version(cls.grch37)
        vep_dict = get_vep_version_from_vcf(SV_VEP_OUTPUT)
        kwargs = vep_dict_to_variant_annotation_version_kwargs(VEPConfig(cls.grch37), vep_dict)
        kwargs.update(genome_build=cls.grch37, gene_annotation_release=gene_annotation_release,
                      status=VariantAnnotationVersion.Status.ACTIVE)
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)
        slowly_create_loci_and_variants_for_vcf(cls.grch37, SV_VEP_OUTPUT, get_variant_id_from_info=True)

        cls.lock = AnnotationRangeLock.objects.create(version=cls.vav, min_variant_id=INV, max_variant_id=DUP, count=3)
        cls.standard_run = cls._finished_run(STANDARD)
        cls.sv_run = cls._finished_run(STRUCTURAL_VARIANT)
        cls.gene = create_fake_gene_version(cls.grch37, "ENSG_FAKE_2786", "FAKE2786",
                                            AnnotationConsortium.ENSEMBL).gene
        cls._annotate(cls.standard_run, [INV, DUP])
        cls._annotate(cls.sv_run, [DEL])

    @classmethod
    def _finished_run(cls, pipeline_type, **kwargs) -> AnnotationRun:
        return AnnotationRun.objects.create(annotation_range_lock=cls.lock, pipeline_type=pipeline_type,
                                            upload_end=timezone.now(), **kwargs)

    @classmethod
    def _annotate(cls, run, variant_ids):
        for klass, table, kwargs in [
            (VariantAnnotation, VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION, {}),
            (VariantTranscriptAnnotation, VariantAnnotationVersion.TRANSCRIPT_ANNOTATION, {}),
            (VariantGeneOverlap, VariantAnnotationVersion.VARIANT_GENE_OVERLAP, {"gene": cls.gene}),
        ]:
            records = [klass(version=cls.vav, variant_id=variant_id, annotation_run=run, **kwargs)
                       for variant_id in variant_ids]
            bulk_create_in_partition(cls.vav, klass, table, records)

    def setUp(self):
        dump_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, dump_dir)
        dump_dir_settings = self.settings(ANNOTATION_VCF_DUMP_DIR=dump_dir)
        dump_dir_settings.enable()
        self.addCleanup(dump_dir_settings.disable)

    def _call(self, **kwargs) -> str:
        out = StringIO()
        call_command(COMMAND, variant_annotation_version=self.vav.pk, stdout=out, **kwargs)
        return out.getvalue()

    def _variant_annotation(self, variant_id) -> VariantAnnotation:
        return VariantAnnotation.objects.get(version=self.vav, variant_id=variant_id)

    def _assert_untouched(self, mock_vep):
        mock_vep.assert_not_called()
        self.assertEqual(self._variant_annotation(INV).annotation_run, self.standard_run)
        self.assertEqual(self._variant_annotation(DUP).annotation_run, self.standard_run)

    def test_reannotates_into_the_sv_run(self, mock_vep):
        output = self._call()
        self.assertIn("2 variants with STANDARD annotation now routed to STRUCTURAL_VARIANT", output)
        mock_vep.assert_called_once()

        dup = self._variant_annotation(DUP)
        self.assertEqual(dup.annotation_run, self.sv_run)
        self.assertEqual(dup.variant_class, VariantClass.DUPLICATION)
        inv = self._variant_annotation(INV)
        self.assertEqual(inv.annotation_run, self.sv_run)
        self.assertEqual(inv.vep_skipped_reason, VEPSkippedReason.TOO_LONG)
        for klass in [VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap]:
            self.assertFalse(klass.objects.filter(annotation_run=self.standard_run).exists(), klass.__name__)
        self.assertEqual(self._variant_annotation(DEL).annotation_run, self.sv_run)

        self.assertEqual(AnnotationRangeLock.objects.filter(version=self.vav).count(), 1)
        self.assertEqual(AnnotationRun.objects.filter(annotation_range_lock__version=self.vav).count(), 2)
        self.sv_run.refresh_from_db()
        self.assertIsNone(self.sv_run.vcf_annotated_filename)

        self.assertIn("0 variants", self._call())

    def test_dry_run_changes_nothing(self, mock_vep):
        output = self._call(dry_run=True)
        self.assertIn(f"2 variants -> STRUCTURAL_VARIANT run {self.sv_run.pk}", output)
        self._assert_untouched(mock_vep)

    def test_lock_without_finished_sv_run(self, mock_vep):
        AnnotationRun.objects.filter(pk=self.sv_run.pk).update(upload_end=None, status=AnnotationStatus.ERROR,
                                                               error_exception="failed")
        self.assertIn("no finished STRUCTURAL_VARIANT run", self._call(dry_run=True))
        self.assertIn("0 in 0 range locks", self._call())
        self._assert_untouched(mock_vep)

    def test_leaves_historical_versions(self, mock_vep):
        VariantAnnotationVersion.objects.filter(pk=self.vav.pk).update(
            status=VariantAnnotationVersion.Status.HISTORICAL)
        self.assertIn("left as they are", self._call())
        self._assert_untouched(mock_vep)

    @override_settings(ANNOTATION_ANNOTSV_ENABLED=True)
    @patch(f"{MODULE}.run_annotsv", return_value=ANNOTSV_TSV)
    def test_annotsv_over_finished_annotsv_run(self, mock_annotsv, _mock_vep):
        pipeline_version = AnnotationPipelineVersion.objects.create(
            pipeline_type=ANNOTSV, genome_build=self.grch37, code_version="3.5", data_version="3.5",
            status=AnnotationPipelineVersion.Status.ACTIVE)
        annotsv_run = self._finished_run(ANNOTSV, pipeline_version=pipeline_version)
        self._call()
        mock_annotsv.assert_called_once()
        self.assertEqual(mock_annotsv.call_args.args[0], annotsv_run)
        self.assertEqual(self._variant_annotation(DUP).annotsv_acmg_score, -0.5)
