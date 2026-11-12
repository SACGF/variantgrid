"""
one_off_backfill_t2t_gnomad_ref_swaps (#2788): fills gnomAD on the variants of a version flagged
backfilled_gnomad_ref_swaps=False from the settings' gnomAD file, via bcftools.
"""
import copy
import gzip
import os
import shutil
import subprocess
import tempfile
import unittest

from django.conf import settings
from django.core.management import CommandError, call_command
from django.test import TestCase, override_settings

from annotation.fake_data import get_fake_vep_version, retire_seeded_annotation_version
from annotation.management.commands.one_off_backfill_t2t_gnomad_ref_swaps import GNOMAD_INFO_FIELDS
from annotation.models import AnnotationVersion, VariantAnnotation, VariantAnnotationVersion
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from genes.enums import AnnotationConsortium
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

GNOMAD_VERSION = "4.1.test"
_INFO_TYPES = {"AC": "Integer", "AN": "Integer", "AF": "Float"}
GNOMAD_HEADER = "##fileformat=VCFv4.2\n##contig=<ID=NC_000001.11>\n" + "".join(
    f'##INFO=<ID={info_field},Number=A,Type={_INFO_TYPES.get(info_field[:2], "String")},Description="">\n'
    for info_field in GNOMAD_INFO_FIELDS.values()
) +"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


class OneOffBackfillGnomADRefSwapsTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        retire_seeded_annotation_version(cls.genome_build)
        kwargs = get_fake_vep_version(cls.genome_build, AnnotationConsortium.ENSEMBL, 5)
        kwargs["gnomad"] = GNOMAD_VERSION
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs, status=VariantAnnotationVersion.Status.ACTIVE,
                                                          backfilled_gnomad_ref_swaps=False)
        AnnotationVersion.latest(cls.genome_build, validate=False)

        # [0] is a ref swap site the file now has, [1] is in no gnomAD, [2] already had its AF
        cls.variants = [slowly_create_test_variant("1", 100 + i, "A", "G", cls.genome_build) for i in range(3)]
        annotation_range_lock = AnnotationRangeLock.objects.create(version=cls.vav, min_variant=cls.variants[0],
                                                                   max_variant=cls.variants[-1],
                                                                   count=len(cls.variants))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock)
        partition = cls.vav.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, partition):
            for variant in cls.variants:
                VariantAnnotation.objects.create(version=cls.vav, variant=variant, annotation_run=annotation_run)
        VariantAnnotation.objects.filter(version=cls.vav, variant=cls.variants[2]).update(gnomad_af=0.25)

    def _settings(self, gnomad_filename: str) -> dict:
        annotation = copy.deepcopy(settings.ANNOTATION)
        annotation[self.genome_build.name]["vep_config"]["gnomad4"] = gnomad_filename
        return {"ANNOTATION": annotation, "ANNOTATION_VCF_DUMP_DIR": os.path.dirname(gnomad_filename)}

    @staticmethod
    def _write_gnomad(gnomad_filename: str, records: list[str]):
        vcf_filename = gnomad_filename.removesuffix(".bgz")
        with open(vcf_filename, "w", encoding="utf-8") as f:
            f.write(GNOMAD_HEADER)
            f.writelines(record + "\n" for record in records)
        subprocess.run(["bgzip", "-c", vcf_filename], stdout=open(gnomad_filename, "wb"), check=True)
        subprocess.run(["tabix", "-p", "vcf", gnomad_filename], check=True)

    def _gnomad_afs(self) -> list:
        qs = VariantAnnotation.objects.filter(version=self.vav, variant__in=self.variants).order_by("variant_id")
        return list(qs.values_list("gnomad_af", flat=True))

    @unittest.skipUnless(shutil.which(settings.BCFTOOLS_COMMAND) and shutil.which("bgzip"), "needs bcftools")
    def test_fills_missing_gnomad_and_sets_flag(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            gnomad_filename = os.path.join(tmp_dir, f"gnomad{GNOMAD_VERSION}_GRCh38_combined_af.vcf.bgz")
            self._write_gnomad(gnomad_filename, [
                "NC_000001.11\t100\t.\tA\tG\t.\tPASS\tAC=500;AF=0.5;AN=1000;grpmax=nfe",
                "NC_000001.11\t102\t.\tA\tG\t.\tPASS\tAC=300;AF=0.3;AN=1000;grpmax=nfe",
            ])
            with override_settings(**self._settings(gnomad_filename)):
                call_command("one_off_backfill_t2t_gnomad_ref_swaps", "--dry-run", stdout=open(os.devnull, "w"))
                self.vav.refresh_from_db()
                self.assertFalse(self.vav.backfilled_gnomad_ref_swaps)
                self.assertEqual(self._gnomad_afs(), [None, None, 0.25])

                call_command("one_off_backfill_t2t_gnomad_ref_swaps", stdout=open(os.devnull, "w"))

        annotation = VariantAnnotation.objects.get(version=self.vav, variant=self.variants[0])
        self.assertEqual((annotation.gnomad_ac, annotation.gnomad_an, annotation.gnomad_popmax), (500, 1000, "NFE"))
        # Only variants with no gnomAD values are dumped - [2]'s AF stays as annotated
        self.assertEqual(self._gnomad_afs(), [0.5, None, 0.25])
        self.vav.refresh_from_db()
        self.assertTrue(self.vav.backfilled_gnomad_ref_swaps)

    def test_refuses_when_settings_gnomad_differs_from_version(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            gnomad_filename = os.path.join(tmp_dir, "gnomad4.1.other_GRCh38_combined_af.vcf.bgz")
            with gzip.open(gnomad_filename, "wt") as f:
                f.write(GNOMAD_HEADER)
            with override_settings(**self._settings(gnomad_filename)):
                with self.assertRaisesRegex(CommandError, "4.1.other"):
                    call_command("one_off_backfill_t2t_gnomad_ref_swaps")
        self.vav.refresh_from_db()
        self.assertFalse(self.vav.backfilled_gnomad_ref_swaps)
