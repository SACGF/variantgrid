import os
from unittest.mock import patch

from django.conf import settings
from django.core.management import call_command
from django.test import TestCase
from django.test.utils import override_settings

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.fake_annotation import get_fake_annotation_settings_dict
from annotation.models import VariantAnnotation, VariantAnnotationPipelineType
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import AnnotationRun, VariantAnnotationVersion
from annotation.sv_conservation import (
    ConservationTrack,
    conservation_sidecar_filename,
    write_conservation_sidecar,
)
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import (
    VEPConfig,
    get_vep_version_from_vcf,
    vep_dict_to_variant_annotation_version_kwargs,
)
from library.genomics.vcf_enums import VariantClass
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf


# CNV stuff only comes in with columns version >= 3
@override_settings(**get_fake_annotation_settings_dict(columns_version=3))
class TestAnnotationVCFCNV(TestCase):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
    TEST_ANNOTATION_VCF_GRCH37 = os.path.join(TEST_DATA_DIR, "test_columns_version3_grch37_sv.vep_annotated.vcf")
    TEST_ANNOTATION_VCF_GRCH38 = os.path.join(TEST_DATA_DIR, "test_columns_version3_grch38_sv.vep_annotated.vcf")

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.variant_annotation_versions_by_build = {}

        VCFS = {
            "GRCh37": cls.TEST_ANNOTATION_VCF_GRCH37,
            "GRCh38": cls.TEST_ANNOTATION_VCF_GRCH38,
        }

        for genome_build_name, vcf in VCFS.items():
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
            vep_config = VEPConfig(genome_build)
            vep_dict = get_vep_version_from_vcf(vcf)
            kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_dict)
            kwargs["genome_build"] = genome_build
            vav, created = VariantAnnotationVersion.objects.get_or_create(**kwargs)
            if not created:
                vav.truncate_related_objects()
            cls.variant_annotation_versions_by_build[genome_build_name] = vav

            slowly_create_loci_and_variants_for_vcf(genome_build, vcf, get_variant_id_from_info=True)

    def test_import_variant_annotations_grch37(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh37')
        vav = self.variant_annotation_versions_by_build[genome_build.name]

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH37)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

        # Verify a few fields - got these from vep_reader --pick
        # 3	127535894	.	G	<INV>	.	.	SVTYPE=INV;SVLEN=1184482;variant_id=202
        va = VariantAnnotation.objects.get(variant_id=202)
        self.assertEqual(va.variant_class, VariantClass.INVERSION)
        self.assertEqual(va.impact, PathogenicityImpact.MODIFIER)
        # #1571 - a ranged inv is HGVS from coordinates alone, so 1.1Mb is no obstacle
        self.assertEqual("NC_000003.11:g.127535894_128720376inv", va.hgvs_g)

        # 17	41236500	.	G	<DEL>	.	.	SVTYPE=DEL;SVLEN=14500;variant_id=203
        va = VariantAnnotation.objects.get(variant_id=203)
        self.assertEqual(va.variant_class, VariantClass.DELETION)
        self.assertEqual(va.impact, PathogenicityImpact.HIGH)

        # 21	36418200	.	G	<DUP>	.	.	SVTYPE=DUP;SVLEN=6000;variant_id=205
        va = VariantAnnotation.objects.get(variant_id=205)
        self.assertEqual(va.variant_class, VariantClass.DUPLICATION)
        self.assertEqual(va.impact, PathogenicityImpact.MODIFIER)

    def test_import_variant_annotations_grch38(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh38')
        vav = self.variant_annotation_versions_by_build[genome_build.name]

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH38)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

        # Verify a few fields
        # 3	128486000	.	C	<DUP>	.	.	SVLEN=9951;SVTYPE=DUP;variant_id=101
        va = VariantAnnotation.objects.get(variant_id=101)
        self.assertEqual(va.variant_class, VariantClass.DUPLICATION)
        self.assertEqual(va.impact, PathogenicityImpact.MODIFIER)

        # 21	35041808	.	G	<DEL>	.	.	SVLEN=10000;SVTYPE=DEL;variant_id=102
        va = VariantAnnotation.objects.get(variant_id=102)
        self.assertEqual(va.variant_class, VariantClass.DELETION)
        self.assertEqual(va.impact, PathogenicityImpact.HIGH)
        # #1571 - the padding base at POS is excluded, so the interval opens at POS + 1
        self.assertEqual("NC_000021.9:g.35041809_35051808del", va.hgvs_g)


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class TestAnnotationVCFCNV4(TestAnnotationVCFCNV):
    """ columns_version 4 = VEP 115 + gnomAD 4.1 + masked SpliceAI. SV conservation _max columns
        (phastCons/phyloP) are computed with pyBigWig instead of VEP --custom bigwig (#1657): the
        annotate stage writes a sidecar TSV next to the annotated VCF and the import path merges it
        into the same columns. Here we write that sidecar directly (scoring accuracy is covered by
        annotation.tests.test_sv_conservation, which needs the real bigWig data) and assert the import
        populates the columns from it. """
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
    TEST_ANNOTATION_VCF_GRCH37 = os.path.join(TEST_DATA_DIR, "test_columns_version4_grch37_sv.vep_annotated.vcf")
    TEST_ANNOTATION_VCF_GRCH38 = os.path.join(TEST_DATA_DIR, "test_columns_version4_grch38_sv.vep_annotated.vcf")

    # {variant_id: {db_column: value}} - written to the conservation sidecar before import.
    CONSERVATION_GRCH37 = {
        202: {"phastcons_100_way_vertebrate": 1.0, "phastcons_46_way_mammalian": 1.0,
              "phylop_100_way_vertebrate": 9.873, "phylop_46_way_mammalian": 2.894},
        203: {"phastcons_100_way_vertebrate": 1.0, "phastcons_46_way_mammalian": 1.0,
              "phylop_100_way_vertebrate": 6.180, "phylop_46_way_mammalian": 2.873},
    }
    CONSERVATION_GRCH38 = {
        101: {"phastcons_100_way_vertebrate": 1.0, "phastcons_30_way_mammalian": 1.0,
              "phylop_100_way_vertebrate": 9.556, "phylop_30_way_mammalian": 1.251},
        102: {"phastcons_100_way_vertebrate": 1.0, "phastcons_30_way_mammalian": 1.0,
              "phylop_100_way_vertebrate": 7.272, "phylop_30_way_mammalian": 1.312},
    }

    def setUp(self):
        super().setUp()
        self._sidecars = []

    def tearDown(self):
        for sidecar in self._sidecars:
            if os.path.exists(sidecar):
                os.remove(sidecar)
        super().tearDown()

    def _write_conservation_sidecar(self, vcf_filename, conservation):
        columns = sorted({c for values in conservation.values() for c in values})
        tracks = [ConservationTrack(name=c, path="", db_column=c) for c in columns]
        sidecar = conservation_sidecar_filename(vcf_filename)
        write_conservation_sidecar(sidecar, conservation, tracks)
        self._sidecars.append(sidecar)

    def test_import_variant_annotations_grch37(self):
        # Write the pyBigWig sidecar next to the annotated VCF so the import path picks it up.
        self._write_conservation_sidecar(self.TEST_ANNOTATION_VCF_GRCH37, self.CONSERVATION_GRCH37)
        super().test_import_variant_annotations_grch37()

        va = VariantAnnotation.objects.get(variant_id=202)
        self.assertAlmostEqual(va.phastcons_100_way_vertebrate, 1.0)
        self.assertAlmostEqual(va.phastcons_46_way_mammalian, 1.0)
        self.assertAlmostEqual(va.phylop_100_way_vertebrate, 9.873, places=3)
        self.assertAlmostEqual(va.phylop_46_way_mammalian, 2.894, places=3)

        va = VariantAnnotation.objects.get(variant_id=203)
        self.assertAlmostEqual(va.phylop_100_way_vertebrate, 6.180, places=3)
        self.assertAlmostEqual(va.phylop_46_way_mammalian, 2.873, places=3)

    def test_import_variant_annotations_grch38(self):
        self._write_conservation_sidecar(self.TEST_ANNOTATION_VCF_GRCH38, self.CONSERVATION_GRCH38)
        super().test_import_variant_annotations_grch38()

        va = VariantAnnotation.objects.get(variant_id=101)
        self.assertAlmostEqual(va.phastcons_100_way_vertebrate, 1.0)
        self.assertAlmostEqual(va.phastcons_30_way_mammalian, 1.0)
        self.assertAlmostEqual(va.phylop_100_way_vertebrate, 9.556, places=3)
        self.assertAlmostEqual(va.phylop_30_way_mammalian, 1.251, places=3)

        va = VariantAnnotation.objects.get(variant_id=102)
        self.assertAlmostEqual(va.phylop_100_way_vertebrate, 7.272, places=3)
        self.assertAlmostEqual(va.phylop_30_way_mammalian, 1.312, places=3)

    def test_import_without_conservation_sidecar_fails(self):
        """ #1657: no sidecar means the conservation columns would import as nulls indistinguishable from
            real values, so the run has to fail (visible, retryable) rather than land wrong data. """
        genome_build = GenomeBuild.get_name_or_alias('GRCh37')
        vav = self.variant_annotation_versions_by_build[genome_build.name]
        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH37)
        with self.assertRaises(FileNotFoundError):
            import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

    def _import_grch37_run(self) -> AnnotationRun:
        genome_build = GenomeBuild.get_name_or_alias('GRCh37')
        vav = self.variant_annotation_versions_by_build[genome_build.name]
        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH37)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)
        return annotation_run

    def test_backfill_sv_conservation(self):
        """ #1657: the backfill re-scores rows the pyBigWig stage wrote - correcting a value the
            off-by-one window got wrong and filling one a failed sidecar left null, while leaving a
            value that already agrees (and any run VEP scored itself) alone. """
        self._write_conservation_sidecar(self.TEST_ANNOTATION_VCF_GRCH37, self.CONSERVATION_GRCH37)
        annotation_run = self._import_grch37_run()
        annotation_run.sv_conservation_pybigwig = True  # our stage wrote these columns, not VEP
        annotation_run.save()

        va = VariantAnnotation.objects.get(variant_id=202)
        va.phylop_100_way_vertebrate = 3.5      # what the off-by-one window scored
        va.phastcons_100_way_vertebrate = None  # what a failed sidecar left
        va.save()

        scored = {202: {"phylop_100_way_vertebrate": 9.873, "phastcons_100_way_vertebrate": 1.0,
                        "phylop_46_way_mammalian": 2.894, "phastcons_46_way_mammalian": 1.0}}
        with patch("annotation.management.commands.backfill_sv_conservation.score_sv_variants",
                   return_value=scored) as mock_score:
            call_command("backfill_sv_conservation", "--genome-build", "GRCh37")

        va.refresh_from_db()
        self.assertAlmostEqual(va.phylop_100_way_vertebrate, 9.873, places=3)
        self.assertAlmostEqual(va.phastcons_100_way_vertebrate, 1.0)

        vav = self.variant_annotation_versions_by_build['GRCh37']
        vav.refresh_from_db()
        self.assertTrue(vav.backfilled_sv_conservation)

        # A run VEP scored itself already holds VEP's values - the backfill doesn't re-score it
        annotation_run.sv_conservation_pybigwig = False
        annotation_run.save()
        with patch("annotation.management.commands.backfill_sv_conservation.score_sv_variants",
                   return_value=scored) as mock_score:
            call_command("backfill_sv_conservation", "--genome-build", "GRCh37")
        mock_score.assert_not_called()

    def test_recalculate_symbolic_hgvs(self):
        """ #1571: rows imported before symbolic DEL/DUP/INV resolved from coordinates hold a
            placeholder message - the backfill turns those into real HGVS """
        self._write_conservation_sidecar(self.TEST_ANNOTATION_VCF_GRCH37, self.CONSERVATION_GRCH37)
        self._import_grch37_run()

        va = VariantAnnotation.objects.get(variant_id=202)
        expected_hgvs_g = va.hgvs_g
        va.hgvs_g = VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE
        va.hgvs_c = VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE
        va.save()

        untouched = VariantAnnotation.objects.get(variant_id=203)
        untouched_hgvs_g = untouched.hgvs_g

        call_command("one_off_recalculate_symbolic_hgvs", "--genome-build", "GRCh37")

        va.refresh_from_db()
        self.assertEqual(expected_hgvs_g, va.hgvs_g)
        # No local transcript sequences in the fixture, so c.HGVS resolves to blank rather than the message
        self.assertIsNone(va.hgvs_c)

        untouched.refresh_from_db()
        self.assertEqual(untouched_hgvs_g, untouched.hgvs_g)
