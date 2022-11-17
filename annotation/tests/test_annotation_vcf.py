import copy
import os
from uuid import uuid4

from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation.annotation_versions import get_variant_annotation_version, \
    get_annotation_range_lock_and_unannotated_count
from annotation.models import VariantAnnotation
from annotation.models.damage_enums import PathogenicityImpact, ALoFTPrediction
from annotation.models.models import AnnotationRun, VariantAnnotationVersion, VariantTranscriptAnnotation
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import vep_parse_version_line, get_vep_version_from_vcf, \
    vep_dict_to_variant_annotation_version_kwargs, VEPVersionMismatchError, VEPConfig
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf

TEST_IMPORT_PROCESSING_DIR = os.path.join(settings.PRIVATE_DATA_ROOT, 'import_processing',
                                          "test", str(uuid4()))

TEST_ANNOTATION = copy.deepcopy(settings.ANNOTATION)
# I didn't have the phasCons/phyloP data on my laptop when I generated the test VCFs, so need to disable
TEST_ANNOTATION[settings.BUILD_GRCH37]["vep_config"].update({
    "phastcons100way": None,
    "phastcons46way": None,
    "phylop100way": None,
    "phylop46way": None,
})
TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"].update({
    "phastcons100way": None,
    "phastcons30way": None,
    "phylop100way": None,
    "phylop30way": None,
})

ANNOTATION_COLUMNS_V1 = copy.deepcopy(TEST_ANNOTATION)
ANNOTATION_COLUMNS_V1[settings.BUILD_GRCH37]["columns_version"] = 1
ANNOTATION_COLUMNS_V1[settings.BUILD_GRCH38]["columns_version"] = 1


@override_settings(IMPORT_PROCESSING_DIR=TEST_IMPORT_PROCESSING_DIR,
                   VARIANT_ZYGOSITY_GLOBAL_COLLECTION="global",
                   ANNOTATION_VEP_FAKE_VERSION=True,
                   ANNOTATION=ANNOTATION_COLUMNS_V1)
class TestAnnotationVCF(TestCase):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
    TEST_ANNOTATION_VCF_GRCH37 = os.path.join(TEST_DATA_DIR, "test_grch37.vep_annotated.vcf")
    TEST_ANNOTATION_VCF_GRCH38 = os.path.join(TEST_DATA_DIR, "test_grch38.vep_annotated.vcf")

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
                print("Truncating!")
                vav.truncate_related_objects()
            cls.variant_annotation_versions_by_build[genome_build_name] = vav

            slowly_create_loci_and_variants_for_vcf(genome_build, vcf, get_variant_id_from_info=True)

    def test_version_mismatch(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh37')
        vav = get_variant_annotation_version(genome_build)  # Will be fake due to ANNOTATION_VEP_FAKE_VERSION
        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH38)

        with self.assertRaises(VEPVersionMismatchError):
            import_vcf_annotations(annotation_run, delete_temp_files=False)

    def test_import_variant_annotations_grch37(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh37')
        vav = self.variant_annotation_versions_by_build[genome_build.name]

        print("Variants: ", Variant.objects.count())
        print("VariantAnnotation: ", VariantAnnotation.objects.count())
        print(f"VariantAnnotationVersion: {vav}")

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH37)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

        # Verify a few fields
        va = VariantAnnotation.objects.get(variant_id=13629762)
        self.assertEqual(va.impact, PathogenicityImpact.MODERATE)
        self.assertEqual(va.dbsnp_rs_id, "rs145886106")
        self.assertTrue("COSV" in va.cosmic_id)
        self.assertTrue("COSM" in va.cosmic_legacy_id)
        self.assertAlmostEqual(va.gnomad_af, 0.000145)

        va = VariantAnnotation.objects.get(variant_id=13638004)
        self.assertEqual(va.dbsnp_rs_id, "rs776172390")
        self.assertEqual(va.consequence, 'stop_gained')
        self.assertEqual(va.symbol, "ARHGAP11A")

        self._test_extra_grch37()

    def _test_extra_grch37(self):
        """ Columns_version specific tests """

        # This is testing columns_version 1
        va = VariantAnnotation.objects.get(variant_id=13629762)
        self.assertEqual(va.predictions_num_pathogenic, 2)
        self.assertEqual(va.predictions_num_benign, 3)

        va = VariantAnnotation.objects.get(variant_id=13638004)
        self.assertEqual(va.predictions_num_pathogenic, 1)
        self.assertEqual(va.predictions_num_benign, 0)

    def test_import_variant_annotations_grch38(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh38')
        vav = self.variant_annotation_versions_by_build[genome_build.name]

        print("Variants: ", Variant.objects.count())
        print("VariantAnnotation: ", VariantAnnotation.objects.count())
        print(f"VariantAnnotationVersion: {vav}")

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH38)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

        # Verify a few fields
        va = VariantAnnotation.objects.get(variant_id=24601)
        self.assertEqual(va.impact, PathogenicityImpact.MODERATE)
        self.assertEqual(va.dbsnp_rs_id, "rs145886106")
        self.assertEqual(va.cosmic_legacy_id, "COSM7286401")  # Test it has collapsed dupes
        self.assertAlmostEqual(va.gnomad_af, 0.000354913)
        self.assertEqual(va.gnomad_filtered, False)  # Test it converted FILTER properly to bool

        va = VariantAnnotation.objects.get(variant_id=42)
        self.assertEqual(va.dbsnp_rs_id, "rs776172390")
        self.assertEqual(va.consequence, 'stop_gained')
        self.assertEqual(va.symbol, "ARHGAP11A")

        # 2 | 20 | -5 | 20 | 0.05 | 0.17 | 0.01 | 0.04 | OR4F5
        va = VariantAnnotation.objects.get(variant_id=131165)
        self.assertEqual(va.spliceai_pred_dp_ag, 2)
        self.assertEqual(va.spliceai_pred_dp_al, 20)
        self.assertEqual(va.spliceai_pred_dp_dg, -5)
        self.assertEqual(va.spliceai_pred_dp_dl, 20)
        self.assertAlmostEqual(va.spliceai_pred_ds_ag, 0.05)
        self.assertAlmostEqual(va.spliceai_pred_ds_al, 0.17)
        self.assertAlmostEqual(va.spliceai_pred_ds_dg, 0.01)
        self.assertAlmostEqual(va.spliceai_pred_ds_dl, 0.04)
        self.assertEqual(va.spliceai_gene_symbol, "OR4F5")

        # RPE65:G197E | 1&1&1
        va = VariantAnnotation.objects.get(variant_id=131167)
        self._test_extra_grch38()

    def _test_extra_grch38(self):
        """ Do columns_version specific stuff """
        # This is testing columns_version 1
        va = VariantAnnotation.objects.get(variant_id=24601)
        self.assertEqual(va.predictions_num_pathogenic, 2)
        self.assertEqual(va.predictions_num_benign, 3)

        va = VariantAnnotation.objects.get(variant_id=42)
        self.assertEqual(va.predictions_num_pathogenic, 1)
        self.assertEqual(va.predictions_num_benign, 0)

        va = VariantAnnotation.objects.get(variant_id=131167)
        self.assertEqual(va.mastermind_count_1_cdna, 1)
        self.assertEqual(va.mastermind_count_2_cdna_prot, 1)
        self.assertEqual(va.mastermind_count_3_aa_change, 1)
        self.assertEqual(va.mastermind_mmid3, "RPE65:G197E")
        self.assertEqual(va.mutation_assessor_pred_most_damaging, 'M')
        self.assertEqual(va.polyphen2_hvar_pred_most_damaging, 'D')


ANNOTATION_COLUMNS_V2 = copy.deepcopy(TEST_ANNOTATION)
ANNOTATION_COLUMNS_V2[settings.BUILD_GRCH37]["columns_version"] = 2
ANNOTATION_COLUMNS_V2[settings.BUILD_GRCH38]["columns_version"] = 2


@override_settings(IMPORT_PROCESSING_DIR=TEST_IMPORT_PROCESSING_DIR,
                   VARIANT_ZYGOSITY_GLOBAL_COLLECTION="global",
                   ANNOTATION_VEP_FAKE_VERSION=True,
                   ANNOTATION=ANNOTATION_COLUMNS_V2)
class TestAnnotationVCF2(TestAnnotationVCF):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
    TEST_ANNOTATION_VCF_GRCH37 = os.path.join(TEST_DATA_DIR, "test_columns_version2_grch37.vep_annotated.vcf")
    TEST_ANNOTATION_VCF_GRCH38 = os.path.join(TEST_DATA_DIR, "test_columns_version2_grch38.vep_annotated.vcf")

    def _test_extra_grch37(self):
        # This is testing columns_version 2
        pass

    def _test_extra_grch38(self):
        # This is testing columns_version 2
        va = VariantAnnotation.objects.get(variant_id=24601)
        self.assertAlmostEqual(va.metalr_rankscore, 0.80456)
        self.assertAlmostEqual(va.revel_rankscore, 0.69527)
        self.assertAlmostEqual(va.vest4_rankscore, 0.56662)
        self.assertAlmostEqual(va.bayesdel_noaf_rankscore, 0.63287)
        self.assertAlmostEqual(va.cadd_raw_rankscore, 0.41304)
        self.assertAlmostEqual(va.clinpred_rankscore, 0.15198)

        va = VariantAnnotation.objects.get(variant_id=42)
        self.assertEqual(va.aloft_high_confidence, True)
        self.assertEqual(va.aloft_pred, ALoFTPrediction.RECESSIVE)
        self.assertAlmostEqual(va.aloft_prob_dominant, 0.13585)
        self.assertAlmostEqual(va.aloft_prob_recessive, 0.81255)
        self.assertAlmostEqual(va.aloft_prob_tolerant, 0.0516)

        vta = VariantTranscriptAnnotation.objects.get(variant_id=42, hgvs_c='NM_199357.3:c.1417C>T')
        self.assertTrue(vta.nmd_escaping_variant)


class TestVEP(TestCase):
    """ Random VEP annotation methods """
    maxDiff = None

    def test_parse_vep_version(self):
        LINES = [
            '##VEP="v97" time="2019-08-12 22:39:56" cache="/data/annotation/VEP/vep_cache/homo_sapiens/97_GRCh38" ensembl=97.378db18 ensembl-variation=97.26a059c ensembl-io=97.dc917e1 ensembl-funcgen=97.24f4d3c 1000genomes="phase3" COSMIC="88" ClinVar="201904" ESP="V2-SSA137" HGMD-PUBLIC="20184" assembly="GRCh38.p12" dbSNP="151" gencode="GENCODE 31" genebuild="2014-07" gnomAD="r2.1" polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"',
            '##VEP="v97" time="2019-08-08 13:51:20" cache="/media/dlawrence/SpinningIron/reference/VEP/vep_cache/homo_sapiens/97_GRCh37" ensembl=97 ensembl-funcgen=97 ensembl-variation=97 ensembl-io=97 1000genomes="phase3" COSMIC="86" ClinVar="201810" ESP="20141103" HGMD-PUBLIC="20174" assembly="GRCh37.p13" dbSNP="151" gencode="GENCODE 19" genebuild="2011-04" gnomAD="r2.1" polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"'
        ]
        for line in LINES:
            try:
                vep_parse_version_line(line)
            except:
                self.fail(f"vep_parse_version_line died on line: {line}")

    def test_aloft_pick_single(self):
        aloft_data = {
            "aloft_prob_tolerant": ".&0.0516&.&.&.&",
            "aloft_prob_recessive": ".&0.81255&.&.&.&",
            "aloft_prob_dominant": ".&0.13585&.&.&.&",
            "aloft_pred": ".&Recessive&.&.&.&",
            "aloft_high_confidence": ".&High&.&.&.&",
            "aloft_ensembl_transcript": "ENST00000565905&ENST00000361627&ENST00000567348&ENST00000563864&ENST00000543522",
        }
        BulkVEPVCFAnnotationInserter._pick_aloft_values(aloft_data)
        expected_aloft = {
            "aloft_prob_tolerant": "0.0516",
            "aloft_prob_recessive": "0.81255",
            "aloft_prob_dominant": "0.13585",
            "aloft_pred": "Recessive",
            "aloft_high_confidence": "High",
            "aloft_ensembl_transcript": "ENST00000361627",
        }
        self.assertDictEqual(aloft_data, expected_aloft)

    def test_aloft_pick_recessive_2_highs(self):
        aloft_data = {
            "aloft_prob_tolerant": "0.123&0.0516&.&.&.&",
            "aloft_prob_recessive": "0.123&0.81255&.&.&.&",
            "aloft_prob_dominant": "0.123&0.13585&.&.&.&",
            "aloft_pred": "Dominant&Recessive&.&.&.&",
            "aloft_high_confidence": "High&High&.&.&.&",
            "aloft_ensembl_transcript": "ENST00000565905&ENST00000361627&ENST00000567348&ENST00000563864&ENST00000543522",
        }
        BulkVEPVCFAnnotationInserter._pick_aloft_values(aloft_data)
        expected_aloft = {
            "aloft_prob_tolerant": "0.0516",
            "aloft_prob_recessive": "0.81255",
            "aloft_prob_dominant": "0.13585",
            "aloft_pred": "Recessive",
            "aloft_high_confidence": "High",
            "aloft_ensembl_transcript": "ENST00000361627",
        }
        self.assertDictEqual(aloft_data, expected_aloft)

    def test_aloft_pick_high_dominant(self):
        aloft_data = {
            "aloft_prob_tolerant": "0.123&0.0516&.&.&.&",
            "aloft_prob_recessive": "0.123&0.81255&.&.&.&",
            "aloft_prob_dominant": "0.123&0.13585&.&.&.&",
            "aloft_pred": "Dominant&Recessive&.&.&.&",
            "aloft_high_confidence": "High&Low&.&.&.&",
            "aloft_ensembl_transcript": "ENST00000565905&ENST00000361627&ENST00000567348&ENST00000563864&ENST00000543522",
        }
        BulkVEPVCFAnnotationInserter._pick_aloft_values(aloft_data)
        expected_aloft = {
            "aloft_prob_tolerant": "0.123",
            "aloft_prob_recessive": "0.123",
            "aloft_prob_dominant": "0.123",
            "aloft_pred": "Dominant",
            "aloft_high_confidence": "High",
            "aloft_ensembl_transcript": "ENST00000565905",
        }
        self.assertDictEqual(aloft_data, expected_aloft)
