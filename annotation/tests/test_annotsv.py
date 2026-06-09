import os
from unittest import mock

from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.annotsv_annotation import (
    annotsv_check_command_line_version_match,
    AnnotSVVersionMismatchError,
    get_annotsv_command,
    run_annotsv,
)
from annotation.fake_annotation import get_fake_annotation_settings_dict
from annotation.models import VariantAnnotation, VariantAnnotationPipelineType
from annotation.models.models import AnnotationRun, VariantAnnotationVersion
from annotation.vcf_files.bulk_annotsv_tsv_inserter import (
    _extract_pathogenic_overlaps,
    _extract_variant_id,
    _row_to_update,
    import_annotsv_tsv,
)
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import (
    get_vep_version_from_vcf,
    vep_dict_to_variant_annotation_version_kwargs,
    VEPConfig,
)
from genes.models_enums import AnnotationConsortium
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf


TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
TEST_ANNOTSV_TSV = os.path.join(TEST_DATA_DIR, "annotsv", "test_grch37_sv.annotated.tsv")
TEST_SV_VCF_GRCH37 = os.path.join(TEST_DATA_DIR, "test_columns_version4_grch37_sv.vep_annotated.vcf")


class TestExtractVariantId(TestCase):
    def test_from_info(self):
        self.assertEqual(_extract_variant_id({"INFO": "SVTYPE=DEL;variant_id=42"}), 42)
        self.assertEqual(_extract_variant_id({"INFO": "variant_id=7;END=100"}), 7)

    def test_from_id_when_info_missing(self):
        self.assertEqual(_extract_variant_id({"INFO": ".", "ID": "13"}), 13)

    def test_returns_none_when_unparseable(self):
        self.assertIsNone(_extract_variant_id({"INFO": ".", "ID": "."}))
        self.assertIsNone(_extract_variant_id({}))


class TestRowToUpdate(TestCase):
    def test_parses_int_and_float(self):
        row = {
            "ACMG_class": "5",
            "AnnotSV_ranking_score": "2.5",
            "RE_gene": "enh_X",
            "Repeat_type_left": "NA",
            "B_gain_AFmax": "0.005",
        }
        update = _row_to_update(row)
        self.assertEqual(update["annotsv_acmg_class"], 5)
        self.assertEqual(update["annotsv_acmg_score"], 2.5)
        self.assertEqual(update["annotsv_re_gene"], "enh_X")
        self.assertNotIn("annotsv_repeat_type_left", update)  # NA is skipped
        self.assertEqual(update["annotsv_b_gain_af_max"], 0.005)

    def test_skips_empty_and_dot(self):
        row = {
            "ACMG_class": "",
            "RE_gene": ".",
            "AnnotSV_ranking_score": "NA",
        }
        self.assertEqual(_row_to_update(row), {})

    def test_parses_new_typed_fields(self):
        row = {
            "AnnotSV_ranking_criteria": "1A,2C",
            "Frameshift": "yes",
            "Exons_spanned": "3",
            "Dist_nearest_SS": "120",
            "Nearest_SS_type": "5'",
            "OMIM_inheritance": "AD",
            "OMIM_morbid": "yes",
            "OMIM_phenotype": "Breast cancer",
            "OMIM_ID": "113705",
        }
        update = _row_to_update(row)
        self.assertEqual(update["annotsv_acmg_criteria"], "1A,2C")
        self.assertIs(update["annotsv_frameshift"], True)
        self.assertEqual(update["annotsv_exons_spanned"], 3)
        self.assertEqual(update["annotsv_dist_nearest_ss"], 120)
        self.assertEqual(update["annotsv_nearest_ss_type"], "5'")
        self.assertEqual(update["annotsv_omim_inheritance"], "AD")
        self.assertIs(update["annotsv_omim_morbid"], True)
        self.assertEqual(update["annotsv_omim_phenotype"], "Breast cancer")
        self.assertEqual(update["annotsv_omim_id"], "113705")

    def test_bool_field_na(self):
        # NA -> not present in update
        update = _row_to_update({"Frameshift": "NA", "OMIM_morbid": "NA"})
        self.assertNotIn("annotsv_frameshift", update)
        self.assertNotIn("annotsv_omim_morbid", update)
        # yes -> True
        self.assertIs(_row_to_update({"OMIM_morbid": "yes"})["annotsv_omim_morbid"], True)
        # no -> False
        self.assertIs(_row_to_update({"OMIM_morbid": "no"})["annotsv_omim_morbid"], False)
        # bogus -> dropped (parser returns None)
        self.assertNotIn("annotsv_omim_morbid",
                         _row_to_update({"OMIM_morbid": "maybe"}))


class TestPathogenicOverlapsAssembly(TestCase):
    def test_drops_empty_keeps_populated(self):
        row = {
            "P_loss_source": "ClinGen",
            "P_loss_phen": "Hereditary breast cancer",
            "P_loss_hpo": "NA",
            "P_loss_coord": "chr17:41200000-41260000",
            "P_inv_source": "ClinVar&dbVar",
            "P_inv_phen": "NA",
            "P_inv_hpo": "HP:0001234",
            "P_inv_coord": "chr3:127500000-128800000",
            "P_gain_source": "NA",
            "P_gain_phen": "NA",
            "P_gain_hpo": "NA",
            "P_gain_coord": "NA",
            "P_ins_source": "",
            "P_ins_phen": ".",
            "P_ins_hpo": "NA",
            "P_ins_coord": "NA",
        }
        overlaps = _extract_pathogenic_overlaps(row)
        self.assertEqual(set(overlaps.keys()), {"loss", "inv"})
        self.assertEqual(overlaps["loss"], {
            "source": "ClinGen",
            "phen": "Hereditary breast cancer",
            "coord": "chr17:41200000-41260000",
        })
        self.assertEqual(overlaps["inv"], {
            "source": "ClinVar&dbVar",
            "hpo": "HP:0001234",
            "coord": "chr3:127500000-128800000",
        })

    def test_returns_none_when_all_na(self):
        row = {
            f"P_{event}_{sub}": "NA"
            for event in ("gain", "loss", "ins", "inv")
            for sub in ("source", "phen", "hpo", "coord")
        }
        self.assertIsNone(_extract_pathogenic_overlaps(row))

    def test_row_to_update_includes_overlaps(self):
        row = {
            "P_loss_source": "ClinGen",
            "P_loss_coord": "chr1:1-2",
        }
        update = _row_to_update(row)
        self.assertEqual(update["annotsv_pathogenic_overlaps"],
                         {"loss": {"source": "ClinGen", "coord": "chr1:1-2"}})

    def test_row_to_update_omits_overlaps_when_none(self):
        update = _row_to_update({"ACMG_class": "5"})
        self.assertNotIn("annotsv_pathogenic_overlaps", update)


class TestGetAnnotsvCommand(TestCase):
    @override_settings(
        ANNOTATION_ANNOTSV_BIN="/fake/AnnotSV",
        ANNOTATION_ANNOTSV_ANNOTATIONS_DIR="/fake/anno",
        ANNOTATION_ANNOTSV_GENOME_BUILD={"GRCh37": "GRCh37", "GRCh38": "GRCh38"},
        ANNOTATION_ANNOTSV_EXTRA_ARGS=["-SVminSize", "50"],
    )
    def test_command_includes_required_flags(self):
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cmd = get_annotsv_command("/tmp/in.vcf", "/tmp/out", genome_build,
                                  AnnotationConsortium.REFSEQ)
        self.assertIn("/fake/AnnotSV", cmd)
        self.assertIn("-SVinputFile", cmd)
        self.assertIn("/tmp/in.vcf", cmd)
        self.assertIn("-genomeBuild", cmd)
        self.assertIn("-tx", cmd)
        # AnnotationConsortium.REFSEQ should map to "RefSeq"
        tx_idx = cmd.index("-tx")
        self.assertEqual(cmd[tx_idx + 1], "RefSeq")
        # extra args appended
        self.assertEqual(cmd[-2:], ["-SVminSize", "50"])

    @override_settings(
        ANNOTATION_ANNOTSV_BIN="/fake/AnnotSV",
        ANNOTATION_ANNOTSV_ANNOTATIONS_DIR="/fake/anno",
        ANNOTATION_ANNOTSV_GENOME_BUILD={"GRCh37": "GRCh37", "GRCh38": "GRCh38"},
        ANNOTATION_ANNOTSV_EXTRA_ARGS=[],
    )
    def test_ensembl_consortium_passes_ENSEMBL(self):
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cmd = get_annotsv_command("/tmp/in.vcf", "/tmp/out", genome_build,
                                  AnnotationConsortium.ENSEMBL)
        tx_idx = cmd.index("-tx")
        self.assertEqual(cmd[tx_idx + 1], "ENSEMBL")


class TestRunAnnotsvSubprocessMocked(TestCase):
    @override_settings(
        ANNOTATION_ANNOTSV_BIN="/fake/AnnotSV",
        ANNOTATION_ANNOTSV_ANNOTATIONS_DIR="/fake/anno",
        ANNOTATION_ANNOTSV_GENOME_BUILD={"GRCh37": "GRCh37"},
        ANNOTATION_ANNOTSV_EXTRA_ARGS=[],
        ANNOTATION_ANNOTSV_TIMEOUT_SECONDS=60,
    )
    def test_subprocess_invocation_and_filename(self):
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        with mock.patch("annotation.annotsv_annotation.subprocess.run") as run_mock, \
                mock.patch("annotation.annotsv_annotation.os.makedirs"), \
                mock.patch("annotation.annotsv_annotation.os.path.exists", return_value=True):
            run_mock.return_value = mock.Mock(returncode=0, stdout="ok", stderr="")
            tsv, rc, _stdout, _stderr = run_annotsv(
                "/tmp/dump_99_structural_variant.vcf", "/tmp/out_dir",
                genome_build, AnnotationConsortium.REFSEQ,
            )
        self.assertEqual(rc, 0)
        self.assertEqual(tsv, "/tmp/out_dir/dump_99_structural_variant.annotated.tsv")
        run_mock.assert_called_once()
        cmd_called = run_mock.call_args.args[0]
        self.assertIn("/tmp/dump_99_structural_variant.vcf", cmd_called)


class TestVersionCheck(TestCase):
    @override_settings(ANNOTATION_ANNOTSV_ENABLED=False)
    def test_skipped_when_disabled(self):
        # Should be a no-op even with a totally bogus VAV stand-in
        annotsv_check_command_line_version_match(mock.Mock(annotsv_code="1.0", annotsv_bundle="X"))

    @override_settings(
        ANNOTATION_ANNOTSV_ENABLED=True,
        ANNOTATION_ANNOTSV_BUNDLE_VERSION="bundle-2024-01",
    )
    def test_raises_on_code_mismatch(self):
        vav = mock.Mock(annotsv_code="3.5.8", annotsv_bundle=None)
        with mock.patch("annotation.annotsv_annotation.get_annotsv_command_line_version",
                        return_value="3.5.7"):
            with self.assertRaises(AnnotSVVersionMismatchError):
                annotsv_check_command_line_version_match(vav)

    @override_settings(
        ANNOTATION_ANNOTSV_ENABLED=True,
        ANNOTATION_ANNOTSV_BUNDLE_VERSION="bundle-old",
    )
    def test_raises_on_bundle_mismatch(self):
        vav = mock.Mock(annotsv_code=None, annotsv_bundle="bundle-new")
        annotsv_check_command_line_version_match  # alias to avoid lint
        with self.assertRaises(AnnotSVVersionMismatchError):
            annotsv_check_command_line_version_match(vav)


@override_settings(**get_fake_annotation_settings_dict(columns_version=4))
class TestImportAnnotsvTsv(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        vep_config = VEPConfig(genome_build)
        vep_dict = get_vep_version_from_vcf(TEST_SV_VCF_GRCH37)
        kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_dict)
        kwargs["genome_build"] = genome_build
        vav, created = VariantAnnotationVersion.objects.get_or_create(**kwargs)
        if not created:
            vav.truncate_related_objects()
        cls.vav = vav
        cls.genome_build = genome_build

        slowly_create_loci_and_variants_for_vcf(genome_build, TEST_SV_VCF_GRCH37,
                                                get_variant_id_from_info=True)

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        cls.annotation_run = AnnotationRun.objects.create(
            annotation_range_lock=annotation_range_lock,
            pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
            vcf_annotated_filename=TEST_SV_VCF_GRCH37,
        )
        import_vcf_annotations(cls.annotation_run, delete_temp_files=False, vep_version_check=False)

    def test_import_full_lines_updates_variant_annotation(self):
        self.annotation_run.annotsv_tsv_filename = TEST_ANNOTSV_TSV
        self.annotation_run.save()

        updated = import_annotsv_tsv(self.annotation_run)
        # The fixture references 3 variants (202, 203, 205) all of which have
        # VariantAnnotation rows from the VEP import above.
        self.assertEqual(updated, 3)

        self.annotation_run.refresh_from_db()
        self.assertTrue(self.annotation_run.annotsv_imported)

        va_202 = VariantAnnotation.objects.get(variant_id=202)
        self.assertEqual(va_202.annotsv_acmg_class, 2)
        self.assertAlmostEqual(va_202.annotsv_acmg_score, 0.123)
        self.assertEqual(va_202.annotsv_pathogenic_overlaps, {
            "inv": {
                "phen": "Inversion test phenotype",
                "hpo": "HP:0001234",
                "source": "ClinVar&dbVar",
                "coord": "chr3:127500000-128800000",
            },
        })

        va_203 = VariantAnnotation.objects.get(variant_id=203)
        self.assertEqual(va_203.annotsv_acmg_class, 5)
        self.assertAlmostEqual(va_203.annotsv_acmg_score, 2.5)
        self.assertEqual(va_203.annotsv_acmg_criteria, "1A,2C")
        self.assertEqual(va_203.annotsv_segdup_left, "segdup_X")
        self.assertAlmostEqual(va_203.annotsv_b_ins_af_max, 0.005)
        self.assertIs(va_203.annotsv_frameshift, True)
        self.assertEqual(va_203.annotsv_exons_spanned, 3)
        self.assertEqual(va_203.annotsv_dist_nearest_ss, 120)
        self.assertEqual(va_203.annotsv_nearest_ss_type, "5'")
        self.assertEqual(va_203.annotsv_omim_inheritance, "AD")
        self.assertIs(va_203.annotsv_omim_morbid, True)
        self.assertEqual(va_203.annotsv_omim_phenotype, "Breast cancer")
        self.assertEqual(va_203.annotsv_omim_id, "113705")
        self.assertEqual(va_203.annotsv_pathogenic_overlaps, {
            "loss": {
                "phen": "Hereditary breast cancer",
                "hpo": "HP:0003002",
                "source": "ClinGen",
                "coord": "chr17:41200000-41260000",
            },
        })

        va_205 = VariantAnnotation.objects.get(variant_id=205)
        self.assertEqual(va_205.annotsv_acmg_class, 1)
        self.assertEqual(va_205.annotsv_re_gene, "enh_RUNX1")
        self.assertEqual(va_205.annotsv_encode_blacklist_left, "blocklist")
        self.assertEqual(va_205.annotsv_encode_blacklist_characteristics_left, "low_mappability")
        self.assertIs(va_205.annotsv_omim_morbid, False)
        self.assertIsNone(va_205.annotsv_pathogenic_overlaps)

    def test_no_op_when_no_tsv_filename(self):
        # Ensure unset filename is a no-op (and does not flip annotsv_imported)
        self.annotation_run.annotsv_tsv_filename = None
        self.annotation_run.annotsv_imported = False
        self.annotation_run.save()
        self.assertEqual(import_annotsv_tsv(self.annotation_run), 0)
        self.annotation_run.refresh_from_db()
        self.assertFalse(self.annotation_run.annotsv_imported)
