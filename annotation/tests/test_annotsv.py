import os
import tempfile
from unittest import mock

from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation.annotation_run_files import get_annotsv_dir, write_qs_to_vcf
from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.annotsv_annotation import get_annotsv_command, get_annotsv_tsv_filename
from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import VariantAnnotation, VariantAnnotationPipelineType
from annotation.models.models import (
    AnnotationPipelineVersion,
    AnnotationRangeLock,
    AnnotationRun,
    VariantAnnotationVersion,
)
from annotation.pipelines.annotsv import AnnotSVRunner
from annotation.vcf_files.bulk_annotsv_tsv_inserter import (
    _extract_pathogenic_overlaps,
    _extract_variant_id,
    _row_to_update,
    import_annotsv_tsv,
)
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import (
    VEPConfig,
    get_vep_version_from_vcf,
    vep_dict_to_variant_annotation_version_kwargs,
)
from genes.models_enums import AnnotationConsortium
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_loci_and_variants_for_vcf,
    slowly_create_test_variant,
)

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


@override_settings(
    ANNOTATION_ANNOTSV_BIN="/fake/AnnotSV",
    ANNOTATION_ANNOTSV_ANNOTATIONS_DIR="/fake/anno",
    ANNOTATION_ANNOTSV_GENOME_BUILD={"GRCh37": "GRCh37"},
    ANNOTATION_ANNOTSV_EXTRA_ARGS=[],
    ANNOTATION_ANNOTSV_TIMEOUT_SECONDS=60,
    ANNOTATION_ANNOTSV_BUNDLE_VERSION="bundle-2024-01",
)
class TestRunAnnotsvSubprocessMocked(TestCase):
    """ What AnnotSVRunner.annotate() hands the tool. AnnotSV's header check rejects a sites-only VCF, so
        the dump is written with a dummy sample (dump_samples) and fed to AnnotSV as-is. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.genome_build)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.genome_build, AnnotationConsortium.REFSEQ, 4)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)

    def _make_annotsv_run(self) -> AnnotationRun:
        annotation_range_lock = AnnotationRangeLock.objects.create(
            version=self.vav, min_variant=self.variants[0], max_variant=self.variants[-1], count=1)
        pipeline_version = AnnotationPipelineVersion.objects.create(
            pipeline_type=VariantAnnotationPipelineType.ANNOTSV, genome_build=self.genome_build,
            code_version="3.5.8", data_version="bundle-2024-01",
            status=AnnotationPipelineVersion.Status.ACTIVE)
        return AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                            pipeline_type=VariantAnnotationPipelineType.ANNOTSV,
                                            pipeline_version=pipeline_version)

    def test_dump_carries_dummy_sample_and_is_what_runs(self):
        annotation_run = self._make_annotsv_run()
        runner = AnnotSVRunner()

        def _execute_cmd(cmd, **kwargs):  # pylint: disable=unused-argument
            output_dir = cmd[cmd.index("-outputDir") + 1]
            svinput = cmd[cmd.index("-SVinputFile") + 1]
            with open(get_annotsv_tsv_filename(svinput, output_dir), "w") as f:
                f.write("AnnotSV_ID\n")
            return 0, "ok", ""

        with tempfile.TemporaryDirectory() as dump_dir:
            with override_settings(ANNOTATION_VCF_DUMP_DIR=dump_dir), \
                    mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                               return_value="3.5.8"), \
                    mock.patch.object(AnnotSVRunner, "get_variants_qs",
                                      return_value=Variant.objects.filter(pk__in=[v.pk for v in self.variants])), \
                    mock.patch("annotation.pipelines.annotsv.execute_cmd",
                               side_effect=_execute_cmd) as execute_cmd_mock:
                runner.annotate(annotation_run)

                dump_filename = annotation_run.vcf_dump_filename
                self.assertEqual(annotation_run.vcf_annotated_filename,
                                 get_annotsv_tsv_filename(dump_filename, get_annotsv_dir(annotation_run)))

                # The dump itself is the AnnotSV input - no separate copy
                cmd_called = execute_cmd_mock.call_args.args[0]
                self.assertIn(dump_filename, cmd_called)

                with open(dump_filename) as f:
                    lines = f.read().splitlines()
                column_line = [line for line in lines if line.startswith("#CHROM")][0]
                self.assertTrue(column_line.endswith("\tFORMAT\tvariantgrid"))
                self.assertIn('##FORMAT=<ID=GT,', "\n".join(lines))
                records = [line for line in lines if not line.startswith("#")]
                self.assertTrue(records)
                for record in records:
                    self.assertTrue(record.endswith("\tGT\t0/1"), record)


class TestDumpSamples(TestCase):
    """ dump_samples plumbs through write_qs_to_vcf - a pipeline that doesn't set it still dumps
        sites-only, which is what VEP takes. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variant = slowly_create_test_variant("1", 100000, 'A', 'T', cls.genome_build)

    def test_no_samples_writes_sites_only(self):
        qs = Variant.objects.filter(pk=self.variant.pk)
        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf_filename = os.path.join(tmp_dir, "sites_only.vcf")
            write_qs_to_vcf(vcf_filename, self.genome_build, qs)
            with open(vcf_filename) as f:
                lines = f.read().splitlines()
        column_line = [line for line in lines if line.startswith("#CHROM")][0]
        self.assertNotIn("FORMAT", column_line)


@override_settings(
    ANNOTATION_ANNOTSV_BIN="/fake/AnnotSV",
    ANNOTATION_ANNOTSV_BUNDLE_VERSION="bundle-2024-01",
)
class TestCurrentAnnotsvVersion(TestCase):
    """ #720: registering the installed AnnotSV, which is how an upgrade reaches the scheduler. """

    def setUp(self):
        super().setUp()
        self.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.runner = AnnotSVRunner()

    def _register(self, code_version):
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        return_value=code_version):
            return self.runner.get_or_create_current_version(self.genome_build)

    def test_reuses_row_for_same_code_and_bundle(self):
        first, created = self._register("3.5.8")
        second, created_again = self._register("3.5.8")
        self.assertEqual(first.pk, second.pk)
        self.assertTrue(created)
        self.assertFalse(created_again)
        self.assertEqual(first.code_version, "3.5.8")
        self.assertEqual(first.data_version, "bundle-2024-01")

    def test_first_registered_version_is_active(self):
        """ Nothing to supersede, so no promote step - enabling AnnotSV is register-and-go. """
        first, _ = self._register("3.5.8")
        self.assertEqual(first.status, AnnotationPipelineVersion.Status.ACTIVE)

    def test_upgraded_binary_registers_as_new(self):
        """ An upgrade waits to be promoted, so swapping the binary never re-annotates by itself. """
        old, _ = self._register("3.5.8")
        new, _ = self._register("3.5.9")
        self.assertNotEqual(old.pk, new.pk)
        self.assertEqual(new.status, AnnotationPipelineVersion.Status.NEW)
        old.refresh_from_db()
        self.assertEqual(old.status, AnnotationPipelineVersion.Status.ACTIVE)
        self.assertEqual(AnnotationPipelineVersion.objects.filter(
            pipeline_type=VariantAnnotationPipelineType.ANNOTSV).count(), 2)

    def test_promoting_retires_the_version_it_replaces(self):
        old, _ = self._register("3.5.8")
        new, _ = self._register("3.5.9")

        new.promote_to_active()

        old.refresh_from_db()
        self.assertEqual(old.status, AnnotationPipelineVersion.Status.HISTORICAL)
        self.assertEqual(AnnotationPipelineVersion.get_active(
            VariantAnnotationPipelineType.ANNOTSV, self.genome_build), new)

    def test_check_tool_version_rejects_a_drifted_binary(self):
        """ A run carries the version it was scheduled against, so writing another tool's output into it
            would silently lose which AnnotSV produced the data. """
        pipeline_version, _ = self._register("3.5.8")
        annotation_run = mock.Mock(pipeline_version=pipeline_version, genome_build=self.genome_build)
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        return_value="3.5.9"):
            with self.assertRaises(ValueError):
                self.runner.check_tool_version(annotation_run)
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        return_value="3.5.8"):
            self.runner.check_tool_version(annotation_run)  # matching binary is fine


# The VEP run here only exists to write the rows AnnotSV updates - its fixture VCF has no #1657
# conservation sidecar, which an SV import otherwise requires.
@override_settings(**get_fake_annotation_settings_dict(columns_version=4),
                   ANNOTATION_VEP_SV_CONSERVATION_PYBIGWIG_ENABLED=False)
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
        cls.vep_run = AnnotationRun.objects.create(
            annotation_range_lock=annotation_range_lock,
            pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
            vcf_annotated_filename=TEST_SV_VCF_GRCH37,
        )
        import_vcf_annotations(cls.vep_run, delete_temp_files=False, vep_version_check=False)
        # #720: AnnotSV is its own run over the same range lock, updating the rows the VEP run wrote
        cls.annotation_run = AnnotationRun.objects.create(
            annotation_range_lock=annotation_range_lock,
            pipeline_type=VariantAnnotationPipelineType.ANNOTSV,
        )

    def test_import_full_lines_updates_variant_annotation(self):
        self.annotation_run.vcf_annotated_filename = TEST_ANNOTSV_TSV
        self.annotation_run.save()

        updated = import_annotsv_tsv(self.annotation_run)
        # The fixture references 3 variants (202, 203, 205) all of which have
        # VariantAnnotation rows from the VEP import above. The AnnotSV run owns none of them - it finds
        # them via the version partition and its range lock (#720).
        self.assertEqual(updated, 3)

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

    def test_has_annotsv_only_once_the_tsv_is_imported(self):
        """ Gates the AnnotSV section of the variant details page - SVs annotated before AnnotSV was
            enabled have VEP rows with none of its columns filled in. """
        va = VariantAnnotation.objects.get(variant_id=203)
        self.assertFalse(va.has_annotsv)

        self.annotation_run.vcf_annotated_filename = TEST_ANNOTSV_TSV
        self.annotation_run.save()
        import_annotsv_tsv(self.annotation_run)

        va.refresh_from_db()
        self.assertTrue(va.has_annotsv)

    def test_no_op_when_no_tsv_filename(self):
        self.annotation_run.vcf_annotated_filename = None
        self.annotation_run.save()
        self.assertEqual(import_annotsv_tsv(self.annotation_run), 0)

    def test_annotsv_run_owns_no_annotation_rows(self):
        """ #720: the annotation_run FK stays on the VEP run that created the rows, which is what makes
            resetting / deleting / re-running an AnnotSV run safe. """
        self.assertEqual(VariantAnnotation.objects.filter(annotation_run=self.annotation_run).count(), 0)
        self.assertTrue(VariantAnnotation.objects.filter(annotation_run=self.vep_run).exists())
        self.annotation_run.delete_related_objects()
        self.assertTrue(VariantAnnotation.objects.filter(annotation_run=self.vep_run).exists())
