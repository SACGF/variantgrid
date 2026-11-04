from django.test import TestCase

from annotation.enums import (
    AlphaMissensePrediction,
    MetaRNNPrediction,
    PathogenicityImpact,
    SIFTPrediction,
    VariantAnnotationPipelineType,
    VEPCustom,
)
from annotation.vep_columns import (
    VEP_COLUMNS,
    all_variant_grid_column_ids,
    filter_for,
)
from snpdb.models import VariantGridColumn


def _build_field_formatters(**filter_kwargs) -> dict:
    """ Mirror BulkVEPVCFAnnotationInserter._add_vep_field_handlers: destination column
        -> formatter, built straight from the active VEPColumnDefs. """
    return {
        vgc: cvf.formatter
        for cvf in filter_for(**filter_kwargs)
        for vgc in cvf.variant_grid_columns
        if cvf.formatter is not None
    }


class VepColumnsRegistryTest(TestCase):

    def test_referenced_variant_grid_columns_exist(self):
        known = set(VariantGridColumn.objects.values_list("pk", flat=True))
        referenced = all_variant_grid_column_ids()
        missing = referenced - known
        self.assertFalse(missing, f"vep_columns.py references unknown VariantGridColumn ids: {sorted(missing)}")

    def test_dedup_keys_unique(self):
        seen = set()
        dupes = []
        for c in VEP_COLUMNS:
            k = (c.source_field, c.vep_plugin, c.vep_custom,
                 c.min_columns_version, c.max_columns_version,
                 c.min_vep_version, c.max_vep_version,
                 c.min_cosmic_version, c.max_cosmic_version,
                 c.genome_builds, c.pipeline_types,
                 c.variant_grid_columns)
            if k in seen:
                dupes.append(k)
            seen.add(k)
        self.assertFalse(dupes, f"Duplicate VEPColumnDef keys: {dupes}")

    def test_gnomad4_minor_version_filter(self):
        common_kwargs = dict(
            genome_build_name="GRCh38",
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=3,
            vep_custom=VEPCustom.GNOMAD_4,
        )
        rows_40 = [c for c in filter_for(gnomad4_minor_version="4.0", **common_kwargs)
                   if "gnomad_filtered" in c.variant_grid_columns]
        rows_41 = [c for c in filter_for(gnomad4_minor_version="4.1", **common_kwargs)
                   if "gnomad_filtered" in c.variant_grid_columns]
        self.assertEqual(len(rows_40), 1)
        self.assertEqual(rows_40[0].source_field, "gnomad_filtered")
        self.assertEqual(len(rows_41), 1)
        self.assertEqual(rows_41[0].source_field, "FILTER")

    def test_cosmic_count_source_field_per_release(self):
        """ #1673: COSMIC renames the sample count INFO field, so the release picks the field - and
            with it the CSQ column the inserter reads. """
        expected_source_fields = {
            95: "CNT",
            97: "CNT",
            99: "SAMPLE_COUNT",
            101: "GENOME_SCREEN_SAMPLE_COUNT",
        }
        for cosmic_version, source_field in expected_source_fields.items():
            with self.subTest(cosmic_version=cosmic_version):
                rows = [c for c in filter_for(genome_build_name="GRCh37",
                                              pipeline_type=VariantAnnotationPipelineType.STANDARD,
                                              columns_version=5,
                                              vep_custom=VEPCustom.COSMIC,
                                              cosmic_version=cosmic_version)
                        if "cosmic_count" in c.variant_grid_columns]
                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0].source_field, source_field)
                self.assertEqual(rows[0].vep_info_field, f"COSMIC_{source_field}")

    @staticmethod
    def _gnomad4_columns(genome_build_name: str, gnomad4_minor_version: str = "4.1") -> list[str]:
        rows = filter_for(
            genome_build_name=genome_build_name,
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=3,
            vep_custom=VEPCustom.GNOMAD_4,
            gnomad4_minor_version=gnomad4_minor_version,
        )
        return [vgc for c in rows for vgc in c.variant_grid_columns]

    def test_gnomad4_faf_is_grch38_only(self):
        """ The faf95/faf99 rows narrow to GRCh38 only; T2T keeps the rest of gnomAD 4, fafmax included """
        grch38 = self._gnomad4_columns("GRCh38")
        t2t = self._gnomad4_columns("T2T-CHM13v2.0")
        for vgc in ("gnomad_faf95", "gnomad_faf99"):
            self.assertIn(vgc, grch38)
            self.assertNotIn(vgc, t2t)
        self.assertIn("gnomad_fafmax_faf95_max", t2t)

    def test_filter_for_gnomad4_fills_each_column_once(self):
        """ The minor version picks one source per column (gnomad_filtered has a 4.0 and a 4.1 def) """
        for genome_build_name in ("GRCh38", "T2T-CHM13v2.0"):
            for gnomad4_minor_version in ("4.0", "4.1"):
                with self.subTest(genome_build_name=genome_build_name, gnomad4_minor_version=gnomad4_minor_version):
                    columns = self._gnomad4_columns(genome_build_name, gnomad4_minor_version)
                    self.assertIn("gnomad_filtered", columns)
                    self.assertEqual(len(columns), len(set(columns)))

    def test_filter_for_gnomad3_excluded_at_v3(self):
        rows = filter_for(
            genome_build_name="GRCh38",
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=3,
            vep_custom=VEPCustom.GNOMAD_3,
        )
        self.assertEqual(rows, ())

    # ----- field_formatters snapshot --------------------------------------------
    # Pin the destination-column -> formatter map the config builds, so a config edit
    # that silently drops or re-points a formatter shows up here.

    def test_field_formatters_snapshot_grch38_v4(self):
        """ Representative modern config (GRCh38, columns_version 4, gnomAD 4.1, VEP 112). """
        formatters = _build_field_formatters(
            genome_build_name="GRCh38",
            pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=4,
            vep_version=112,
            gnomad4_minor_version="4.1",
        )
        expected_keys = {
            "af_1kg", "af_uk10k",
            "aloft_ensembl_transcript", "aloft_high_confidence", "aloft_pred",
            "aloft_prob_dominant", "aloft_prob_recessive", "aloft_prob_tolerant",
            "alphamissense_pred", "alphamissense_score",
            "bayesdel_noaf_score", "cadd_phred", "cadd_raw", "canonical", "clinpred_score",
            "cosmic_count", "cosmic_id", "cosmic_legacy_id", "dbsnp_rs_id",
            "denovo_db_case_count", "denovo_db_control_count",
            "gnomad2_liftover_af", "gnomad_filtered", "gnomad_non_par", "gnomad_popmax",
            "hgnc_id", "impact", "interpro_domain",
            "mastermind_count_1_cdna", "mastermind_count_2_cdna_prot", "mastermind_count_3_aa_change",
            "mavedb_score", "metarnn_pred", "metarnn_score", "mpc_score",
            "mutpred2_score", "mutpred2_top5_mechanisms", "nmd_escaping_variant",
            "phastcons_100_way_vertebrate", "phastcons_30_way_mammalian",
            "phylop_100_way_vertebrate", "phylop_30_way_mammalian",
            "primateai_pred", "primateai_score", "revel_score",
            "sift", "somatic", "topmed_af", "variant_class",
            "varity_er_score", "varity_r_score", "vest4_score",
        }
        self.assertEqual(set(formatters), expected_keys)

        # Behaviour spot-checks across every formatter family, matching the old inserter dict.
        self.assertEqual(formatters["af_1kg"]("0.6764&0.2433"), 0.6764)           # pick highest float
        self.assertEqual(formatters["cosmic_count"]("3&10"), 10)                  # pick highest int
        self.assertEqual(formatters["denovo_db_case_count"]("1&2&3"), 6)          # sum int
        self.assertEqual(formatters["mavedb_score"]("NA&-1.5&2.0"), -1.5)         # pick lowest float, NA-aware
        self.assertEqual(formatters["gnomad_popmax"]("nfe"), "NFE")              # str.upper
        self.assertFalse(formatters["gnomad_filtered"]("PASS"))                   # FILTER -> bool
        self.assertTrue(formatters["gnomad_filtered"]("RF"))
        self.assertTrue(formatters["canonical"]("YES"))
        self.assertFalse(formatters["canonical"](""))
        self.assertIsNone(formatters["gnomad_non_par"]("."))                      # empty -> None
        self.assertEqual(formatters["gnomad_non_par"]("Y"), "Y")
        self.assertEqual(formatters["hgnc_id"]("HGNC:55"), 55)
        self.assertEqual(formatters["cosmic_id"]("COSV123&rs456"), "COSV123")
        self.assertEqual(formatters["dbsnp_rs_id"]("COSV123&rs456"), "rs456")
        self.assertEqual(formatters["sift"]("deleterious"), SIFTPrediction.DAMAGING)
        self.assertTrue(formatters["somatic"]("0&1"))
        self.assertEqual(formatters["mastermind_count_2_cdna_prot"]("3&5&7"), 5)
        self.assertEqual(formatters["impact"]("HIGH"), PathogenicityImpact.HIGH)
        self.assertEqual(formatters["metarnn_pred"]("T&D"),
                         MetaRNNPrediction.get_most_damaging(["T", "D"]))
        self.assertEqual(formatters["alphamissense_pred"]("B&LP"),
                         AlphaMissensePrediction.LIKELY_PATHOGENIC)
        self.assertTrue(formatters["nmd_escaping_variant"]("NMD_escaping_variant"))
        self.assertEqual(formatters["interpro_domain"]("just_one"), "just_one")   # remove empty multiples

        # Rankscores and clinpred_pred carried no formatter in the old dict - still none.
        self.assertNotIn("revel_rankscore", formatters)
        self.assertNotIn("clinpred_pred", formatters)
        self.assertNotIn("gerp_pp_rs", formatters)

    def test_field_formatters_gnomad_filtered_gated_by_source(self):
        """ gnomad_filtered gets gnomad_filtered_func only where it is FILTER-sourced:
            gnomAD3 (GRCh38 <=v2) and gnomAD4 4.1. The gnomad_filtered=1 defs (gnomAD4 4.0,
            gnomAD2 GRCh37) auto-convert the bool, so carry no formatter. """
        # gnomAD3 (GRCh38, columns_version 2) - FILTER text needs parsing
        f_v2 = _build_field_formatters(
            genome_build_name="GRCh38", pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=2,
        )
        self.assertIn("gnomad_filtered", f_v2)
        self.assertFalse(f_v2["gnomad_filtered"]("PASS"))
        self.assertTrue(f_v2["gnomad_filtered"]("RF"))

        # gnomAD4 4.0 - gnomad_filtered=1 auto-converts, no formatter
        f_40 = _build_field_formatters(
            genome_build_name="GRCh38", pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=4, gnomad4_minor_version="4.0",
        )
        self.assertNotIn("gnomad_filtered", f_40)

        # gnomAD2 (GRCh37) - gnomad_filtered=1 auto-converts, no formatter
        f_37 = _build_field_formatters(
            genome_build_name="GRCh37", pipeline_type=VariantAnnotationPipelineType.STANDARD,
            columns_version=4,
        )
        self.assertNotIn("gnomad_filtered", f_37)

    def test_field_formatters_no_vgc_maps_to_two_formatters(self):
        """ Across every active config, a destination column never resolves to two different
            formatters (the build-mutually-exclusive gnomAD families guarantee this). """
        builds = ["GRCh37", "GRCh38", "T2T-CHM13v2.0"]
        for build in builds:
            for cv in range(1, 6):
                for vep_version in (111, 112, 116):
                    for gnomad4_minor in ("4.0", "4.1", None):
                        for pt in (VariantAnnotationPipelineType.STANDARD,
                                   VariantAnnotationPipelineType.STRUCTURAL_VARIANT):
                            seen = {}
                            for cvf in filter_for(genome_build_name=build, pipeline_type=pt,
                                                  columns_version=cv, vep_version=vep_version,
                                                  gnomad4_minor_version=gnomad4_minor):
                                if cvf.formatter is None:
                                    continue
                                for vgc in cvf.variant_grid_columns:
                                    prev = seen.get(vgc)
                                    self.assertTrue(prev is None or prev is cvf.formatter,
                                                    f"{vgc} maps to two formatters for "
                                                    f"{build}/cv{cv}/vep{vep_version}/{gnomad4_minor}/{pt.name}")
                                    seen[vgc] = cvf.formatter
