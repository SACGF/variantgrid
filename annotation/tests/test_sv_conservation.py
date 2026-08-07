"""
Tests for the pyBigWig SV conservation scorer (#1657).

Scores the structural variants in an existing VEP-annotated fixture VCF with pyBigWig and asserts the 4
conservation _max columns (phastCons/phyloP) agree with the values VEP wrote into the fixture's CSQ,
within tolerance. Requires the real conservation bigWig data files - skipped if they aren't installed.
"""
import os
import re
import unittest

from django.test import TestCase

from annotation.sv_conservation import (
    get_sv_conservation_tracks,
    score_sv_vcf,
    sv_conservation_window,
)
from snpdb.models import GenomeBuild

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
TOLERANCE = 0.05

# Build -> VEP-annotated SV fixture whose CSQ carries the phastCons*/phyloP* _max ground-truth values.
BUILD_FIXTURES = {
    "GRCh37": "test_columns_version4_grch37_sv.vep_annotated.vcf",
    "GRCh38": "test_columns_version4_grch38_sv.vep_annotated.vcf",
}


def _parse_fixture(filename):
    """ -> ({variant_id: {vep_field: float}}, csq_format_fields). Reads the representative (first)
        CSQ record per variant and pulls out every *_max conservation field. """
    csq_format = None
    expected = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_format = re.search(r'Format: ([^"]+)', line).group(1).split("|")
            if line.startswith("#"):
                continue
            info = line.rstrip("\n").split("\t")[7]
            kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            variant_id = int(kv["variant_id"])
            csq = dict(zip(csq_format, kv["CSQ"].split(",")[0].split("|")))
            maxes = {k: float(v) for k, v in csq.items()
                     if k.endswith("_max") and v not in (None, "")}
            expected[variant_id] = maxes
    return expected, csq_format


class SVConservationWindowTests(TestCase):
    def test_window_del_dup_uses_svlen(self):
        # POS-1 (0-based inclusive) .. POS+|SVLEN| (0-based exclusive)
        self.assertEqual(sv_conservation_window(pos=1000, end=None, svlen=500), (999, 1500))

    def test_window_negative_svlen_abs(self):
        self.assertEqual(sv_conservation_window(pos=1000, end=None, svlen=-500), (999, 1500))

    def test_window_ins_uses_svlen_over_small_end(self):
        # INS: tiny END (POS+1) but large SVLEN -> footprint is POS+|SVLEN|
        self.assertEqual(sv_conservation_window(pos=1000, end=1001, svlen=800), (999, 1800))


class SVConservationScoringTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.builds = {}
        for build_name in BUILD_FIXTURES:
            try:
                cls.builds[build_name] = GenomeBuild.get_name_or_alias(build_name)
            except GenomeBuild.DoesNotExist:
                pass

    def _check_build(self, build_name):
        genome_build = self.builds.get(build_name)
        if genome_build is None:
            self.skipTest(f"{build_name} genome build not available")

        tracks = get_sv_conservation_tracks(genome_build)
        if not tracks:
            self.skipTest(f"No conservation bigWig tracks configured for {build_name}")
        missing = [t.path for t in tracks if not os.path.exists(t.path)]
        if missing:
            self.skipTest(f"Conservation bigWig data files not installed: {missing}")

        fixture = os.path.join(TEST_DATA_DIR, BUILD_FIXTURES[build_name])
        expected, _ = _parse_fixture(fixture)

        results = score_sv_vcf(fixture, genome_build)
        self.assertTrue(results, "Scorer returned no results")

        # track.name is the VEP custom short_name; its CSQ field is "<name>_max".
        checks = 0
        for variant_id, expected_maxes in expected.items():
            variant_results = results.get(variant_id)
            self.assertIsNotNone(variant_results, f"No pyBigWig result for variant {variant_id}")
            for track in tracks:
                vep_field = f"{track.name}_max"
                if vep_field not in expected_maxes:
                    continue
                got = variant_results.get(track.db_column)
                self.assertIsNotNone(
                    got, f"variant {variant_id} track {track.name} ({track.db_column}) not scored")
                self.assertLessEqual(
                    abs(got - expected_maxes[vep_field]), TOLERANCE,
                    f"variant {variant_id} {track.name}: pyBigWig={got} VEP={expected_maxes[vep_field]}")
                checks += 1
        self.assertGreater(checks, 0, "No conservation columns were compared")

    def test_grch37_sv_conservation_matches_vep(self):
        self._check_build("GRCh37")

    def test_grch38_sv_conservation_matches_vep(self):
        self._check_build("GRCh38")


if __name__ == "__main__":
    unittest.main()
