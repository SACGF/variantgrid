import os
import tempfile

from django.test import SimpleTestCase

from annotation.vep_warnings import (
    annotation_run_shortfall,
    parse_skipped_variants_file,
    parse_vep_warnings,
)

# VEP >= 112 (BaseVEP::skipped_variant_msg) - "line <N> skipped (<input line, cut to 50 chars>): <desc>"
MODERN_CONTIG_WARNING = ("WARNING: line 3 skipped (27 1000 . G T . . variant_id=99...): "
                         "Chromosome 27 not found in annotation sources or synonyms; "
                         "chromosome 27 does not overlap any features")
MODERN_OTHER_WARNING = ("WARNING: line 7 skipped (1 200 . GG A . . variant_id=1...): "
                        "Length of reference allele (GG length 2) does not match coordinates 200-200")
MODERN_UNSEEN_REASON_WARNING = "WARNING: line 11 skipped (1 300 . A <FOO> . . ...): a brand new reason"

# Pre-112 wording
LEGACY_CONTIG_WARNING = "WARNING: Chromosome 27 not found in annotation sources on line 4"
LEGACY_INCOMPLETE_WARNING = "WARNING: VCF line on line 5 looks incomplete, skipping:"
LEGACY_INCOMPLETE_INPUT_LINE = "1\t500\t.\tA\t.\t.\t.\tvariant_id=12345"


class ParseVEPWarningsTests(SimpleTestCase):
    def test_empty_and_none(self):
        for vep_warnings in (None, "", "\n\n"):
            parsed = parse_vep_warnings(vep_warnings)
            self.assertEqual(parsed.skipped_count, 0)
            self.assertEqual(parsed.skipped_contigs, set())
            self.assertEqual(parsed.incomplete_variant_ids, set())

    def test_modern_contig_not_found(self):
        parsed = parse_vep_warnings(MODERN_CONTIG_WARNING)
        self.assertEqual(parsed.skipped_line_numbers, {3})
        self.assertEqual(parsed.skipped_contigs, {"27"})

    def test_modern_count_does_not_depend_on_description(self):
        """ The count pattern stops at 'skipped' - a reason nobody has enumerated still counts """
        parsed = parse_vep_warnings(f"{MODERN_OTHER_WARNING}\n{MODERN_UNSEEN_REASON_WARNING}")
        self.assertEqual(parsed.skipped_line_numbers, {7, 11})
        self.assertEqual(parsed.skipped_contigs, set())

    def test_legacy_contig_not_found(self):
        parsed = parse_vep_warnings(LEGACY_CONTIG_WARNING)
        self.assertEqual(parsed.skipped_contigs, {"27"})
        self.assertEqual(parsed.skipped_line_numbers, {4})

    def test_legacy_incomplete_recovers_variant_id(self):
        blob = f"{LEGACY_INCOMPLETE_WARNING}\n{LEGACY_INCOMPLETE_INPUT_LINE}"
        parsed = parse_vep_warnings(blob)
        self.assertEqual(parsed.incomplete_variant_ids, {12345})
        self.assertEqual(parsed.skipped_line_numbers, {5})

    def test_mixed_blob(self):
        blob = "\n".join([
            MODERN_CONTIG_WARNING,
            MODERN_OTHER_WARNING,
            LEGACY_CONTIG_WARNING,
            LEGACY_INCOMPLETE_WARNING,
            LEGACY_INCOMPLETE_INPUT_LINE,
            "2026-08-10 - some unrelated log line",
        ])
        parsed = parse_vep_warnings(blob)
        self.assertEqual(parsed.skipped_line_numbers, {3, 4, 5, 7})
        self.assertEqual(parsed.skipped_count, 4)
        self.assertEqual(parsed.skipped_contigs, {"27"})
        self.assertEqual(parsed.incomplete_variant_ids, {12345})

    def test_duplicate_line_numbers_counted_once(self):
        parsed = parse_vep_warnings(f"{MODERN_OTHER_WARNING}\n{MODERN_OTHER_WARNING}")
        self.assertEqual(parsed.skipped_count, 1)


# printf "Line %-5s\t%-40.40s\t%s\n" - both leading fields blank-padded
SKIPPED_VARIANTS_FILE = (
    "[VEP skipped the following variants from /data/dump.vcf]\n"
    "Line 3    \t27 1000 . G T . . variant_id=99          \tChromosome 27 not found in annotation sources\n"
    "Line 7    \t1 200 . GG A . . variant_id=1            \tLength of reference allele does not match\n"
    "Line 11   \t1 300 . A <FOO> . . variant_id=2         \ta brand new reason\n"
)


class ParseSkippedVariantsFileTests(SimpleTestCase):
    def _parse(self, contents):
        with tempfile.TemporaryDirectory() as tmp_dir:
            filename = os.path.join(tmp_dir, "skipped_variants.tsv")
            with open(filename, "w") as f:
                f.write(contents)
            return parse_skipped_variants_file(filename)

    def test_parses_rows(self):
        skipped = self._parse(SKIPPED_VARIANTS_FILE)
        self.assertEqual([s.line_number for s in skipped], [3, 7, 11])
        self.assertEqual(skipped[0].line, "27 1000 . G T . . variant_id=99")
        self.assertEqual(skipped[2].description, "a brand new reason")

    def test_header_only_is_no_skips(self):
        self.assertEqual(self._parse("[VEP skipped the following variants from /data/dump.vcf]\n"), [])

    def test_empty_file(self):
        self.assertEqual(self._parse(""), [])


class AnnotationRunShortfallTests(SimpleTestCase):
    def test_complete_run(self):
        self.assertEqual(annotation_run_shortfall(100, 100, None), 0)

    def test_skips_accounted_for(self):
        blob = f"{MODERN_CONTIG_WARNING}\n{MODERN_OTHER_WARNING}"
        self.assertEqual(annotation_run_shortfall(100, 98, blob), 0)

    def test_truncated(self):
        self.assertEqual(annotation_run_shortfall(12501, 0, None), 12501)
