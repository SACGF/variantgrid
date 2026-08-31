"""
Adversarial unit tests for the genes app.
Focus: find bugs via edge cases. Tests are written to assert correct behaviour;
failures indicate real bugs in the production code.

No DB needed for most tests (pure Python / mock-based).
"""
import collections
import types

from django.test import TestCase

from genes.gene_matching import tokenize_gene_symbols
from genes.hgvs.hgvs import (
    HGVSComponents,
    HGVSDisplay,
    HGVSDiff,
    hgvs_nomen_equivalent,
    hgvs_diff_description,
)
from genes.hgvs.phgvs import PHGVS
from genes.models import TranscriptVersion
from genes.transcripts_utils import get_refseq_type, transcript_is_lrg

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tv(genome_build_data, data=None):
    """Return a SimpleNamespace that looks enough like TranscriptVersion for
    method calls that only touch genome_build_data / data / tags."""
    obj = types.SimpleNamespace()
    obj.genome_build_data = genome_build_data
    obj.data = data if data is not None else {}
    return obj


# ---------------------------------------------------------------------------
# PHGVS.parse()
# ---------------------------------------------------------------------------

class TestPHGVSParse(TestCase):

    def test_basic_snp_1letter(self):
        p = PHGVS.parse("p.G123A")
        self.assertEqual(p.aa_from, "Gly")
        self.assertEqual(p.codon, "123")
        self.assertEqual(p.aa_to, "Ala")
        self.assertTrue(p.is_confirmed)
        self.assertIsNone(p.fallback)

    def test_three_letter_lowercase_input_normalized(self):
        # The regex is re.IGNORECASE, and protein_letters_1to3_extended is keyed on single
        # uppercase letters, so 3-letter names come back via the .capitalize() fallback.
        p = PHGVS.parse("p.gly123ala")
        self.assertEqual(p.aa_from, "Gly",
                         "3-letter lowercase aa_from should be normalised to title-case")
        self.assertEqual(p.aa_to, "Ala",
                         "3-letter lowercase aa_to should be normalised to title-case")

    def test_with_transcript_prefix(self):
        p = PHGVS.parse("NM_001145661.2:p.Gly123Ala")
        self.assertEqual(p.transcript, "NM_001145661.2")
        self.assertEqual(p.aa_from, "Gly")
        self.assertEqual(p.aa_to, "Ala")

    def test_unconfirmed_parentheses(self):
        p = PHGVS.parse("p.(Gly123Ala)")
        self.assertFalse(p.is_confirmed)
        self.assertEqual(p.p_dot, "p.(Gly123Ala)")

    def test_intronic_p_question(self):
        p = PHGVS.parse("p.?")
        self.assertTrue(p.intron)
        self.assertEqual(p.p_dot, "p.?")
        self.assertFalse(p.is_confirmed)

    def test_override_is_confirmed_false_to_true(self):
        p = PHGVS.parse("p.(Gly123Ala)", override_is_confirmed_to=True)
        self.assertTrue(p.is_confirmed)
        self.assertEqual(p.p_dot, "p.Gly123Ala")  # parens dropped when confirmed

    def test_override_confirmed_on_intronic_has_no_effect_on_output(self):
        # override_is_confirmed_to=True sets is_confirmed=True on the object,
        # but p.? output is controlled by intron=True, not is_confirmed.
        p = PHGVS.parse("p.?", override_is_confirmed_to=True)
        self.assertTrue(p.is_confirmed)
        self.assertEqual(p.p_dot, "p.?")  # intronic notation unchanged

    def test_synonymous_equals(self):
        p = PHGVS.parse("p.Gly123=")
        self.assertFalse(p.intron)
        self.assertIn("=", p.p_dot)

    def test_no_p_dot_is_fallback(self):
        p = PHGVS.parse("c.123A>G")
        self.assertEqual(p.fallback, "c.123A>G")
        self.assertFalse(p.intron)

    def test_without_transcript_strips_correctly(self):
        p = PHGVS.parse("NM_001:p.Gly123Ala")
        wt = p.without_transcript
        self.assertEqual(wt.transcript, "")
        self.assertEqual(wt.p_dot, "p.Gly123Ala")

    def test_full_p_hgvs_roundtrip_with_transcript(self):
        raw = "NM_001145661.2:p.Gly123Ala"
        p = PHGVS.parse(raw)
        self.assertEqual(str(p), raw)


# ---------------------------------------------------------------------------
# HGVSComponents.__init__() / parsing
# ---------------------------------------------------------------------------

class TestHGVSComponentsInit(TestCase):

    def test_basic_parse_with_gene(self):
        c = HGVSComponents("NM_001145661.2(GATA2):c.1121G>A")
        self.assertEqual(c.transcript, "NM_001145661.2")
        self.assertEqual(c.gene_symbol, "GATA2")
        self.assertEqual(c.nomen, "c.1121G>A")

    def test_basic_parse_no_gene(self):
        c = HGVSComponents("NM_001145661.2:c.1121G>A")
        self.assertEqual(c.transcript, "NM_001145661.2")
        self.assertIsNone(c.gene_symbol)
        self.assertEqual(c.nomen, "c.1121G>A")

    def test_uppercase_kind_not_matched(self):
        # HGVS_REGEX requires lowercase kind letter (c. not C.)
        # Uppercase kind silently falls through: no gene/transcript extracted.
        c = HGVSComponents("NM_001:C.123A>G")
        # nomen should be the full string (no match)
        self.assertIsNone(c.transcript)

    def test_versioned_transcript_param_overrides(self):
        """Versioned transcript param takes precedence over HGVS string transcript."""
        c = HGVSComponents("NM_001.5(BRCA1):c.123A>G", transcript="NM_009.3")
        self.assertEqual(c.transcript, "NM_009.3")
        self.assertEqual(c.gene_symbol, "BRCA1")

    def test_unversioned_transcript_param_replaced_by_hgvs_string(self):
        """An unversioned transcript param is overwritten by the HGVS string's transcript."""
        c = HGVSComponents("NM_001.5(BRCA1):c.123A>G", transcript="NM_009")
        # No '.' in "NM_009" → hgvs string wins
        self.assertEqual(c.transcript, "NM_001.5")

    def test_transcript_parts_with_version(self):
        c = HGVSComponents("NM_001145661.2:c.123A>G")
        parts = c.transcript_parts
        self.assertEqual(parts.identifier, "NM_001145661")
        self.assertEqual(parts.version, 2)

    def test_transcript_parts_without_version(self):
        c = HGVSComponents("NM_001145661:c.123A>G")
        parts = c.transcript_parts
        self.assertEqual(parts.identifier, "NM_001145661")
        self.assertIsNone(parts.version)

    def test_without_transcript_version_strips_version(self):
        c = HGVSComponents("NM_001.5(BRCA1):c.123A>G")
        stripped = c.without_transcript_version
        self.assertEqual(str(stripped), "NM_001(BRCA1):c.123A>G")

    def test_with_gene_symbol_adds_gene(self):
        c = HGVSComponents("NM_001.5:c.123A>G")
        result = c.with_gene_symbol("BRCA1")
        self.assertEqual(str(result), "NM_001.5(BRCA1):c.123A>G")

    def test_with_gene_symbol_no_transcript_returns_self(self):
        # When transcript is None (no regex match), returns self unchanged
        c = HGVSComponents(None)
        result = c.with_gene_symbol("BRCA1")
        self.assertIs(result, c)

    def test_with_transcript_version_adds_version(self):
        c = HGVSComponents("NM_001(BRCA1):c.123A>G")
        new_c = c.with_transcript_version(7)
        self.assertEqual(str(new_c), "NM_001.7(BRCA1):c.123A>G")

    def test_sort_str_numerical_order(self):
        # Numerical part of sort_str should sort c.9 before c.100
        low = HGVSComponents("NM_001.1:c.9A>G")
        high = HGVSComponents("NM_001.1:c.100A>G")
        self.assertLess(low, high)

    def test_eq_same_string(self):
        a = HGVSComponents("NM_001.2:c.123A>G")
        b = HGVSComponents("NM_001.2:c.123A>G")
        self.assertEqual(a, b)

    def test_hash_and_eq_agree(self):
        # Components identify on the string alone, so equal objects always collide in a set
        a = HGVSComponents("NM_001.2:c.123A>G")
        b = HGVSComponents("NM_001.2:c.123A>G")
        self.assertEqual(hash(a), hash(b))
        self.assertEqual(len({a, b}), 1)

    def test_sort_str_ignores_display_state(self):
        # Components order on the HGVS itself, nothing a view happened to set
        low = HGVSComponents("NM_001.1:c.9A>G")
        high = HGVSComponents("NM_001.1:c.100A>G")
        self.assertLess(low, high)
        # the display wrapper is where being normalised gets to win
        self.assertLess(HGVSDisplay(high, is_normalised=True), HGVSDisplay(low, is_normalised=False))

    def test_sort_str_no_transcript(self):
        # A bare nomen has no transcript to append
        self.assertEqual(HGVSComponents("c.123A>G").sort_str, "0000000123A>G")

    def test_genome_build_only_breaks_ties(self):
        """ The build must not be concatenated onto sort_str - a shorter key is a prefix of a longer
            one whenever a transcript is unversioned, and the suffix would flip that comparison """
        build = collections.namedtuple("FakeBuild", "pk")("GRCh38")
        unversioned = HGVSDisplay.parse("NM_033213:c.3+46C>G", genome_build=build)
        versioned = HGVSDisplay.parse("NM_033213.5:c.3+46C>G", genome_build=build)
        self.assertLess(unversioned, versioned)

    def test_genome_build_breaks_an_actual_tie(self):
        grch37 = collections.namedtuple("FakeBuild", "pk")("GRCh37")
        grch38 = collections.namedtuple("FakeBuild", "pk")("GRCh38")
        # same transcript and nomen, so sort_str ties and only the build separates them
        a = HGVSDisplay.parse("NM_001.2(BRCA1):c.123A>G", genome_build=grch37)
        b = HGVSDisplay.parse("NM_001.2:c.123A>G", genome_build=grch38)
        self.assertEqual(a.sort_str, b.sort_str)
        self.assertLess(a, b)

    def test_display_eq_includes_view_state(self):
        components = HGVSComponents("NM_001.2:c.123A>G")
        a = HGVSDisplay(components)
        b = HGVSDisplay(components, is_normalised=True)
        self.assertNotEqual(a, b)
        self.assertEqual(len({a, b, HGVSDisplay(components)}), 2)


# ---------------------------------------------------------------------------
# HGVSComponents.diff() — bug hunting
# ---------------------------------------------------------------------------

class TestHGVSComponentsDiff(TestCase):

    def test_identical_same(self):
        a = HGVSComponents("NM_001.2(BRCA1):c.123A>G")
        b = HGVSComponents("NM_001.2(BRCA1):c.123A>G")
        self.assertEqual(a.diff(b), HGVSDiff.SAME)

    def test_diff_transcript_id(self):
        a = HGVSComponents("NM_001.2:c.123A>G")
        b = HGVSComponents("NM_002.2:c.123A>G")
        self.assertIn(HGVSDiff.DIFF_TRANSCRIPT_ID, a.diff(b))
        self.assertNotIn(HGVSDiff.DIFF_TRANSCRIPT_VER, a.diff(b))

    def test_diff_transcript_ver_both_have_versions(self):
        a = HGVSComponents("NM_001.2:c.123A>G")
        b = HGVSComponents("NM_001.5:c.123A>G")
        diff = a.diff(b)
        self.assertIn(HGVSDiff.DIFF_TRANSCRIPT_VER, diff)
        self.assertNotIn(HGVSDiff.DIFF_TRANSCRIPT_ID, diff)

    def test_diff_gene_case_insensitive_not_flagged(self):
        a = HGVSComponents("NM_001.2(BRCA1):c.123A>G")
        b = HGVSComponents("NM_001.2(brca1):c.123A>G")
        diff = a.diff(b)
        self.assertNotIn(HGVSDiff.DIFF_GENE, diff)

    def test_diff_gene_actually_different(self):
        a = HGVSComponents("NM_001.2(BRCA1):c.123A>G")
        b = HGVSComponents("NM_001.2(BRCA2):c.123A>G")
        self.assertIn(HGVSDiff.DIFF_GENE, a.diff(b))

    def test_diff_gene_missing_on_one_side_not_flagged(self):
        # One side has no gene symbol → guard in diff() prevents DIFF_GENE
        a = HGVSComponents("NM_001.2(BRCA1):c.123A>G")
        b = HGVSComponents("NM_001.2:c.123A>G")
        self.assertNotIn(HGVSDiff.DIFF_GENE, a.diff(b))

    def test_diff_nomen_significant(self):
        a = HGVSComponents("NM_001.2:c.123A>G")
        b = HGVSComponents("NM_001.2:c.456A>G")
        diff = a.diff(b)
        self.assertIn(HGVSDiff.DIFF_NOMEN, diff)
        self.assertNotIn(HGVSDiff.DIFF_NOMEN_EXPANDED, diff)

    def test_diff_nomen_expanded_explicit_implicit(self):
        # "c.123del" vs "c.123delA" → only expanded diff, not significant
        a = HGVSComponents("NM_001.2:c.123del")
        b = HGVSComponents("NM_001.2:c.123delA")
        diff = a.diff(b)
        self.assertIn(HGVSDiff.DIFF_NOMEN_EXPANDED, diff)
        self.assertNotIn(HGVSDiff.DIFF_NOMEN, diff)

    def test_diff_symmetric(self):
        a = HGVSComponents("NM_001.2:c.123A>G")
        b = HGVSComponents("NM_001.5:c.123A>G")
        self.assertEqual(a.diff(b), b.diff(a))

    def test_diff_multiple_flags(self):
        a = HGVSComponents("NM_001.2(BRCA1):c.123A>G")
        b = HGVSComponents("NM_001.5(BRCA2):c.456T>C")
        diff = a.diff(b)
        self.assertIn(HGVSDiff.DIFF_TRANSCRIPT_VER, diff)
        self.assertIn(HGVSDiff.DIFF_GENE, diff)
        self.assertIn(HGVSDiff.DIFF_NOMEN, diff)


# ---------------------------------------------------------------------------
# hgvs_nomen_equivalent()
# ---------------------------------------------------------------------------

class TestHGVSNomenEquivalent(TestCase):

    def test_identical_strings(self):
        self.assertTrue(hgvs_nomen_equivalent("c.123del", "c.123del"))

    def test_del_explicit_vs_implicit_equivalent(self):
        self.assertTrue(hgvs_nomen_equivalent("c.123del", "c.123delA"))

    def test_del_count_vs_explicit_match(self):
        self.assertTrue(hgvs_nomen_equivalent("c.123del3", "c.123delACG"))

    def test_del_count_vs_explicit_mismatch(self):
        self.assertFalse(hgvs_nomen_equivalent("c.123del3", "c.123delAC"))

    def test_del_count_vs_count_different(self):
        self.assertFalse(hgvs_nomen_equivalent("c.123del3", "c.123del4"))

    def test_del_different_positions(self):
        self.assertFalse(hgvs_nomen_equivalent("c.123del", "c.456del"))

    def test_del_vs_dup_not_equivalent(self):
        self.assertFalse(hgvs_nomen_equivalent("c.123del", "c.123dup"))

    def test_dup_explicit_vs_implicit(self):
        self.assertTrue(hgvs_nomen_equivalent("c.123dup", "c.123dupA"))

    def test_ins_explicit_vs_implicit(self):
        # "c.123_124ins" vs "c.123_124insATG" — one has no sequence (implicit)
        self.assertTrue(hgvs_nomen_equivalent("c.123_124ins", "c.123_124insATG"))

    def test_ins_count_vs_explicit(self):
        self.assertTrue(hgvs_nomen_equivalent("c.123_124ins3", "c.123_124insATG"))

    def test_ins_count_vs_explicit_mismatch(self):
        self.assertFalse(hgvs_nomen_equivalent("c.123_124ins4", "c.123_124insATG"))

    def test_delins_vs_del_not_equivalent(self):
        # One is delins (has ins portion), other is plain del → not equivalent
        self.assertFalse(hgvs_nomen_equivalent("c.123del", "c.123delinsATG"))

    def test_snp_different(self):
        self.assertFalse(hgvs_nomen_equivalent("c.123A>G", "c.123A>T"))

    def test_one_none_input(self):
        # One None → regex won't match → falls to final return False
        self.assertFalse(hgvs_nomen_equivalent("c.123del", None))


# ---------------------------------------------------------------------------
# tokenize_gene_symbols()
# ---------------------------------------------------------------------------

class TestTokenizeGeneSymbols(TestCase):

    def test_comma_separated(self):
        result = tokenize_gene_symbols("BRCA1,BRCA2,TP53")
        self.assertEqual(result.valid, {"BRCA1", "BRCA2", "TP53"})

    def test_space_separated(self):
        result = tokenize_gene_symbols("BRCA1 BRCA2 TP53")
        self.assertEqual(result.valid, {"BRCA1", "BRCA2", "TP53"})

    def test_mixed_delimiters(self):
        result = tokenize_gene_symbols("BRCA1, BRCA2;TP53\tATM")
        self.assertEqual(result.valid, {"BRCA1", "BRCA2", "TP53", "ATM"})

    def test_lowercase_uppercased(self):
        result = tokenize_gene_symbols("brca1 tp53")
        self.assertEqual(result.valid, {"BRCA1", "TP53"})

    def test_duplicates_collapsed(self):
        result = tokenize_gene_symbols("BRCA1 BRCA1 BRCA1")
        self.assertEqual(len(result.valid), 1)
        self.assertIn("BRCA1", result.valid)

    def test_whitespace_only_returns_empty_set(self):
        result = tokenize_gene_symbols("   \t\n  ")
        self.assertEqual(result.valid, set())
        self.assertEqual(result.oversized, [])

    def test_none_input_does_not_crash(self):
        # clean_string(None) returns "" so this should return empty set
        result = tokenize_gene_symbols(None)
        self.assertEqual(result.valid, set())
        self.assertEqual(result.oversized, [])

    def test_token_exactly_100_chars_included(self):
        token = "A" * 100
        result = tokenize_gene_symbols(token)
        self.assertIn(token, result.valid)
        self.assertEqual(result.oversized, [])

    def test_token_101_chars_excluded(self):
        token = "A" * 101
        result = tokenize_gene_symbols(token)
        self.assertNotIn(token, result.valid)
        self.assertEqual(result.valid, set())
        self.assertEqual(result.oversized, [token])


# ---------------------------------------------------------------------------
# get_refseq_type() / transcript_is_lrg()
# ---------------------------------------------------------------------------

class TestTranscriptUtils(TestCase):

    # --- get_refseq_type ---

    def test_nm_type_is_mrna(self):
        self.assertEqual(get_refseq_type("NM_001145661"), "mRNA")

    def test_nc_type_is_genomic(self):
        self.assertEqual(get_refseq_type("NC_000001.11"), "genomic")

    def test_unknown_prefix_returns_none(self):
        self.assertIsNone(get_refseq_type("BRCA1"))

    def test_enst_prefix_returns_none_from_get_refseq_type(self):
        # ENST is not in the RefSeq prefix list
        self.assertIsNone(get_refseq_type("ENST00000001"))

    # --- transcript_is_lrg ---

    def test_lrg_uppercase_true(self):
        self.assertTrue(transcript_is_lrg("LRG_1"))
        self.assertTrue(transcript_is_lrg("LRG_199t1"))

    def test_lrg_lowercase_false(self):
        # transcript_is_lrg is case-sensitive: lowercase LRG not detected
        self.assertFalse(transcript_is_lrg("lrg_1"))


# ---------------------------------------------------------------------------
# TranscriptVersion._validate_cdna_match() — with fabricated data
# ---------------------------------------------------------------------------

class TestTranscriptVersionValidateCdnaMatch(TestCase):
    """
    Tests _validate_cdna_match() by passing a SimpleNamespace as 'self',
    populated with genome_build_data matching the method's access patterns.

    Exon tuple format: (genomic_start, genomic_end, exon_id, cdna_start, cdna_end, gap)
    """

    def _validate(self, strand, exons):
        tv = _tv({"strand": strand, "exons": exons})
        return TranscriptVersion._validate_cdna_match(tv)

    def test_valid_plus_strand_single_exon(self):
        exons = [[1000, 1100, 0, 1, 101, None]]
        self.assertEqual(self._validate("+", exons), [])

    def test_valid_plus_strand_two_exons_continuous(self):
        exons = [
            [1000, 1100, 0, 1,   101, None],
            [2000, 2100, 1, 102, 201, None],
        ]
        self.assertEqual(self._validate("+", exons), [])

    def test_first_exon_cdna_not_starting_at_1(self):
        exons = [[1000, 1100, 0, 5, 105, None]]
        errors = self._validate("+", exons)
        self.assertTrue(errors, "Should report error when first exon starts at cDNA 5 not 1")
        self.assertTrue(any("not 1" in e for e in errors))

    def test_gap_in_cdna_between_exons(self):
        exons = [
            [1000, 1100, 0, 1,   100, None],
            [2000, 2100, 1, 102, 201, None],  # position 101 skipped
        ]
        errors = self._validate("+", exons)
        self.assertTrue(errors, "Should report missing cDNA positions")
        self.assertTrue(any("missing" in e for e in errors))

    def test_cdna_overlap_between_exons(self):
        # Exon 2 starts before exon 1 ends → missing will be negative → still non-zero
        exons = [
            [1000, 1100, 0, 1,  100, None],
            [2000, 2100, 1, 99, 198, None],  # overlaps
        ]
        errors = self._validate("+", exons)
        self.assertTrue(errors, "Should report error on cDNA overlap")

    def test_valid_minus_strand_two_exons(self):
        # Minus strand: exons stored low→high genomically.
        # After reversal by _validate_cdna_match, exon_id=0 (highest genomic) comes first.
        exons = [
            [1000, 1100, 1, 102, 201, None],  # genomically first = transcript-second
            [2000, 2100, 0, 1,   101, None],  # genomically last  = transcript-first
        ]
        self.assertEqual(self._validate("-", exons), [])

    def test_minus_strand_gap_detected(self):
        exons = [
            [1000, 1100, 1, 103, 201, None],  # transcript-second; cDNA 103-201
            [2000, 2100, 0, 1,   101, None],  # transcript-first; cDNA 1-101; gap at 102
        ]
        errors = self._validate("-", exons)
        self.assertTrue(errors, "Should detect cDNA gap on minus strand")

    def test_no_exons_key_returns_no_errors(self):
        tv = _tv({"strand": "+"})  # no 'exons' key
        errors = TranscriptVersion._validate_cdna_match(tv)
        self.assertEqual(errors, [])


# ---------------------------------------------------------------------------
# TranscriptVersion.tags / canonical_tag — with fabricated data
# ---------------------------------------------------------------------------

class TestTranscriptVersionTags(TestCase):
    """
    Tests tags and canonical_tag cached_property logic using SimpleNamespace.
    """

    def _get_tags(self, genome_build_data, data=None):
        # tags reads the cdot data itself, so build the nesting it expects
        tv = _tv(genome_build_data, data)
        tv.genome_build = types.SimpleNamespace(name="GRCh38")
        tv.data = dict(tv.data, genome_builds={"GRCh38": genome_build_data})
        return TranscriptVersion.tags.func(tv)

    def _get_canonical_tag(self, tags_list):
        tv = types.SimpleNamespace()
        tv.tags = tags_list
        tv.CANONICAL_SCORES = TranscriptVersion.CANONICAL_SCORES
        return TranscriptVersion.canonical_tag.func(tv)

    def test_mane_select_in_tags(self):
        tags = self._get_tags({"tag": "MANE Select,RefSeq Select"})
        self.assertIn("MANE Select", tags)
        self.assertIn("RefSeq Select", tags)

    def test_basic_tag_stripped(self):
        tags = self._get_tags({"tag": "MANE Select,basic"})
        self.assertNotIn("basic", tags)
        self.assertIn("MANE Select", tags)

    def test_tag_with_space_after_comma_matched(self):
        # "MANE Select, RefSeq Select" splits with a leading space on the second tag,
        # which would otherwise never match CANONICAL_SCORES.
        tags = self._get_tags({"tag": "MANE Select, RefSeq Select"})
        self.assertIn("RefSeq Select", tags,
                      "Tag with space after comma should still match 'RefSeq Select'")

    def test_tag_in_genome_build_data_preferred(self):
        tags = self._get_tags(
            genome_build_data={"tag": "MANE Select"},
            data={"tag": "RefSeq Select"},
        )
        # genome_build_data takes precedence (via `or`)
        self.assertIn("MANE Select", tags)

    def test_empty_genome_build_tag_falls_back_to_data(self):
        # cdot 0.2.12 - 0.2.13 put 'tag' on the transcript, 0.2.14+ in the genome build data
        tags = self._get_tags(
            genome_build_data={"tag": ""},
            data={"tag": "RefSeq Select"},
        )
        self.assertIn("RefSeq Select", tags,
                      "Empty genome_build_data tag should fall back to data tag")

    def test_canonical_tag_mane_select_wins_over_refseq(self):
        canonical = self._get_canonical_tag(["MANE Select", "RefSeq Select"])
        self.assertEqual(canonical, "MANE Select")

    def test_canonical_tag_refseq_select_when_no_mane(self):
        canonical = self._get_canonical_tag(["RefSeq Select"])
        self.assertEqual(canonical, "RefSeq Select")

    def test_mane_underscore_variant_recognised(self):
        # CANONICAL_SCORES has both "MANE Select" and "MANE_Select"
        canonical = self._get_canonical_tag(["MANE_Select"])
        self.assertEqual(canonical, "MANE_Select")


# ---------------------------------------------------------------------------
# hgvs_diff_description()
# ---------------------------------------------------------------------------

class TestChgvsDiffDescription(TestCase):

    def test_same_returns_empty(self):
        self.assertEqual(hgvs_diff_description(HGVSDiff.SAME), [])

    def test_expanded_excluded_without_include_minor(self):
        desc = hgvs_diff_description(HGVSDiff.DIFF_NOMEN_EXPANDED, include_minor=False)
        self.assertEqual(desc, [])

    def test_expanded_included_with_include_minor(self):
        desc = hgvs_diff_description(HGVSDiff.DIFF_NOMEN_EXPANDED, include_minor=True)
        self.assertTrue(len(desc) > 0)

    def test_combined_flags_produces_multiple_descriptions(self):
        combined = HGVSDiff.DIFF_TRANSCRIPT_VER | HGVSDiff.DIFF_GENE
        desc = hgvs_diff_description(combined)
        self.assertEqual(len(desc), 2)

