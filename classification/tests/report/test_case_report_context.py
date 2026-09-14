from datetime import timedelta

from django.test import TestCase
from django.utils import timezone

from classification.enums import SpecialEKeys
from classification.enums.classification_enums import SomaticClinicalSignificance
from classification.report.case_report_context import (
    Alteration,
    ReportContext,
    ReportVariant,
    ReportVariantKind,
    amp_tier,
    build_gene_groups,
    build_tier_groups,
    context_as_dict,
    sort_report_variants,
)
from classification.report.template_validation import FIXTURE_CONTEXT
from classification.tests.report.fake_report_variants import (
    FakeModification,
    fake_gene_level_variant,
    fake_report_variant,
)

TIER_1 = SomaticClinicalSignificance.TIER_1
TIER_2 = SomaticClinicalSignificance.TIER_2
TIER_3 = SomaticClinicalSignificance.TIER_3


class AmpTierTest(TestCase):
    """ The sub-tier the printed report needs is the tier and the AMP level together, and where the
        levels don't say which sub-tier it is the report gets the bare tier plus a warning """

    @staticmethod
    def _amp_tier(tier, *levels):
        values = {SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE: tier}
        for level in levels:
            values[f"amp:level_{level.lower()}"] = "met"
        return amp_tier(FakeModification(values=values))

    def test_tier_1_with_level_a_is_ia(self):
        self.assertEqual(self._amp_tier(TIER_1, "A"), ("IA", []))

    def test_tier_1_with_level_b_only_is_ib(self):
        self.assertEqual(self._amp_tier(TIER_1, "B"), ("IB", []))

    def test_level_a_wins_over_b(self):
        self.assertEqual(self._amp_tier(TIER_1, "A", "B"), ("IA", []))

    def test_tier_2_with_level_c_is_iic(self):
        self.assertEqual(self._amp_tier(TIER_2, "C"), ("IIC", []))

    def test_tier_2_with_level_d_only_is_iid(self):
        self.assertEqual(self._amp_tier(TIER_2, "D"), ("IID", []))

    def test_tier_3_and_4_have_no_sub_tier(self):
        self.assertEqual(self._amp_tier(TIER_3), ("III", []))
        self.assertEqual(self._amp_tier(SomaticClinicalSignificance.TIER_4), ("IV", []))

    def test_tier_1_with_no_matching_level_warns(self):
        tier, warnings = self._amp_tier(TIER_1, "C")
        self.assertEqual(tier, "I")
        self.assertEqual(len(warnings), 1)

    def test_tier_2_with_no_matching_level_warns(self):
        tier, warnings = self._amp_tier(TIER_2, "A")
        self.assertEqual(tier, "II")
        self.assertEqual(len(warnings), 1)

    def test_tier_1_or_2_warns(self):
        tier, warnings = self._amp_tier(SomaticClinicalSignificance.TIER_1_OR_2)
        self.assertEqual(tier, "I/II")
        self.assertEqual(len(warnings), 1)

    def test_untiered_is_blank_with_no_warning(self):
        self.assertEqual(amp_tier(FakeModification(values={})), ("", []))


class AmplificationMagnitudeTest(TestCase):
    """ copy_number is the caller's absolute count, fold_change its ratio against the normal - an
        amplification may carry either or both, and a template prints what it has """

    def test_both_reach_the_context_the_templates_read(self):
        amplification = fake_report_variant("ZGENE", copy_number=12, fold_change=5.22,
                                            variant=fake_gene_level_variant("<GAIN:HGNC:9>"))

        as_dict = context_as_dict(ReportContext(
            source_level="S", variants=[amplification], kind_groups=[], tier_groups=[],
            gene_groups=[]))["variants"][0]

        self.assertEqual(12, as_dict["copy_number"])
        self.assertEqual(5.22, as_dict["fold_change"])

    def test_a_ratio_only_call_has_no_copy_number(self):
        amplification = fake_report_variant("ZGENE", fold_change=4.31428,
                                            variant=fake_gene_level_variant("<GAIN:HGNC:9>"))
        self.assertIsNone(amplification.copy_number)
        self.assertEqual(4.31428, amplification.fold_change)


class ReportVariantKindTest(TestCase):

    def test_gene_level_alts_say_what_the_event_is(self):
        fusion = fake_report_variant("GENE1", variant=fake_gene_level_variant(
            "<FUSION:HGNC:1>", gene_symbols=["GENE1", "GENE2"]))
        self.assertEqual(fusion.kind, ReportVariantKind.FUSION)
        self.assertEqual(fusion.alteration, Alteration.FUSION)

        gain = fake_report_variant("GENE3", variant=fake_gene_level_variant("<GAIN:HGNC:3>"))
        self.assertEqual(gain.kind, ReportVariantKind.COPY_NUMBER)
        self.assertEqual(gain.alteration, Alteration.AMPLIFICATION)

    def test_a_fusion_is_filed_under_its_five_prime_partner(self):
        """ Both partners travel so a template can print GENE1::GENE2, but the report sorts on the
            anchor - the gene the Variant itself is filed under """
        fusion = fake_report_variant("GENE1", variant=fake_gene_level_variant(
            "<FUSION:HGNC:1>", gene_symbols=["GENE1", "GENE2"]))
        self.assertEqual(fusion.gene_symbol, "GENE1")
        self.assertEqual(fusion.gene_symbols, ["GENE1", "GENE2"])
        self.assertEqual(fusion.gene_label, "GENE1::GENE2")

    def test_an_unmatched_classification_falls_back_to_variant_class(self):
        record = FakeModification(values={SpecialEKeys.GENE_SYMBOL: "GENE5",
                                          SpecialEKeys.VARIANT_CLASS: "copy_number_gain"})
        variant = ReportVariant.build(record, user=None, evidence={})
        self.assertEqual(variant.kind, ReportVariantKind.COPY_NUMBER)
        self.assertEqual(variant.alteration, Alteration.AMPLIFICATION)

    def test_anything_else_is_a_small_variant(self):
        variant = fake_report_variant("GENE6")
        self.assertEqual(variant.kind, ReportVariantKind.SMALL_VARIANT)
        self.assertEqual(variant.alteration, Alteration.VARIANT)


class ReportOrderTest(TestCase):

    def test_kind_then_tier_then_gene_then_vaf(self):
        small_low = fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"], vaf=0.1, pk=1)
        small_high = fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"], vaf=0.6, pk=2)
        small_tier_2 = fake_report_variant("AGENE", tier=TIER_2, amp_levels=["C"], vaf=0.9, pk=3)
        small_b_gene = fake_report_variant("BGENE", tier=TIER_1, amp_levels=["A"], vaf=0.9, pk=4)
        amplification = fake_report_variant("ZGENE", tier=TIER_2, amp_levels=["C"], copy_number=11,
                                            variant=fake_gene_level_variant("<GAIN:HGNC:9>"), pk=5)
        fusion = fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"], pk=6,
                                     variant=fake_gene_level_variant("<FUSION:HGNC:1>",
                                                                     gene_symbols=["AGENE", "BGENE"]))

        ordered = sort_report_variants([fusion, amplification, small_tier_2, small_b_gene,
                                        small_low, small_high])

        self.assertEqual([v.modification.pk for v in ordered],
                         [small_high.modification.pk, small_low.modification.pk,
                          small_b_gene.modification.pk, small_tier_2.modification.pk,
                          amplification.modification.pk, fusion.modification.pk])

    def test_the_bigger_fold_change_prints_first_when_neither_has_a_copy_number(self):
        """ A caller that writes only a ratio still orders the amplifications in a gene """
        smaller = fake_report_variant("ZGENE", tier=TIER_2, amp_levels=["C"], fold_change=3.6,
                                      variant=fake_gene_level_variant("<GAIN:HGNC:9>"), pk=1)
        bigger = fake_report_variant("ZGENE", tier=TIER_2, amp_levels=["C"], fold_change=5.2,
                                     variant=fake_gene_level_variant("<GAIN:HGNC:9>"), pk=2)

        ordered = sort_report_variants([smaller, bigger])

        self.assertEqual([v.modification.pk for v in ordered],
                         [bigger.modification.pk, smaller.modification.pk])

    def test_sub_tier_orders_ahead_of_the_bare_tier(self):
        ia = fake_report_variant("ZGENE", tier=TIER_1, amp_levels=["A"], pk=1)
        ib = fake_report_variant("AGENE", tier=TIER_1, amp_levels=["B"], pk=2)
        iic = fake_report_variant("AGENE", tier=TIER_2, amp_levels=["C"], pk=3)
        untiered = fake_report_variant("AGENE", pk=4)

        ordered = sort_report_variants([untiered, iic, ib, ia])
        self.assertEqual([v.amp_tier for v in ordered], ["IA", "IB", "IIC", ""])


class TierGroupTest(TestCase):

    def test_an_amplification_prints_under_its_tier_beside_the_small_variants(self):
        small = fake_report_variant("AGENE", tier=TIER_2, amp_levels=["C"], vaf=0.4, pk=1)
        amplification = fake_report_variant("ZGENE", tier=TIER_2, amp_levels=["C"], copy_number=11,
                                            variant=fake_gene_level_variant("<GAIN:HGNC:9>"), pk=2)

        tier_groups = {tg.tier: tg for tg in build_tier_groups([small, amplification])}
        tier_2_genes = [gene.gene_symbol for gene in tier_groups[TIER_2].genes]
        self.assertEqual(tier_2_genes, ["AGENE", "ZGENE"])

    def test_the_standard_tiers_are_always_printed(self):
        tier_groups = build_tier_groups([fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"])])
        self.assertEqual([tg.tier for tg in tier_groups], [TIER_1, TIER_2, TIER_3])
        self.assertEqual([tg.genes for tg in tier_groups][1:], [[], []])

    def test_an_unreported_variant_is_left_out_but_counted(self):
        """ The tier prints "no reportable variants" rather than "none detected" when the case has
            variants in it that are not going on the report """
        unreported = fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"], reported=False)
        tier_1 = build_tier_groups([unreported])[0]
        self.assertEqual(tier_1.genes, [])
        self.assertEqual(tier_1.unreported_count, 1)


class GeneGroupTest(TestCase):

    def test_the_most_recently_curated_gene_paragraph_wins_and_the_rest_are_named(self):
        older = fake_report_variant("AGENE", gene_summary="older text", pk=1,
                                    modified=timezone.now() - timedelta(days=1))
        newer = fake_report_variant("AGENE", gene_summary="newer text", pk=2, modified=timezone.now())

        gene_group = build_gene_groups([older, newer])[0]
        self.assertEqual(gene_group.gene_summary, "newer text")
        self.assertEqual(gene_group.gene_summary_source, newer.modification.pk)
        self.assertEqual(len(gene_group.warnings), 1)

    def test_agreeing_gene_paragraphs_raise_no_warning(self):
        first = fake_report_variant("AGENE", gene_summary="same text", pk=1)
        second = fake_report_variant("AGENE", gene_summary="same text", pk=2)
        self.assertEqual(build_gene_groups([first, second])[0].warnings, [])

    def test_gene_groups_expose_the_evidence_dicts_templates_read(self):
        """ gene_groups is the shape the sapath#246 templates were written against """
        gene_group = build_gene_groups([fake_report_variant("AGENE", c_hgvs="c.1A>G")])[0]
        self.assertEqual(gene_group.classifications[0]["c_hgvs"]["value"], "c.1A>G")


class FixtureContextTest(TestCase):
    """ FIXTURE_CONTEXT is a hand written stand-in for context_as_dict, because the validation
        module cannot import the model that imports it - so it has to be kept in step by hand """

    def test_the_validation_fixture_has_the_same_keys_as_a_real_context(self):
        real = context_as_dict(ReportContext(source_level="S", variants=[], kind_groups=[],
                                             tier_groups=[], gene_groups=[]))
        self.assertEqual(sorted(FIXTURE_CONTEXT), sorted(real))
