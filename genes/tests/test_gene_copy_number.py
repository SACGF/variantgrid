"""Whole-gene copy number identity - the words accepted, the alt encoding, and lookup."""
from django.core.exceptions import ValidationError

from genes.gene_copy_number import (
    find_gene_copy_number_events_for_string,
    parse_gene_copy_number_string,
    resolve_gene_copy_number_string,
)
from genes.models import GeneCopyNumberEvent, GeneCopyNumberEventKind, GeneLevelId
from genes.tests.gene_level_test_utils import create_gene_copy_number_event
from genes.tests.test_gene_fusions import GeneFusionTestCase
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt


class TestCopyNumberString(GeneFusionTestCase):
    """ What a lab or a search box writes, and the one string written back out """

    def test_the_words_accepted_for_each_direction(self):
        expected = {
            "amplification": GeneCopyNumberEventKind.GAIN,
            "amp": GeneCopyNumberEventKind.GAIN,
            "gain": GeneCopyNumberEventKind.GAIN,
            "loss": GeneCopyNumberEventKind.LOSS,
            "deletion": GeneCopyNumberEventKind.LOSS,
            "del": GeneCopyNumberEventKind.LOSS,
        }
        for word, kind in expected.items():
            self.assertEqual(("EGFR", kind), parse_gene_copy_number_string(f"EGFR {word.upper()}"), word)

    def test_written_out_as_amplification_and_loss(self):
        """ 'deletion' is input only - as output it reads as a coordinate event """
        for imported, canonical in [("EGFR gain", "EGFR amplification"),
                                    ("EGFR deletion", "EGFR loss")]:
            self.assertEqual(canonical, resolve_gene_copy_number_string(imported).canonical_str)

    def test_an_old_symbol_resolves_to_the_approved_one(self):
        self.assertEqual("SEPTIN14 amplification",
                         resolve_gene_copy_number_string("SEPT14 amp").canonical_str)

    def test_a_gene_we_do_not_know_mints_nothing(self):
        """ Otherwise any two words would become a copy number event """
        self.assertIsNone(resolve_gene_copy_number_string("NOTAGENE amplification"))
        self.assertIsNone(resolve_gene_copy_number_string("not a copy number event"))

    def test_resolves_to_the_coordinate_the_loader_writes(self):
        """ A classification target and the loader put the same coordinate through the pipeline, so
            both arrive at one Variant """
        from_loader = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN,
                                                    resolver=self.resolver)
        resolved = resolve_gene_copy_number_string("EGFR amplification")
        self.assertEqual(from_loader.variant.coordinate, resolved.variant_coordinate)

    def test_resolving_a_string_creates_no_variant(self):
        """ The Variant is the insert pipeline's to make - @see snpdb.gene_level_variants """
        before = GeneCopyNumberEvent.objects.count()
        self.assertIsNotNone(resolve_gene_copy_number_string("EGFR amplification"))
        self.assertEqual(before, GeneCopyNumberEvent.objects.count())


class TestCopyNumberVariant(GeneFusionTestCase):

    def test_the_alt_repeats_the_gene_the_position_names(self):
        event = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN,
                                              resolver=self.resolver)
        egfr = self.hgnc_ids["EGFR"]
        self.assertEqual(egfr, event.variant.locus.position)
        self.assertEqual(f"<{GeneLevelSymbolicAlt.GAIN}:{GeneIdNamespace.HGNC}:{egfr}>",
                         event.variant.alt.seq)

    def test_gain_and_loss_of_one_gene_are_different_variants(self):
        gain = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN,
                                             resolver=self.resolver)
        loss = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.LOSS,
                                             resolver=self.resolver)
        self.assertNotEqual(gain.variant, loss.variant)

    def test_a_gene_hgnc_lacks_uses_the_local_namespace(self):
        event = create_gene_copy_number_event("RP11-458D21.5", GeneCopyNumberEventKind.GAIN,
                                              resolver=self.resolver)
        self.assertGreaterEqual(event.gene.pk, GeneLevelId.CUSTOM_ID_START)
        self.assertIn(f":{GeneIdNamespace.GENE}:", event.variant.alt.seq)

    def test_clean_rejects_an_alt_that_says_something_else(self):
        event = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN,
                                              resolver=self.resolver)
        event.kind = GeneCopyNumberEventKind.LOSS
        self.assertRaises(ValidationError, event.clean)


class TestCopyNumberLookup(GeneFusionTestCase):
    """ find_gene_copy_number_events_for_string is what search runs on, so it must never create
        anything """

    def test_finds_an_existing_event(self):
        event = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN,
                                              resolver=self.resolver)
        self.assertEqual([event], find_gene_copy_number_events_for_string("EGFR amp"))

    def test_the_other_direction_is_not_it(self):
        create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN, resolver=self.resolver)
        self.assertEqual([], find_gene_copy_number_events_for_string("EGFR loss"))

    def test_resolves_the_alias(self):
        event = create_gene_copy_number_event("SEPTIN14", GeneCopyNumberEventKind.LOSS,
                                              resolver=self.resolver)
        self.assertEqual([event], find_gene_copy_number_events_for_string("SEPT14 deletion"))

    def test_creates_nothing(self):
        before = (GeneCopyNumberEvent.objects.count(), GeneLevelId.objects.count())
        self.assertEqual([], find_gene_copy_number_events_for_string("EGFR amplification"))
        self.assertEqual([], find_gene_copy_number_events_for_string("NOTAGENE amplification"))
        self.assertEqual(before, (GeneCopyNumberEvent.objects.count(), GeneLevelId.objects.count()))
