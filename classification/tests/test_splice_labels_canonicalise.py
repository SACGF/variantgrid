"""Re-pointing splice Variants loaded before labels were canonical - @see splice_labels_canonicalise."""
from django.core.management import call_command
from django.test import TestCase

from classification.models import ImportedAlleleInfo
from genes.gene_level_resolver import GeneLevelNameResolver
from genes.gene_splice import ResolvedSpliceEvent, get_splice_event_variant
from genes.models import HGNC, GeneSymbol, HGNCImport
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import create_gene_level_variant
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion, Variant

AR_HGNC_ID = 644


class SpliceLabelsCanonicaliseTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=AR_HGNC_ID, gene_symbol_id="AR", hgnc_import=HGNCImport.objects.create(),
                            status=HGNCStatus.APPROVED, approved_name="androgen receptor")

    @staticmethod
    def _variant_with_label(label: str) -> Variant:
        """ A splice Variant as it was stored before labels were canonical """
        resolved_gene = GeneLevelNameResolver().resolve_gene("AR")
        event = ResolvedSpliceEvent(gene=resolved_gene.gene_level_id, label=label)
        return create_gene_level_variant(event.variant_coordinate)

    def test_the_variant_keeps_its_pk(self):
        """ Only the alt moves, so classifications, VariantAllele and samples follow it """
        variant = self._variant_with_label("V7")
        allele_info = ImportedAlleleInfo.objects.create(
            imported_c_hgvs="ARV7",
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                GenomeBuild.grch37()),
            variant_coordinate=f"GENE_LEVEL:{AR_HGNC_ID}-{AR_HGNC_ID} <SPLICE:HGNC:{AR_HGNC_ID}:V7>")

        call_command("splice_labels_canonicalise")

        variant.refresh_from_db()
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", variant.alt.seq)
        self.assertEqual("v_7", get_splice_event_variant(variant).label)
        allele_info.refresh_from_db()
        self.assertIn(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", allele_info.variant_coordinate,
                      "a re-match has to find the variant we have")

    def test_a_dry_run_changes_nothing(self):
        variant = self._variant_with_label("EX14SKIP")
        call_command("splice_labels_canonicalise", "--dry-run")
        variant.refresh_from_db()
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:EX14SKIP>", variant.alt.seq)

    def test_breakpoints_with_no_vcf_to_take_a_build_from_are_left(self):
        """ The build is part of what breakpoints name, and only the VCF says which one """
        variant = self._variant_with_label("X_66905968_66914514")
        call_command("splice_labels_canonicalise")
        variant.refresh_from_db()
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:X_66905968_66914514>", variant.alt.seq)

    def test_a_canonical_variant_is_left_alone(self):
        variant = self._variant_with_label("v_7")
        alt_id = variant.alt_id
        call_command("splice_labels_canonicalise")
        variant.refresh_from_db()
        self.assertEqual(alt_id, variant.alt_id)
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", variant.alt.seq)

    def test_a_canonical_variant_already_there_is_reported_rather_than_merged(self):
        old = self._variant_with_label("V7")
        canonical = self._variant_with_label("v_7")
        call_command("splice_labels_canonicalise")
        old.refresh_from_db()
        canonical.refresh_from_db()
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V7>", old.alt.seq, "left for the user to merge")
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", canonical.alt.seq)
