""" Gene identity on the MME wire - resolution both ways, and the published forms. """
from types import SimpleNamespace
from unittest.mock import patch

from django.test import TestCase, override_settings

from genes.models import GeneSymbol, GeneSymbolAlias
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource
from mme.genes import (
    canonical_gene_key,
    clear_gene_identity_cache,
    gene_identity_from_id,
    resolve_gene,
)
from mme.matching import find_matches
from mme.serializers.patient_profile import build_patient_profile, classification_genomic_feature
from mme.tests.fakes import (
    FakeAlleleInfo,
    FakeClassification,
    FakeResolvedVariantInfo,
    FakeSubmission,
    make_gene_version,
)

MYH11_ENSEMBL = "ENSG00000133392"
MYH11_ENTREZ = "4629"


class GeneIdentityTestCase(TestCase):
    """ MYH11 as both consortia hold it, so a symbol resolves to an ENSG and back. """

    def setUp(self):
        clear_gene_identity_cache()
        self.addCleanup(clear_gene_identity_cache)
        make_gene_version(MYH11_ENSEMBL, "MYH11", AnnotationConsortium.ENSEMBL)
        make_gene_version(MYH11_ENTREZ, "MYH11", AnnotationConsortium.REFSEQ)

    def assertIsMYH11(self, identity):
        self.assertIsNotNone(identity)
        self.assertEqual(identity.ensembl_gene_id, MYH11_ENSEMBL)
        self.assertEqual(identity.entrez_gene_id, MYH11_ENTREZ)
        self.assertEqual(identity.symbol, "MYH11")

    def test_symbol_resolves_to_ensembl_and_entrez(self):
        self.assertIsMYH11(gene_identity_from_id("MYH11"))

    def test_ensembl_id_resolves_to_the_same_identity(self):
        self.assertIsMYH11(gene_identity_from_id(MYH11_ENSEMBL))

    def test_versioned_ensembl_id_resolves(self):
        self.assertIsMYH11(gene_identity_from_id(f"{MYH11_ENSEMBL}.12"))

    def test_bare_entrez_id_resolves(self):
        self.assertIsMYH11(gene_identity_from_id(MYH11_ENTREZ))

    def test_prefixed_entrez_id_resolves(self):
        self.assertIsMYH11(gene_identity_from_id(f"GeneID:{MYH11_ENTREZ}"))

    def test_every_form_shares_a_match_key(self):
        forms = ["MYH11", "myh11", MYH11_ENSEMBL, MYH11_ENTREZ, f"GeneID:{MYH11_ENTREZ}"]
        key_sets = [gene_identity_from_id(form).match_keys() for form in forms]
        for keys in key_sets[1:]:
            self.assertEqual(keys, key_sets[0])

    def test_superseded_symbol_reaches_the_current_gene(self):
        # the same Ensembl gene, renamed at a later annotation version
        make_gene_version(MYH11_ENSEMBL, "SMMHC", AnnotationConsortium.ENSEMBL, version=2)
        identity = gene_identity_from_id("SMMHC")
        self.assertEqual(identity.ensembl_gene_id, MYH11_ENSEMBL)
        self.assertIn("MYH11", identity.match_keys())

    def test_safe_alias_reaches_the_same_identity(self):
        alias_symbol = GeneSymbol.objects.create(symbol="FAA1")
        GeneSymbolAlias.objects.create(alias="FAA1", gene_symbol=GeneSymbol.objects.get(symbol="MYH11"),
                                       source=GeneSymbolAliasSource.NCBI)
        # FAA1 holds no genes of its own, so it is a subset of MYH11 - safe to expand
        self.assertFalse(GeneSymbol.objects.get(symbol="MYH11").has_different_genes(alias_symbol))
        self.assertIn("FAA1", gene_identity_from_id("MYH11").match_keys())

    def test_alias_to_a_different_gene_stays_out_of_the_expansion(self):
        make_gene_version("ENSG00000141510", "TP53", AnnotationConsortium.ENSEMBL)
        GeneSymbolAlias.objects.create(alias="TP53", gene_symbol=GeneSymbol.objects.get(symbol="MYH11"),
                                       source=GeneSymbolAliasSource.NCBI)
        # each symbol has a gene the other does not, so this is not safe to expand
        self.assertNotIn("TP53", gene_identity_from_id("MYH11").match_keys())
        self.assertNotIn("ENSG00000141510", gene_identity_from_id("MYH11").match_keys())

    def test_unknown_identifier_does_not_resolve(self):
        self.assertIsNone(gene_identity_from_id("NOT_A_GENE"))
        self.assertIsNone(gene_identity_from_id(""))

    def test_unknown_ensembl_id_still_self_compares(self):
        identity = gene_identity_from_id("ENSG00000999999")
        self.assertEqual(identity.match_keys(), {"ENSG00000999999"})


class ResolveGeneTestCase(TestCase):
    """ Which symbol a classification is taken to be about. """

    def setUp(self):
        clear_gene_identity_cache()
        self.addCleanup(clear_gene_identity_cache)
        make_gene_version(MYH11_ENSEMBL, "MYH11", AnnotationConsortium.ENSEMBL)
        make_gene_version(MYH11_ENTREZ, "MYH11", AnnotationConsortium.REFSEQ)

    def test_refseq_curated_transcript_reaches_the_ensembl_gene(self):
        """ Labs curate against RefSeq, so the resolved symbol is re-queried for the ENSG. """
        allele_info = FakeAlleleInfo(
            grch38=FakeResolvedVariantInfo(gene_symbol=GeneSymbol.objects.get(symbol="MYH11")))
        identity = resolve_gene(FakeClassification(allele_info=allele_info))
        self.assertEqual(identity.ensembl_gene_id, MYH11_ENSEMBL)

    def test_resolved_symbol_wins_over_the_evidence_key(self):
        allele_info = FakeAlleleInfo(
            grch38=FakeResolvedVariantInfo(gene_symbol=GeneSymbol.objects.get(symbol="MYH11")))
        identity = resolve_gene(FakeClassification(gene_symbol="TP53", allele_info=allele_info))
        self.assertEqual(identity.symbol, "MYH11")

    def test_a_build_resolved_on_37_only_is_still_used(self):
        allele_info = FakeAlleleInfo(
            grch37=FakeResolvedVariantInfo(gene_symbol=GeneSymbol.objects.get(symbol="MYH11")))
        identity = resolve_gene(FakeClassification(allele_info=allele_info))
        self.assertEqual(identity.ensembl_gene_id, MYH11_ENSEMBL)

    def test_evidence_key_used_when_nothing_is_resolved(self):
        identity = resolve_gene(FakeClassification(gene_symbol="MYH11"))
        self.assertEqual(identity.ensembl_gene_id, MYH11_ENSEMBL)

    def test_unresolvable_gene_is_an_ordinary_none(self):
        self.assertIsNone(resolve_gene(FakeClassification(gene_symbol="NOT_A_GENE")))
        self.assertIsNone(resolve_gene(FakeClassification()))


@override_settings(MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"}, MME_TEST=True,
                   MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
class OutboundGeneTestCase(TestCase):
    """ What we publish in genomicFeatures[].gene. """

    def setUp(self):
        clear_gene_identity_cache()
        self.addCleanup(clear_gene_identity_cache)
        make_gene_version(MYH11_ENSEMBL, "MYH11", AnnotationConsortium.ENSEMBL)
        make_gene_version(MYH11_ENTREZ, "MYH11", AnnotationConsortium.REFSEQ)
        patcher = patch("mme.serializers.patient_profile.assert_mme_eligible")
        patcher.start()
        self.addCleanup(patcher.stop)

    @override_settings(MME_GENE_ID_PREFERENCE="ensembl")
    def test_ensembl_published_as_id_with_the_alternatives_alongside(self):
        gf = classification_genomic_feature(FakeClassification(gene_symbol="MYH11"))
        self.assertEqual(gf[0]["gene"], {
            "id": MYH11_ENSEMBL,
            "_ensemblGeneID": MYH11_ENSEMBL,
            "_entrezGeneID": MYH11_ENTREZ,
            "_geneName": "MYH11",
        })

    @override_settings(MME_GENE_ID_PREFERENCE="symbol")
    def test_symbol_preference_reproduces_the_pre_v1_1_output(self):
        gf = classification_genomic_feature(FakeClassification(gene_symbol="MYH11"))
        self.assertEqual(gf[0]["gene"], {"id": "MYH11"})

    @override_settings(MME_GENE_ID_PREFERENCE="ensembl")
    def test_unresolvable_gene_falls_back_to_the_bare_symbol(self):
        gf = classification_genomic_feature(FakeClassification(gene_symbol="NOT_A_GENE"))
        self.assertEqual(gf[0]["gene"], {"id": "NOT_A_GENE"})

    @override_settings(MME_GENE_ID_PREFERENCE="ensembl")
    def test_unresolvable_gene_still_submits(self):
        submission = FakeSubmission(classification=FakeClassification(gene_symbol="NOT_A_GENE"))
        profile = build_patient_profile(submission)
        self.assertEqual(profile["genomicFeatures"][0]["gene"], {"id": "NOT_A_GENE"})

    @override_settings(MME_GENE_ID_PREFERENCE="ensembl")
    def test_gene_with_no_ensembl_annotation_publishes_the_symbol(self):
        make_gene_version("672", "BRCA1", AnnotationConsortium.REFSEQ)
        gf = classification_genomic_feature(FakeClassification(gene_symbol="BRCA1"))
        self.assertEqual(gf[0]["gene"], {"id": "BRCA1", "_entrezGeneID": "672",
                                         "_geneName": "BRCA1"})


@override_settings(MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"},
                   MME_GENE_ID_PREFERENCE="ensembl",
                   MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
class InboundGeneMatchTestCase(TestCase):
    """ The regression this work exists for: a peer naming a gene by an identifier form we
        do not publish used to score 0 on the heaviest term in the match. """

    def setUp(self):
        clear_gene_identity_cache()
        self.addCleanup(clear_gene_identity_cache)
        make_gene_version(MYH11_ENSEMBL, "MYH11", AnnotationConsortium.ENSEMBL)
        make_gene_version(MYH11_ENTREZ, "MYH11", AnnotationConsortium.REFSEQ)

    def _find(self, gene: dict, our_symbol: str = "MYH11"):
        fake_cm = SimpleNamespace(classification=FakeClassification(gene_symbol=our_symbol))
        fake_qs = SimpleNamespace(iterator=lambda: iter([fake_cm]))
        with patch("mme.matching.mme_eligible_classifications", return_value=fake_qs):
            return find_matches({"genomicFeatures": [{"gene": gene}]})

    def test_ensembl_query_matches_a_classification_we_hold_by_symbol(self):
        matches = self._find({"id": MYH11_ENSEMBL})
        self.assertEqual(len(matches), 1)
        self.assertGreater(matches[0].score, 0)

    def test_entrez_query_matches(self):
        self.assertGreater(self._find({"id": f"GeneID:{MYH11_ENTREZ}"})[0].score, 0)

    def test_symbol_query_still_matches(self):
        self.assertGreater(self._find({"id": "MYH11"})[0].score, 0)

    def test_identifier_form_does_not_change_the_score(self):
        by_symbol = self._find({"id": "MYH11"})[0].score
        by_ensembl = self._find({"id": MYH11_ENSEMBL})[0].score
        self.assertEqual(by_symbol, by_ensembl)

    def test_extension_fields_are_read_when_a_peer_supplies_them(self):
        matches = self._find({"id": "SOMETHING_ELSE", "_ensemblGeneID": MYH11_ENSEMBL})
        self.assertGreater(matches[0].score, 0)

    def test_a_different_gene_still_does_not_match(self):
        make_gene_version("ENSG00000141510", "TP53", AnnotationConsortium.ENSEMBL)
        self.assertEqual(self._find({"id": MYH11_ENSEMBL}, our_symbol="TP53"), [])

    def test_unresolvable_genes_still_self_compare(self):
        self.assertGreater(self._find({"id": "NOT_A_GENE"}, our_symbol="NOT_A_GENE")[0].score, 0)


class CanonicalGeneKeyTestCase(TestCase):
    """ Metrics count a gene once, whichever form it was published in. """

    def setUp(self):
        clear_gene_identity_cache()
        self.addCleanup(clear_gene_identity_cache)
        make_gene_version(MYH11_ENSEMBL, "MYH11", AnnotationConsortium.ENSEMBL)
        make_gene_version(MYH11_ENTREZ, "MYH11", AnnotationConsortium.REFSEQ)

    def test_every_published_form_gives_one_key(self):
        keys = {canonical_gene_key(gene) for gene in (
            {"id": "MYH11"},
            {"id": MYH11_ENSEMBL},
            {"id": MYH11_ENSEMBL, "_geneName": "MYH11"},
            {"id": f"GeneID:{MYH11_ENTREZ}"},
        )}
        self.assertEqual(keys, {MYH11_ENSEMBL})

    def test_unresolvable_gene_keys_on_the_published_id(self):
        self.assertEqual(canonical_gene_key({"id": "not_a_gene"}), "NOT_A_GENE")
