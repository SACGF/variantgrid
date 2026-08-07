from types import SimpleNamespace
from unittest.mock import patch, MagicMock

from django.test import TestCase, override_settings

from mme.client import MMEClient
from mme.serializers.patient_profile import classification_ontology_slots
from mme.tests.fakes import FakeClassification, make_term
from ontology.models import OntologyService


@override_settings(MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class OntologySnakeTestCase(TestCase):

    def _mondo_classification(self):
        mondo = make_term("MONDO:0007739", OntologyService.MONDO, 7739, "Huntington disease")
        return FakeClassification(terms=[mondo])

    @override_settings(MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
    def test_exact_alias_on_produces_omim_disorder(self):
        omim = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        with patch("mme.serializers.patient_profile.OntologyTermRelation.as_omim", return_value=omim):
            _features, disorders = classification_ontology_slots(self._mondo_classification())
        self.assertEqual(disorders, [{"id": "MIM:143100"}])

    @override_settings(MME_ONTOLOGY_SNAKE_EXACT=False, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
    def test_exact_alias_off_yields_no_disorders_for_mondo_only(self):
        with patch("mme.serializers.patient_profile.OntologyTermRelation.as_omim") as mock_as_omim:
            _features, disorders = classification_ontology_slots(self._mondo_classification())
        self.assertEqual(disorders, [])
        mock_as_omim.assert_not_called()

    @override_settings(MME_ONTOLOGY_SNAKE_EXACT=False, MME_ONTOLOGY_PHENOTYPE_EXPANSION=True)
    def test_phenotype_expansion_on_adds_derived_features(self):
        mondo = make_term("MONDO:0007739", OntologyService.MONDO, 7739, "Huntington disease")
        hpo_leaf = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        snake = SimpleNamespace(leaf_term=hpo_leaf)
        with patch("mme.serializers.patient_profile.OntologySnake.snake_from", return_value=[snake]):
            features, _disorders = classification_ontology_slots(FakeClassification(terms=[mondo]))
        self.assertEqual(len(features), 1)
        entry = features[0]
        self.assertEqual(entry["_derivedFrom"], "MONDO:0007739")
        self.assertIn("(via MONDO:0007739)", entry["label"])
        self.assertEqual(entry["observed"], "no")

    @override_settings(MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
    def test_phenotype_expansion_off_adds_no_features(self):
        omim = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        with patch("mme.serializers.patient_profile.OntologySnake.snake_from") as mock_snake:
            features, _disorders = classification_ontology_slots(FakeClassification(terms=[omim]))
        self.assertEqual(features, [])
        mock_snake.assert_not_called()

    @override_settings(MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False)
    def test_observed_patient_hpo_never_derived_and_dedups(self):
        # Condition HPO term also present as observed patient phenotype -> single entry, observed, no _derivedFrom
        seizure = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        condition_term = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")

        phenotype_description = MagicMock()
        phenotype_description.get_ontology_term_ids.return_value = [seizure.pk]
        ptp = SimpleNamespace(phenotype_description=phenotype_description)
        patient = SimpleNamespace(patient_text_phenotype=ptp)
        sample = SimpleNamespace(patient=patient)

        classification = FakeClassification(terms=[condition_term])
        classification.sample = sample

        with patch("mme.serializers.patient_profile.OntologyTerm.objects.filter", return_value=[seizure]):
            features, _disorders = classification_ontology_slots(classification)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["id"], "HP:0001250")
        self.assertEqual(features[0]["observed"], "yes")
        self.assertNotIn("_derivedFrom", features[0])

    @override_settings(MME_NODES={"testnode": {"base_url": "https://node.test", "token": "tok", "api_version": "1.1"}},
                       MME_DISCLAIMER="please confirm before acting")
    def test_request_envelope_includes_disclaimer(self):
        client = MMEClient("testnode")
        captured = {}

        def fake_post(url, json, headers, timeout):
            captured["body"] = json
            return SimpleNamespace(raise_for_status=lambda: None, json=lambda: {"results": []})

        with patch("mme.client.requests.post", side_effect=fake_post):
            client.match({"id": "vg:1", "features": [{"id": "HP:0001250"}]})
        self.assertEqual(captured["body"]["disclaimer"], "please confirm before acting")
