from django.test import RequestFactory, TestCase, override_settings

from classification.enums import EvidenceKeyValueType, SubmissionSource
from classification.models import Classification, EvidenceKey
from classification.tests.models.test_utils import ClassificationTestUtils
from classification.views.classification_view import ClassificationView


class ClassificationTestCaseViews(TestCase):

    def setUp(self):
        EvidenceKey.objects.create(
            key='foo',
            value_type=EvidenceKeyValueType.SELECT,
            options=[
                {'key': 'bar', 'aliases': 'BAAA'},
                {'key': 'van'}
            ]
        )

        ClassificationTestUtils.setUp()
        self.requestFactory = RequestFactory()

    def tearDown(self):
        EvidenceKey.objects.filter(key='foo').delete()
        ClassificationTestUtils.tearDown()

    def request_post(self, data):
        _, user = ClassificationTestUtils.lab_and_user()
        request = self.requestFactory.get('/classification/api/classifications/v2/record/')
        request.method = 'POST'
        request.content_type = 'application/json'
        request.user = user
        request.data = data
        return ClassificationView().post(request=request)

    @staticmethod
    def _codes(messages) -> set[str]:
        return {m["code"] for m in messages}

    @override_settings(CLASSIFICATION_MATCH_VARIANTS=False)
    def test_return_data(self):
        lab, _user = ClassificationTestUtils.lab_and_user()
        response = self.request_post({
            "id": f"{lab.group_name}/test_123456",
            "data": {
                "c_hgvs": "NM_000059.3(BRCA2):c.3G>C",
                "genome_build": "GRCh37.p13",
                "molecular_consequence": "coding_sequence_variant, oops",
                "condition": "xxx",
                "clinical_significance": "VOUS",
                "esp_af": "0.1"
            },
            "return_data": True,
            "publish": "lab"
        }).data
        data = response["data"]

        self.assertEqual(data["clinical_significance"]["value"], "VUS")  # VOUS is an alias
        self.assertEqual(data["esp_af"]["value"], 0.1)
        self.assertEqual(data["molecular_consequence"]["value"], ["coding_sequence_variant", "oops"])
        self.assertEqual(self._codes(data["molecular_consequence"]["validation"]), {"invalid_value"})
        # Derived from c_hgvs, so the submitter can't override them
        self.assertEqual(data["gene_symbol"], {"value": "BRCA2", "immutable": "variantgrid"})
        self.assertEqual(data["refseq_transcript_id"], {"value": "NM_000059.3", "immutable": "variantgrid"})
        for mandatory_key in ("zygosity", "allele_origin"):
            self.assertIn("mandatory", self._codes(data[mandatory_key]["validation"]))
        self.assertIn("share_failure", self._codes(response["patch_messages"]))

    @override_settings(CLASSIFICATION_MATCH_VARIANTS=False)
    def test_test_mode(self):
        lab, _user = ClassificationTestUtils.lab_and_user()
        response = self.request_post({
            "id": f"{lab.group_name}/test_123456",
            "test": True,
            "data": {
                "c_hgvs": "NM_000059.3(BRCA2):c.3G>C",
                "genome_build": "GRCh37"
            },
            "publish": "lab"
        }).data

        self.assertIsNone(response["id"])
        self.assertIn("test_mode", self._codes(response["patch_messages"]))
        self.assertFalse(Classification.objects.filter(lab__group_name=lab.group_name,
                                                       lab_record_id="test_123456").exists())

    @override_settings(CLASSIFICATION_MATCH_VARIANTS=False)
    def test_bulk(self):
        lab, _user = ClassificationTestUtils.lab_and_user()
        response = self.request_post({"records": [{
            "id": f"{lab.group_name}/test_1",
            "test": True,
            "data": {
                "c_hgvs": "NM_000059.3(BRCA2):c.3G>C",
                "genome_build": "GRCh37"
            },
            "publish": "lab"
        }, {
            "id": f"{lab.group_name}/test_2",
            "test": True,
            "data": {
                "c_hgvs": "NM_000059.3(BRCA2):c.3G>C",
                "genome_build": "GRCh37"
            },
            "publish": "lab"
        }]})
        results = response.data.get('results')
        self.assertEqual(results[0]['lab_record_id'], "test_1")
        self.assertEqual(results[1]['lab_record_id'], "test_2")

    @override_settings(CLASSIFICATION_MATCH_VARIANTS=False)
    def test_basic_update(self):
        lab, _user = ClassificationTestUtils.lab_and_user()
        # all test that aliases work re foo BAAA -> bar
        response_1 = self.request_post({
            "id": f"{lab.group_name}/test_10",
            "test": False,
            "source": SubmissionSource.API.value,
            "data": {
                "c_hgvs": "xxx",
                "genome_build": "GRCh37",
                "p_hgvs": "abc",
                "consent": "zzz"
            },
            "publish": "lab"
        })
        vc_id = response_1.data.get('id')
        vc = Classification.objects.get(pk=vc_id)
        self.assertEqual(vc.evidence["p_hgvs"].get('immutable'), 'api')

        response_2 = self.request_post({
            "id": f"{lab.group_name}/test_10",
            "test": False,
            "source": SubmissionSource.FORM.value,
            "data": {
                "c_hgvs": "xxx",
                "consent": {"note": "nutty"},
                "p_hgvs": "cant change api immutable"
            },
            "publish": "lab"
        })

        vc.refresh_from_db()
        self.assertEqual(vc.get('consent'), 'zzz')

        self.assertEqual(vc.evidence["p_hgvs"].get('immutable'), 'api')
        self.assertEqual(vc.evidence["consent"].get('note'), 'nutty')
        self.assertEqual(vc.evidence["consent"].get('immutable'), 'api')
