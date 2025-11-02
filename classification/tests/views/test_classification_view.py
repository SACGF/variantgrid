import json

from deepdiff import DeepDiff
from django.test import TestCase, RequestFactory, override_settings

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

    @override_settings(CLASSIFICATION_MATCH_VARIANTS=False)
    @override_settings(ALLELE_ORIGIN_NOT_PROVIDED_BUCKET="U")
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
        # pop all the values which are database id/time based
        # all remaining values should be the same everytime
        pop_me = ["id", "flag_collection", "lab_record_id", "last_edited", "published_version", "title", "version", "resolved_condition", "allele"]
        for pop_key in pop_me:
            response.pop(pop_key)

        expected = {
            "institution_name": "InstX",
            "lab_id": "instx/labby",
            "cr_lab_id": "test_123456",
            "org_name": "InstX",
            "lab_name": "Labby",
            "publish_level": "lab",
            "version_publish_level": "lab",
            "version_is_published": None,
            "is_last_published": True,
            "can_write": True,
            "can_write_latest": True,
            "clinical_context": None,
            "data": {
                "owner": {
                    "value": "joejoe"
                },
                "c_hgvs": {
                    "value": "NM_000059.3(BRCA2):c.3G>C",
                    "immutable": "variantgrid"
                },
                "esp_af": {
                    "value": 0.1,
                    "immutable": "api"
                },
                "zygosity": {
                    "immutable": "api",
                    "validation": [
                        {
                            "code": "mandatory",
                            "message": "Missing mandatory value",
                            "severity": "error"
                        }
                    ]
                },
                "condition": {
                    "value": "xxx",
                    "immutable": "api"
                },
                "gene_symbol": {
                    "value": "BRCA2",
                    "immutable": "variantgrid"
                },
                "genome_build": {
                    "value": "GRCh37.p13",
                    "immutable": "variantgrid"
                },
                "allele_origin": {
                    "immutable": "api",
                    "validation": [
                        {
                            "code": "mandatory",
                            "message": "Missing mandatory value",
                            "severity": "error"
                        }
                    ]
                },
                "refseq_transcript_id": {
                    "value": "NM_000059.3",
                    "immutable": "variantgrid"
                },
                "clinical_significance": {
                    "value": "VUS",
                    "immutable": "api"
                },
                "molecular_consequence": {
                    "value": [
                        "coding_sequence_variant",
                        "oops"
                    ],
                    "immutable": "api",
                    "validation": [
                        {
                            "code": "invalid_value",
                            "message": "Illegal value (oops)",
                            "severity": "error"
                        }
                    ]
                }
            },
            "withdrawn": False,
            "allele_origin_bucket": "U",
            "condition_text_match": None,
            "has_changes": False,
            "messages": [
                {
                    "code": "mandatory",
                    "message": "Missing mandatory value",
                    "severity": "error",
                    "key": "allele_origin"
                },
                {
                    "code": "invalid_value",
                    "message": "Illegal value (oops)",
                    "severity": "error",
                    "key": "molecular_consequence"
                },
                {
                    "code": "mandatory",
                    "message": "Missing mandatory value",
                    "severity": "error",
                    "key": "zygosity"
                }
            ],
            "patch_messages": [
                {
                    "code": "patched",
                    "message": "Patched changed values for allele_origin, c_hgvs, clinical_significance, condition, esp_af, gene_symbol, genome_build, molecular_consequence, owner, refseq_transcript_id, zygosity"
                },
                {
                    "code": "share_failure",
                    "message": "Cannot share record with errors"
                }
            ]
        }

        diffs = DeepDiff(t1=expected, t2=response)
        if diffs:
            print(json.dumps(response))
        self.assertFalse(diffs)

    @override_settings(CLASSIFICATION_MATCH_VARIANTS=False)
    @override_settings(ALLELE_ORIGIN_NOT_PROVIDED_BUCKET="U")
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
        })
        # now in test mode we always return all data (for the sake of useful information when testing)
        response_json = response.data

        for ignore in ['data', 'allele', 'cr_lab_id']:
            response_json.pop(ignore)

        expected = {
            "id": None,
            "lab_record_id": "test_123456",
            "institution_name": "InstX",
            "lab_id": "instx/labby",
            "org_name": "InstX",
            "lab_name": "Labby",
            "title": "instx/labby/test_123456",
            "publish_level": "lab",
            "version_publish_level": "lab",
            "version_is_published": None,
            "is_last_published": None,
            "published_version": None,
            "can_write": False,
            "can_write_latest": False,
            "clinical_context": None,
            "resolved_condition": None,
            "withdrawn": False,
            "allele_origin_bucket": "U",
            "condition_text_match": None,
            "messages": [
                {
                    "severity": "error",
                    "code": "mandatory",
                    "message": "Missing mandatory value",
                    "key": "allele_origin"
                },
                {
                    "key": "clinical_significance",
                    "severity": "error",
                    "code": "missing_significance",
                    "message": "Classification requires a value"
                },
                {
                    "severity": "error",
                    "code": "mandatory",
                    "message": "Missing mandatory value",
                    "key": "condition"
                },
                {
                    "severity": "error",
                    "code": "mandatory",
                    "message": "Missing mandatory value",
                    "key": "zygosity"
                }
            ],
            "patch_messages": [
                {
                    "code": "patched",
                    "message": "Patched changed values for allele_origin, c_hgvs, clinical_significance, condition, gene_symbol, genome_build, owner, refseq_transcript_id, zygosity"
                },
                {
                    "code": "test_mode",
                    "message": "Test mode on, no changes have been saved"
                }
            ]
        }

        diffs = DeepDiff(t1=expected, t2=response_json)
        if diffs:
            print(json.dumps(response_json, indent=4))
        self.assertFalse(diffs)

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
