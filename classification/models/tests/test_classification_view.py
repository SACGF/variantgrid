from django.test import TestCase, RequestFactory, override_settings

from classification.enums import EvidenceKeyValueType, SubmissionSource
from classification.models import Classification, EvidenceKey
from classification.models.tests.test_utils import ClassificationTestUtils
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

    @override_settings(VARIANT_CLASSIFICATION_MATCH_VARIANTS=False)
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
            "can_write": True,
            "can_write_latest": True,
            "clinical_context": None,
            "data": {
                "c_hgvs": {
                    "immutable": "variantgrid",
                    "value": "NM_000059.3(BRCA2):c.3G>C"
                },
                "clinical_significance": {
                    "immutable": "api",
                    "value": "VUS"
                },
                "condition": {
                    "immutable": "api",
                    "value": "xxx"
                },
                "esp_af": {
                    "immutable": "api",
                    "value": 0.1
                },
                "gene_symbol": {
                    "immutable": "variantgrid",
                    "value": "BRCA2"
                },
                "genome_build": {
                    "immutable": "variantgrid",
                    "value": "GRCh37.p13"
                },
                "molecular_consequence": {
                    "immutable": "api",
                    "validation": [
                        {
                            "code": "invalid_value",
                            "message": "Illegal value (oops)",
                            "severity": "error"
                        }
                    ],
                    "value": [
                        "coding_sequence_variant",
                        "oops"
                    ]
                },
                "owner": {
                    "value": "joejoe"
                },
                "refseq_transcript_id": {
                    "immutable": "variantgrid",
                    "value": "NM_000059.3"
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
                }
            },
            #"flag_collection": 1069,
            "has_changes": False,
            #"id": 1068,
            "institution_name": "InstX",
            "lab_id": "instx/labby",
            "lab_name": "Labby",
            "org_name": "InstX",
            #"lab_record_id": "test_123456",
            #"last_edited": 1590471912.520569,
            "messages": [
                {
                    "code": "invalid_value",
                    "key": "molecular_consequence",
                    "message": "Illegal value (oops)",
                    "severity": "error"
                },
                {
                    "code": "mandatory",
                    "key": "zygosity",
                    "message": "Missing mandatory value",
                    "severity": "error"
                }
            ],
            "patch_messages": [
                {
                    "code": "share_failure",
                    "message": "Cannot share record with errors"
                }
            ],
            "publish_level": "lab",
            #"published_version": 1590471912.520569,
            #"title": "instx/labby/test_123456",
            #"version": 1590471912.520569,
            'is_last_published': True,
            "version_is_published": None,
            "version_publish_level": "lab",
            "withdrawn": False
        }

        self.maxDiff = None
        self.assertEqual(response, expected)

    @override_settings(VARIANT_CLASSIFICATION_MATCH_VARIANTS=False)
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
        response_json.pop('data')
        response_json.pop('allele')  # allele data changes a bit, should test elsewhere
        expected = {
            'id': None,
            'lab_record_id': 'test_123456',
            'institution_name': 'InstX',
            'lab_id': 'instx/labby',
            'lab_name': 'Labby',
            'org_name': 'InstX',
            'title': 'instx/labby/test_123456',
            'publish_level': 'lab',
            'published_version': None,
            'resolved_condition': None,
            'version_is_published': None,
            'version_publish_level': 'lab',
            'is_last_published': None,
            'can_write': False,
            'can_write_latest': False,
            'clinical_context': None,
            'withdrawn': False,
            'messages': [
                {'severity': 'error', 'code': 'mandatory', 'message': 'Missing mandatory value', 'key': 'clinical_significance'},
                {'severity': 'error', 'code': 'mandatory', 'message': 'Missing mandatory value', 'key': 'condition'},
                {'severity': 'error', 'code': 'mandatory', 'message': 'Missing mandatory value', 'key': 'zygosity'}
            ],
            'patch_messages': [{'code': 'test_mode', 'message': 'Test mode on, no changes have been saved'}]}

        self.maxDiff = None
        self.assertEqual(response_json, expected)

    @override_settings(VARIANT_CLASSIFICATION_MATCH_VARIANTS=False)
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

    @override_settings(VARIANT_CLASSIFICATION_MATCH_VARIANTS=False)
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
