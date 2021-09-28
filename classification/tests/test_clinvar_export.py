from typing import Optional
from unittest import mock

from django.test import TestCase, override_settings

from classification.enums import SpecialEKeys, SubmissionSource, ShareLevel
from classification.models import Classification, ClinVarExport, ClinVarExportBatch, ClinVarExportStatus, \
    ClinVarExportRequestType, ClinVarExportRequest, ClinVarExportBatchStatus
from classification.models.clinvar_export_prepare import ClinvarAlleleExportPrepare
from classification.models.clinvar_export_sync import clinvar_export_sync, ClinVarResponseOutcome
from classification.models.tests.test_utils import ClassificationTestUtils
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild, ClinVarKey
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant, create_mock_allele
from uicore.json.json_types import JsonObjType


def mock_send_data(
        batch: ClinVarExportBatch,
        request_type: ClinVarExportRequestType,
        url: str,
        json_data: Optional[JsonObjType] = None):

    response_code = 200
    response_json = None

    if request_type == ClinVarExportRequestType.INITIAL_SUBMISSION:
        response_json={
            "id": "SUB999999-1",
            "responses": [],
            "status": "submitted",
            "targetDb": "clinvar",
            "updated": "2021-03-19T17:24:24.384085Z"
        }
    if request_type == ClinVarExportRequestType.POLLING_SUBMISSION:

        poll_count = batch.clinvarexportrequest_set.filter(request_type=ClinVarExportRequestType.POLLING_SUBMISSION).count()
        if poll_count == 0:
            response_json={
                "actions": [
                    {
                        "id": "SUB999999-1",
                        "responses": [
                            {
                                "files": [],
                                "message": None,
                                "objects": [
                                    {
                                        "accession": None,
                                        "content": {
                                            "clinvarProcessingStatus": "In processing",
                                            "clinvarReleaseStatus": "Not released"
                                        },
                                        "targetDb": "clinvar"
                                    }
                                ],
                                "status": "processing"
                            }
                        ],
                        "status": "processing",
                        "targetDb": "clinvar",
                        "updated": "2021-03-19T12:33:09.243072Z"
                    }
                ]
            }
        else:
            response_json={
                "actions": [
                    {
                        "id": "SUB999999-1",
                        "responses": [
                            {
                                "status": "processed",
                                "message": {
                                    "errorCode": None,
                                    "severity": "info",
                                    "text": "Your ClinVar submission processing status is \"Success\". Please find the details in the file referenced by actions[0].responses[0].files[0].url."
                                },
                                "files": [
                                    {
                                        "url": "https://dsubmit.ncbi.nlm.nih.gov/api/2.0/files/xxxxxxxx/sub999999-summary-report.json/?format=attachment"
                                    }
                                ],
                                "objects": []
                            }
                        ],
                        "status": "processed",
                        "targetDb": "clinvar",
                        "updated": "2021-03-24T04:22:04.101297Z"
                    }
                ]
            }
    elif request_type == ClinVarExportRequestType.RESPONSE_FILES:
        response_json = {
            "submissionName": "SUB673156",
            "submissionDate": "2021-03-25",
            "batchProcessingStatus": "Success",
            "batchReleaseStatus": "Not released",
            "totalCount": 1,
            "totalErrors": 0,
            "totalSuccess": 1,
            "totalPublic": 0,
            "submissions": [
                {
                    "identifiers": {
                        "localID": "ALLELE_1",  # might need to force this ID somehow
                        "clinvarLocalKey": "instx/labby/x42",
                        "localKey": "instx/labby/x42",
                        "clinvarAccession": "SCV000839746"
                    },
                    "processingStatus": "Success"
                }
            ]
        }

    return ClinVarExportRequest.objects.create(
        url=url,
        submission_batch=batch,
        request_type=request_type,
        request_json=json_data,
        response_status_code=response_code,
        response_json=response_json
    )


class TestClinVarExport(TestCase):

    @mock.patch('classification.models.clinvar_export_sync.ClinVarExportSync._send_data', side_effect=mock_send_data)
    @override_settings(VARIANT_CLASSIFICATION_MATCH_VARIANTS=False, CLINVAR_EXPORT={"enabled": True, "mode": "prod", "api_key": "ABC123"})
    def test_clinvar_setup(self, mocked_send_data):

        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        variant = slowly_create_test_variant("3", 128198980, 'A', 'T', grch37)
        allele = create_mock_allele(variant, grch37)

        ClassificationTestUtils.setUp()
        lab, user = ClassificationTestUtils.lab_and_user()

        clinvar_key = ClinVarKey.objects.create(id="test_c_key")
        lab.clinvar_key = clinvar_key
        lab.save()

        c = Classification.create(
            user=user,
            lab=lab,
            lab_record_id="x42",
            data={
                SpecialEKeys.C_HGVS: {'value': 'NM_000001.2(MADEUP):c.1913G>A'},
                SpecialEKeys.INTERPRETATION_SUMMARY: {'value': 'I have an interpretation summary'},
                SpecialEKeys.ASSERTION_METHOD: {'value': 'acmg'},
                SpecialEKeys.MODE_OF_INHERITANCE: {'value': ['autosomal_dominant']},
                SpecialEKeys.AFFECTED_STATUS: {'value': 'yes'},
                SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                SpecialEKeys.GENOME_BUILD: {'value': 'GRCh37'}
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False)

        c.publish_latest(user=admin_bot(), share_level=ShareLevel.PUBLIC)

        c.variant = variant
        c.chgvs_grch37 = "NM_000001.2(MADEUP):c.1913G>A"
        c.chgvs_grch37_full = "NM_000001.2(MADEUP):c.1913G>A"
        c.condition_resolution = {"sort_text": "ataxia-telangiectasia with generalized skin pigmentation and early death",
                                  "display_text": "MONDO:0008841 ataxia-telangiectasia with generalized skin pigmentation and early death",
                                  "resolved_join": None,
                                  "resolved_terms": [{"name": "ataxia-telangiectasia with generalized skin pigmentation and early death","term_id": "MONDO:0008841"}]}
        c.save()

        export_prepare = ClinvarAlleleExportPrepare(allele=allele)
        report = export_prepare.update_export_records()
        for report_line in report:
            print(report_line)

        clinvar_export: ClinVarExport = ClinVarExport.objects.first()
        print(clinvar_export.submission_body_validated)
        self.assertIsNotNone(clinvar_export)
        self.assertEqual(clinvar_export.status, ClinVarExportStatus.NEW_SUBMISSION)

        batches = ClinVarExportBatch.create_batches(ClinVarExport.objects.all(), force_update=True)
        self.assertEqual(len(batches), 1)

        batch = batches[0]
        _, outcome = clinvar_export_sync.next_request(batch)
        self.assertEqual(outcome, ClinVarResponseOutcome.ASK_AGAIN_LATER)
        self.assertEqual(batch.submission_identifier, "SUB999999-1")

        # next response is set to be "processing"
        _, outcome = clinvar_export_sync.next_request(batch)
        self.assertEqual(outcome, ClinVarResponseOutcome.ASK_AGAIN_LATER)

        # next response is should be to be "processed"
        _, outcome = clinvar_export_sync.next_request(batch)
        self.assertEqual(outcome, ClinVarResponseOutcome.ASK_AGAIN_NOW)
        self.assertEqual(batch.file_url, 'https://dsubmit.ncbi.nlm.nih.gov/api/2.0/files/xxxxxxxx/sub999999-summary-report.json/?format=attachment')

        _, outcome = clinvar_export_sync.next_request(batch)
        self.assertEqual(outcome, ClinVarResponseOutcome.COMPLETE)
        self.assertEqual(batch.status, ClinVarExportBatchStatus.SUBMITTED)
        clinvar_export.refresh_from_db()
        self.assertEqual(clinvar_export.scv, "SCV000839746")
