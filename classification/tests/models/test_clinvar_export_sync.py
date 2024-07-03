from unittest import TestCase

from classification.models import ClinVarExportRequest
from classification.models.clinvar_export_models import ClinVarExportBatch
from classification.models.clinvar_export_sync import ClinVarExportSync, ClinVarResponseOutcome
from snpdb.models import ClinVarKey


class TestClinVarExportSync(TestCase):

    def test_submitted(self):
        clinvar_key, _ = ClinVarKey.objects.get_or_create(
            id="clinvar_test_key"
        )
        batch = ClinVarExportBatch.objects.create(
            clinvar_key=clinvar_key,
            submission_version=1
        )
        request = ClinVarExportRequest.objects.create(
            submission_batch=batch,
            url="",
            response_status_code=200,
            response_json={
                "actions": [{
                    "id": "SUB14374415-1",
                    "status": "submitted",
                    "updated": "2024-04-12T05:31:07.474814Z",
                    "targetDb": "clinvar",
                    "responses": []
                }]
            }
        )

        status = ClinVarExportSync()._handle_polling(request)
        self.assertEquals(status, ClinVarResponseOutcome.ASK_AGAIN_LATER)
