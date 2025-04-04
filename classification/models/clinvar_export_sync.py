from enum import Enum, auto
from functools import cached_property
from typing import TypedDict, Optional, Tuple

import requests
from django.conf import settings
from django.db import transaction
from requests import Response
from classification.models import ClinVarExportBatch, ClinVarExportRequest, ClinVarExportRequestType, \
    ClinVarExportBatchStatus, ClinVarExportSubmission, ClinVarExportSubmissionStatus, ClinVarExportDeleteStatus
from library.constants import MINUTE_SECS
from library.log_utils import report_message
from library.utils import JsonObjType

"""
This code is responsible for sending our data to ClinVar API (after everything else has run to work out what
data we want to send)
"""


class _ClinVarExportConfigDic(TypedDict):
    """
    How the configuration should be stored in the secrets file
    """
    mode: str
    api_key: str
    org_id: str


class ClinVarRequestExceptionType(Enum):
    OUR_DATA_ISSUE = auto()
    API_KEY_ISSUE = auto()
    TOO_MANY_REQUESTS = auto()
    CLINVAR_SERVER_ISSUE = auto()
    RESPONSE_MISSING_KEY_DATA = auto()
    NOT_SUPPORTED_YET = auto()


class ClinVarRequestException(Exception):

    def __init__(self, exception_type: ClinVarRequestExceptionType, message: Optional[str] = None):
        super().__init__(message)
        self.exception_type = exception_type
        self.message = message

    def __str__(self):
        if self.exception_type == ClinVarRequestExceptionType.OUR_DATA_ISSUE:
            return "Our data format was rejected by ClinVar - internal admin review required"
        elif self.exception_type == ClinVarRequestExceptionType.API_KEY_ISSUE:
            return "API key is not provided correctly or not authorised - internal admin review required"
        elif self.exception_type == ClinVarRequestExceptionType.TOO_MANY_REQUESTS:
            return "Too many requests, need to wait before trying again"
        elif self.exception_type == ClinVarRequestExceptionType.CLINVAR_SERVER_ISSUE:
            return "ClinVar server returned an error code, maybe wait before trying again"
        elif self.exception_type == ClinVarRequestExceptionType.RESPONSE_MISSING_KEY_DATA:
            return self.message
        elif self.exception_type == ClinVarRequestExceptionType.NOT_SUPPORTED_YET:
            return self.message

    @staticmethod
    def raise_for_status_code(code: int):
        if code == 400:
            raise ClinVarRequestException(exception_type=ClinVarRequestExceptionType.OUR_DATA_ISSUE)
        elif code == 401:
            raise ClinVarRequestException(exception_type=ClinVarRequestExceptionType.API_KEY_ISSUE)
        elif code == 429:
            raise ClinVarRequestException(exception_type=ClinVarRequestExceptionType.TOO_MANY_REQUESTS)
        elif 500 <= code <= 599:
            raise ClinVarRequestException(exception_type=ClinVarRequestExceptionType.CLINVAR_SERVER_ISSUE)
        # all good
        return True


class ClinVarResponseOutcome(Enum):
    ASK_AGAIN_NOW = auto()  # if we can immediately do the next thing in the queue for this batch
    ASK_AGAIN_LATER = auto()
    COMPLETE = auto()


class ClinVarExportSync:
    """
    In charge of
    """

    @cached_property
    def _config(self) -> _ClinVarExportConfigDic:
        return settings.CLINVAR_EXPORT or {"mode": None, "api_key": None}

    @property
    def is_enabled(self) -> bool:
        return self._config.get('mode') is not None

    @property
    def is_test(self) -> bool:
        return self._config.get('mode') != "prod"

    @property
    def submission_url(self) -> str:
        if self.is_test:
            return "https://submit.ncbi.nlm.nih.gov/apitest/v1/submissions/"
        else:
            return "https://submit.ncbi.nlm.nih.gov/api/v1/submissions/"

    @property
    def api_key(self) -> Optional[str]:
        return self._config.get('api_key')

    @property
    def org_id(self) -> Optional[str]:
        return self._config.get('org_id')

    def _send_data(self,
                   batch: ClinVarExportBatch,
                   request_type: ClinVarExportRequestType,
                   url: str,
                   json_data: Optional[JsonObjType] = None) -> ClinVarExportRequest:

        if not self.is_enabled:
            raise ValueError("ClinVarExports is not enabled in this environment")
        if not batch.clinvar_key.api_key:
            raise ValueError(f"ClinVarKey {batch.clinvar_key} does not have an API Key")

        response: Response
        headers = {
            "Content-type": "application/json",
            "SP-API-KEY": batch.clinvar_key.api_key
        }

        if json_data:
            response = requests.post(
                url=url,
                headers=headers,
                json=json_data,
                timeout=MINUTE_SECS,
            )
        else:
            response = requests.get(
                url=url,
                headers=headers,
                timeout=MINUTE_SECS,
            )

        response_json: Optional[JsonObjType] = None
        try:
            response_json = response.json()
        except ValueError:
            # ValueError implies
            pass

        return ClinVarExportRequest.objects.create(
            submission_batch=batch,
            request_type=request_type,
            request_json=json_data,
            url=url,
            response_json=response_json,
            response_status_code=response.status_code
        )

    def next_request(self, batch: ClinVarExportBatch) -> Tuple[ClinVarExportRequest, ClinVarResponseOutcome]:
        clinvar_request: ClinVarExportRequest

        # typical order would be:
        # Send our full submission to get a submission identifier
        # Send a polling request using our submission identifier until we get a file URL (likely to be told it's not ready yet)
        # When we eventually get the file URL, Request that file URL

        if unhandled_request := batch.clinvarexportrequest_set.filter(handled=False).first():
            # did we error out after getting a request before, and it's not marked as handled?
            # if so, try to handle it again now (typically after a bug fix has been deployed)
            clinvar_request = unhandled_request

        elif not batch.submission_identifier:
            # we haven't uploaded anything at all for this batch
            batch.status = ClinVarExportBatchStatus.UPLOADING
            batch.save()

            clinvar_request = self._send_data(
                batch=batch,
                request_type=ClinVarExportRequestType.INITIAL_SUBMISSION,
                json_data=batch.to_json(),
                url=self.submission_url)

        elif not batch.file_url:
            # we haven't been assigned a file yet, keep poling until we get one
            clinvar_request = self._send_data(
                batch=batch,
                request_type=ClinVarExportRequestType.POLLING_SUBMISSION,
                url=f"{self.submission_url}{batch.submission_identifier}/actions/")
        else:
            # we have a file URL, request it
            clinvar_request = self._send_data(
                batch=batch,
                request_type=ClinVarExportRequestType.RESPONSE_FILES,
                url=batch.file_url)

        handle_outcome = self.handle_request_response(clinvar_request)
        batch.refresh_from_db()  # just in case we've indirectly modified it
        return clinvar_request, handle_outcome

    def _handle_initial_submission(self, clinvar_request: ClinVarExportRequest) -> ClinVarResponseOutcome:
        status_code = clinvar_request.response_status_code
        if status_code in {200, 201}:
            if (response_json := clinvar_request.response_json) and \
                    (submission_identifier := response_json.get('id')):
                clinvar_request.submission_batch.submission_identifier = submission_identifier
                clinvar_request.submission_batch.save()
                return ClinVarResponseOutcome.ASK_AGAIN_LATER
            else:
                raise ValueError(
                    f"Unexpected response {clinvar_request.pk} code {clinvar_request.response_status_code} - {clinvar_request.response_json}")

        elif status_code == 204:
            # means it's just test? not really sure what the best thing to do about this is
            clinvar_request.submission_batch.status = ClinVarExportBatchStatus.SUBMITTED
            clinvar_request.submission_batch.save()
            return ClinVarResponseOutcome.COMPLETE
        else:
            raise ValueError(f"Processing of JSON for request {clinvar_request.pk} returned error code {clinvar_request.response_status_code} submission ID")

    def _handle_polling(self, clinvar_request: ClinVarExportRequest) -> ClinVarResponseOutcome:
        if (response_json := clinvar_request.response_json) and (actions := response_json.get("actions")):
            for action_json in actions:
                responses = action_json.get("responses")
                if responses:
                    response_json = responses[0]
                    status = response_json.get('status')
                    if status == "processing":
                        return ClinVarResponseOutcome.ASK_AGAIN_LATER
                    elif files := response_json.get("files"):
                        for file_json in files:
                            if url := file_json.get('url'):
                                clinvar_request.submission_batch.file_url = url
                                clinvar_request.submission_batch.save()
                                # don't need to wait before requesting the file
                                return ClinVarResponseOutcome.ASK_AGAIN_NOW
                else:
                    # only expect responses to be an empty array if we're submitting to test
                    if action_json.get("targetDb") == "clinvar-test":
                        return ClinVarResponseOutcome.COMPLETE
                    if status := action_json.get('status'):
                        if status == "submitted":
                            return ClinVarResponseOutcome.ASK_AGAIN_LATER

        # couldn't find "processing" or still in progress, this is unexpected
        raise ValueError(f"Processing of JSON for request {clinvar_request.pk} couldn't find status of processing, submitted or file")

    def _handle_response_file(self, clinvar_request: ClinVarExportRequest) -> ClinVarResponseOutcome:
        if response_json := clinvar_request.response_json:
            success_count: int = response_json.get("totalSuccess")
            success_delete_count: int = response_json.get("totalDeleteCount")

            submission_set = clinvar_request.submission_batch.clinvarexportsubmission_set

            # first check deletions if we have any
            if deletion_json := response_json.get("deletions"):
                for deletion in deletion_json:
                    if identifiers := deletion.get("identifiers"):
                        if delete_scv := identifiers.get("clinvarAccession"):
                            clinvar_export_submissions = list(submission_set.filter(clinvar_export__scv=delete_scv))
                            for clinvar_export_submission in clinvar_export_submissions:
                                clinvar_export_submission: ClinVarExportSubmission
                                clinvar_export_submission.response_json = identifiers
                                clinvar_export_submission.status = ClinVarExportSubmissionStatus.SUCCESS
                                clinvar_export_submission.save()

                                # note we only apply the DELETED status to the SCV for the record that was in the submission
                                # if that SCV is duplicated across our records, the others wont be updated
                                # it's up to other code and admins to stop duplicates
                                clinvar_export = clinvar_export_submission.clinvar_export
                                clinvar_export.classification_based_on = None
                                clinvar_export.delete_status = ClinVarExportDeleteStatus.DELETED
                                clinvar_export.save()

            if submissions_json := response_json.get('submissions'):
                for submission_json in submissions_json:
                    identifiers_json = submission_json.get("identifiers")

                    clinvar_export_submission: ClinVarExportSubmission
                    local_key = identifiers_json.get("clinvarLocalKey")
                    if clinvar_export_submission := submission_set.filter(localKey=local_key).first():
                        clinvar_export = clinvar_export_submission.clinvar_export
                        # TODO, do we want to do anything with the submission? e.g. around status?

                        if scv := identifiers_json.get("clinvarAccession"):
                            clinvar_export.scv = scv
                            clinvar_export.save()

                            clinvar_export_submission.scv = scv

                        new_status: Optional[ClinVarExportSubmissionStatus] = None
                        if processing_status := submission_json.get("processingStatus"):
                            if processing_status == "Success":
                                new_status = ClinVarExportSubmissionStatus.SUCCESS
                            elif processing_status == "Error":
                                new_status = ClinVarExportSubmissionStatus.ERROR
                            # expect statuses of Success or Error

                        clinvar_export_submission.response_json = submission_json
                        if new_status:
                            clinvar_export_submission.status = new_status
                        else:
                            report_message("ClinVarExportSubmission - could not determine status",
                                           level='error', extra_data={"target": clinvar_export_submission.pk})
                        clinvar_export_submission.save()
                    else:
                        report_message("ClinVarExportRequest - could not find local_key",
                                       level='error', extra_data={"target": clinvar_request.pk, "local_key": local_key})

            clinvar_request.submission_batch.status = ClinVarExportBatchStatus.SUBMITTED if (success_count > 0 or success_delete_count > 0) else ClinVarExportBatchStatus.REJECTED
            clinvar_request.submission_batch.save()

            return ClinVarResponseOutcome.COMPLETE
        else:
            raise ValueError(f"Unexpected response {clinvar_request.pk} code {clinvar_request.response_status_code} - {clinvar_request.response_json}")

    def handle_request_response(self, clinvar_request: ClinVarExportRequest) -> ClinVarResponseOutcome:
        if clinvar_request.handled:
            raise ValueError("ClinVarExport has already been handled")
        try:
            ClinVarRequestException.raise_for_status_code(clinvar_request.response_status_code)
        except ClinVarRequestException as clinvar_except:
            if clinvar_except.exception_type == ClinVarRequestExceptionType.OUR_DATA_ISSUE:
                submission_batch = clinvar_request.submission_batch
                submission_batch.status = ClinVarExportBatchStatus.REJECTED
                submission_batch.save()
            raise

        return self._handle_request_response_extra(clinvar_request)

    @transaction.atomic()
    def _handle_request_response_extra(self, clinvar_request: ClinVarExportRequest) -> ClinVarResponseOutcome:
        try:
            if clinvar_request.request_type == ClinVarExportRequestType.INITIAL_SUBMISSION:
                return self._handle_initial_submission(clinvar_request)

            elif clinvar_request.request_type == ClinVarExportRequestType.POLLING_SUBMISSION:
                polling = self._handle_polling(clinvar_request)
                if self.is_test:
                    if polling == ClinVarResponseOutcome.COMPLETE:
                        batch = clinvar_request.submission_batch
                        # TODO, make a ClinVarExportBatchStatus of TESTED?
                        batch.status = ClinVarExportBatchStatus.SUBMITTED
                        batch.save()
                return polling

            elif clinvar_request.request_type == ClinVarExportRequestType.RESPONSE_FILES:
                return self._handle_response_file(clinvar_request)

            else:
                raise ValueError(f'Unexpected request type "{clinvar_request.request_type}"')
        finally:
            clinvar_request.handled = True
            clinvar_request.save()


clinvar_export_sync = ClinVarExportSync()
