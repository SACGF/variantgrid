"""
ClinGen Allele Registry API client + the exceptions it raises.

Separate from clingen_allele.py (which uses it) so the client can be swapped for an implementation
serving recorded responses - @see snpdb.tests.utils.mock_clingen_api
"""
import hashlib
import itertools
import json
import logging
import time
import uuid
from collections.abc import Callable
from functools import lru_cache

import requests
from django.conf import settings

from library.constants import MINUTE_SECS
from library.django_utils.django_file_utils import get_import_processing_filename
from library.utils import iter_fixed_chunks
from snpdb.models import ClinGenAllele
from snpdb.models.models_enums import ClinGenAlleleExternalRecordType


class ClinGenAlleleServerException(ClinGenAllele.ClinGenAlleleRegistryException):
    """ Could not contact server, or response != 200 """
    def __init__(self, url, method, status_code, response_json):
        json_str = ", ".join([f"{k}: {v}" for k, v in response_json.items()])
        msg = f"Error contacting ClinGen Allele Registry. {url=}, {method=}, {status_code=}. JSON: {json_str}"
        super().__init__(msg)
        self.status_code = status_code
        self.response_json = response_json
        self.description = json_str

    @property
    def is_unknown_reference(self):
        """ e.g. error is they don't have that particular transcript - could retry """
        if self.status_code == 500:
            if message := self.response_json.get("message"):
                return "Unknown reference" in message
        return False

    def get_fake_api_response(self):
        msg = self.args[0]
        api_response = {
            'message': msg,
            'errorType': ClinGenAllele.CLINGEN_ALLELE_SERVER_ERROR_TYPE,
            'description': self.description
        }
        return api_response


class ClinGenAlleleAPIException(ClinGenAllele.ClinGenAlleleRegistryException):
    """ API returned 200 OK, but was an error """


class ClinGenAlleleTooLargeException(ClinGenAllele.ClinGenAlleleRegistryException):
    """ Too big for ClinGen Allele Registry  """


class ClinGenAlleleRegistryUnavailableException(ClinGenAllele.ClinGenAlleleRegistryException):
    """ Registry down or unreachable (gateway error, connection error, timeout) - the same call may work later """
    def __init__(self, cause: Exception):
        super().__init__(f"ClinGen Allele Registry is unavailable, it may work if you try again later ({cause})")

    def get_fake_api_response(self):
        """ Stored as a VariantAllele.clingen_error, the server error type means it's retried next time """
        message = str(self)
        return {
            'message': message,
            'errorType': ClinGenAllele.CLINGEN_ALLELE_SERVER_ERROR_TYPE,
            'description': message,
        }


class ClinGenAlleleRegistryAPI:
    """ Manages API connections to ClinGen Allele Registry """

    override_class = None  # Tests sub in a recorded-response implementation - @see variantgrid.test_runner
    # Gateway / unavailable. 500 is also used for input errors (e.g. "Unknown reference") so isn't retried
    RETRY_STATUS_CODES = {502, 503, 504}

    @classmethod
    def instance(cls, **kwargs) -> 'ClinGenAlleleRegistryAPI':
        return (cls.override_class or cls)(**kwargs)

    def __init__(self, api_failure_output_filename=None, max_attempts=3, retry_backoff_secs=5):
        """ max_attempts/retry_backoff_secs apply to PUTs. Batch jobs can wait out an outage, a user waiting on
            a web request should pass max_attempts=1 """
        self.login = settings.CLINGEN_ALLELE_REGISTRY_LOGIN
        self.password = settings.CLINGEN_ALLELE_REGISTRY_PASSWORD
        self.max_attempts = max_attempts
        self.retry_backoff_secs = retry_backoff_secs
        # Left None unless a caller names one - the path is only worked out when there's a failure to dump.
        # Computing it here minted an empty import_processing dir per instance (#928)
        self.api_failure_output_filename = api_failure_output_filename

    def _get_api_failure_output_filename(self) -> str:
        if self.api_failure_output_filename is None:
            self.api_failure_output_filename = get_import_processing_filename(
                "failures", f"{uuid.uuid4()}.json", prefix="clingen_allele_registry")
        return self.api_failure_output_filename

    @staticmethod
    def check_api_response(api_response):
        """ Throws ClinGenAlleleAPIException if 'errorType' set """
        if error_type := api_response.get('errorType'):
            description = api_response['description']
            input_line = api_response['inputLine']
            message = f"ClinGeneAllele API Error: {error_type} ({description}) for input '{input_line}'"
            raise ClinGenAlleleAPIException(message)

    @staticmethod
    def _check_response(response: requests.Response):
        """ Throws Exception if response status code is not 200 OK """
        if response.status_code != 200:
            try:
                response_json = response.json()
            except requests.JSONDecodeError:
                # Gateway errors (502/504) come back as HTML
                response_json = {"text": response.text[:1000]}
            raise ClinGenAlleleServerException(response.url, response.request.method,
                                               response.status_code, response_json)

    @classmethod
    def _is_retryable(cls, e: Exception) -> bool:
        if isinstance(e, ClinGenAlleleServerException):
            return e.status_code in cls.RETRY_STATUS_CODES
        return isinstance(e, (requests.ConnectionError, requests.Timeout))

    @classmethod
    def _json_with_retry(cls, send_request: Callable[[], requests.Response], max_attempts=1, retry_backoff_secs=0):
        """ Retries registry-side failures (gateway errors, connection errors, timeouts) with exponential backoff,
            raising ClinGenAlleleRegistryUnavailableException once attempts run out.
            send_request is called for each attempt, so the PUT is re-signed with the current time.
            Worst case is max_attempts x the request timeout, plus the backoff between attempts. """
        for attempt in range(1, max_attempts + 1):
            try:
                response = send_request()
                cls._check_response(response)
                return response.json()
            except Exception as e:
                if not cls._is_retryable(e):
                    raise
                if attempt == max_attempts:
                    raise ClinGenAlleleRegistryUnavailableException(e) from e
                backoff = retry_backoff_secs * 2 ** (attempt - 1)
                logging.warning("ClinGen Allele Registry call failed (attempt %d/%d), retrying in %ds: %s",
                                attempt, max_attempts, backoff, e)
                time.sleep(backoff)

    def _signed_url(self, url) -> str:
        # copy/pasted from page 5 of https://reg.clinicalgenome.org/doc/AlleleRegistry_1.01.xx_api_v1.pdf
        identity = hashlib.sha1((self.login + self.password).encode('utf-8')).hexdigest()
        gb_time = str(int(time.time()))
        token = hashlib.sha1((url + identity + gb_time).encode('utf-8')).hexdigest()
        return url + '&gbLogin=' + self.login + '&gbTime=' + gb_time + '&gbToken=' + token

    def _put(self, url, data, chunk_size=None):
        """ Registers the alleles (or returns the existing ones) - repeating a PUT gives the same
            canonical allele IDs, so it's safe to retry """
        if chunk_size > settings.CLINGEN_ALLELE_REGISTRY_MAX_RECORDS:
            raise ValueError(f"ClinGen accepts a max of {settings.CLINGEN_ALLELE_REGISTRY_MAX_RECORDS} records")
        logging.debug("Calling ClinGen API")
        default_timeout = 2 * MINUTE_SECS
        if chunk_size:
            timeout = 2 * MINUTE_SECS * chunk_size / 1000
            timeout = max(default_timeout, timeout)
        else:
            timeout = default_timeout

        def send_request():
            return requests.put(self._signed_url(url), data=data, timeout=timeout)

        try:
            return self._json_with_retry(send_request, max_attempts=self.max_attempts,
                                         retry_backoff_secs=self.retry_backoff_secs)
        except Exception as e:
            api_failure = {
                "request": url,
                "timeout": timeout,
                "data": data,
            }
            api_failure_output_filename = self._get_api_failure_output_filename()
            with open(api_failure_output_filename, "w") as f:
                json.dump(api_failure, f)

            msg = f"API call failed, debug info written to '{api_failure_output_filename}'"
            if isinstance(e, ClinGenAlleleRegistryUnavailableException):
                logging.error(msg)
                raise
            raise ClinGenAllele.ClinGenAlleleRegistryException(msg) from e

    @classmethod
    def get_code(cls, code):
        url = settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + f"/allele/{code}"
        return cls.get(url)

    @classmethod
    def get_external_code(cls, er_type: ClinGenAlleleExternalRecordType, external_code):
        suffix = f"/alleles?{er_type.value}={external_code}"
        url = settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + suffix
        return cls.get(url)

    @classmethod
    @lru_cache(maxsize=1000)
    def get_hgvs(cls, hgvs_string: str):
        suffix = f"/allele?hgvs={hgvs_string}"
        url = settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + suffix
        return cls.get(url)

    @classmethod
    def get(cls, url):
        """ Single attempt - GETs back interactive lookups such as search, so fail fast """
        return cls._json_with_retry(lambda: requests.get(url, timeout=MINUTE_SECS))

    def _clingen_hgvs_put_iter(self, hgvs_iter, file_type="hgvs"):
        """ Calls ClinGen in batches
            file_type = {hgvs, id, MyVariantInfo_hg19.id, MyVariantInfo_hg38.id, ExAC.id, gnomAD.id}
         """
        url = settings.CLINGEN_ALLELE_REGISTRY_DOMAIN + f"/alleles?file={file_type}"
        chunk_size = settings.CLINGEN_ALLELE_REGISTRY_BATCH_SIZE

        for hgvs_chunk in iter_fixed_chunks(hgvs_iter, chunk_size):
            data = "\n".join(hgvs_chunk)
            yield self._put(url, data, chunk_size=chunk_size)

    def hgvs_put(self, hgvs_iter, file_type="hgvs"):
        return itertools.chain.from_iterable(self._clingen_hgvs_put_iter(hgvs_iter, file_type=file_type))
