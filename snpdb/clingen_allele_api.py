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


class ClinGenAlleleRegistryAPI:
    """ Manages API connections to ClinGen Allele Registry """

    override_class = None  # Tests sub in a recorded-response implementation - @see variantgrid.test_runner

    @classmethod
    def instance(cls, **kwargs) -> 'ClinGenAlleleRegistryAPI':
        return (cls.override_class or cls)(**kwargs)

    def __init__(self, api_failure_output_filename=None):
        self.login = settings.CLINGEN_ALLELE_REGISTRY_LOGIN
        self.password = settings.CLINGEN_ALLELE_REGISTRY_PASSWORD
        if api_failure_output_filename is None:
            api_failure_output_filename = get_import_processing_filename(uuid.uuid4(),
                                                                         "api_failure_output_filename.json",
                                                                         prefix="clingen_allele_registry")
        self.api_failure_output_filename = api_failure_output_filename

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
            raise ClinGenAlleleServerException(response.url, response.request.method,
                                               response.status_code, response.json())

    def _put(self, url, data, chunk_size=None):
        if chunk_size > settings.CLINGEN_ALLELE_REGISTRY_MAX_RECORDS:
            raise ValueError(f"ClinGen accepts a max of {settings.CLINGEN_ALLELE_REGISTRY_MAX_RECORDS} records")
        logging.debug("Calling ClinGen API")
        # copy/pasted from page 5 of https://reg.clinicalgenome.org/doc/AlleleRegistry_1.01.xx_api_v1.pdf
        identity = hashlib.sha1((self.login + self.password).encode('utf-8')).hexdigest()
        gb_time = str(int(time.time()))
        token = hashlib.sha1((url + identity + gb_time).encode('utf-8')).hexdigest()
        request = url + '&gbLogin=' + self.login + '&gbTime=' + gb_time + '&gbToken=' + token
        default_timeout = 2 * MINUTE_SECS
        if chunk_size:
            timeout = 2 * MINUTE_SECS * chunk_size / 1000
            timeout = max(default_timeout, timeout)
        else:
            timeout = default_timeout

        try:
            response = requests.put(request, data=data, timeout=timeout)
            self._check_response(response)
            return response.json()
        except Exception as e:
            if self.api_failure_output_filename:
                api_failure = {
                    "request": request,
                    "timeout": timeout,
                    "data": data,
                }
                with open(self.api_failure_output_filename, "w") as f:
                    json.dump(api_failure, f)

                msg = f"API call failed, debug info written to '{self.api_failure_output_filename}'"
                raise ClinGenAllele.ClinGenAlleleRegistryException(msg) from e
            else:
                raise e

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
        response = requests.get(url, timeout=MINUTE_SECS)
        cls._check_response(response)
        return response.json()

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
