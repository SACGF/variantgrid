from datetime import timedelta
from unittest import mock

from django.test import TestCase
from django.utils import timezone

from annotation.clinvar_fetch_request import ClinVarFetchRequest
from annotation.clinvar_xml_parser import ClinVarXmlParserOutput
from annotation.clinvar_xml_parser_via_vcv import ClinVarXmlParserViaVCV
from annotation.models import ClinVarRecordCollection

CLINVAR_VARIATION_ID = 987654321


class ClinVarFetchRequestTest(TestCase):
    """ The ClinVar fetch runs outside the transaction holding the collection's row lock """

    @staticmethod
    def _output() -> ClinVarXmlParserOutput:
        return ClinVarXmlParserOutput(all_records=[], urls=["https://example.com/vcv"])

    def test_fetch_stores_response(self):
        with mock.patch.object(ClinVarXmlParserViaVCV, "load_from_clinvar_id", return_value=self._output()):
            collection = ClinVarFetchRequest(clinvar_variation_id=CLINVAR_VARIATION_ID).fetch()
        collection.refresh_from_db()
        self.assertEqual(collection.urls, ["https://example.com/vcv"])
        self.assertEqual(collection.parser_version, ClinVarXmlParserViaVCV.PARSER_VERSION)
        self.assertIsNotNone(collection.last_loaded)

    def test_newer_simultaneous_fetch_not_overwritten(self):
        newer_fetch_date = timezone.now() + timedelta(minutes=1)

        def newer_fetch_stored_meanwhile(clinvar_variation_id):
            ClinVarRecordCollection.objects.filter(clinvar_variation_id=clinvar_variation_id).update(
                last_loaded=newer_fetch_date, urls=["newer"])
            return self._output()

        with mock.patch.object(ClinVarXmlParserViaVCV, "load_from_clinvar_id", side_effect=newer_fetch_stored_meanwhile):
            collection = ClinVarFetchRequest(clinvar_variation_id=CLINVAR_VARIATION_ID).fetch()
        self.assertEqual(collection.urls, ["newer"])
        self.assertEqual(collection.last_loaded, newer_fetch_date)
