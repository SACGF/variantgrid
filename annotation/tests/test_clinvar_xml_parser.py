import os

from django.conf import settings
from django.test import TestCase

from annotation.clinvar_xml_parser_via_vcv import ClinVarXmlParserViaVCV
from annotation.models import ClinVarRecord, ClinVarRecordCollection

TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")


class TestClinVarXmlParser(TestCase):

    def _check_valid(self, clinvar_variation_id:int):
        filename = f"clinvar_variant_{clinvar_variation_id}.xml"
        xml_file = os.path.join(TEST_DATA_DIR, filename)
        records: list[ClinVarRecord]
        with open(xml_file, "rb") as xml_input:
            records = ClinVarXmlParserViaVCV.load_from_input(xml_input)

        return records

    def test_parsing_17662(self):
        records = self._check_valid(clinvar_variation_id=17662)
        clinvar_record_collection = ClinVarRecordCollection.objects.create(clinvar_variation_id=17662)
        for record in records:
            print(record)
            record.clinvar_record_collection = clinvar_record_collection
            record.save()
        # TODO validation

    def test_parsing_97006(self):
        records = self._check_valid(clinvar_variation_id=97006)
        clinvar_record_collection = ClinVarRecordCollection.objects.create(clinvar_variation_id=17662)
        for record in records:
            print(record)
            record.clinvar_record_collection = clinvar_record_collection
            record.save()

    # def test_doing_it_live_41650(self):
    #     records = ClinVarXmlParser.load_from_clinvar_id(14650).all_records
    #     clinvar_record_collection = ClinVarRecordCollection.objects.create(clinvar_variation_id=14650)
    #     for record in records:
    #         print(record)
    #         record.clinvar_record_collection = clinvar_record_collection
    #         record.save()
