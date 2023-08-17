import os
from typing import List
from django.test import TestCase

from django.conf import settings

from annotation.clinvar_xml_parser import ClinVarXmlParser
from annotation.models import ClinVarRecord, ClinVarRecordCollection

TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")


class TestClinVarXmlParser(TestCase):

    def test_parsing(self):
        xml_file = os.path.join(TEST_DATA_DIR, "clinvar_variant_432.xml")
        records: List[ClinVarRecord]
        with open(xml_file, "rb") as xml_input:
            records = ClinVarXmlParser.load_from_input(xml_input)

        clinvar_record_collection = ClinVarRecordCollection.objects.create(clinvar_variation_id=432)
        for record in records:
            record.clinvar_record_collection = clinvar_record_collection
            record.save()