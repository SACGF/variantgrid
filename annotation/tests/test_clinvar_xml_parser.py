import os

from django.conf import settings
from django.test import TestCase

from annotation.clinvar_xml_parser_via_vcv import ClinVarXmlParserViaVCV
from annotation.models import ClinVarRecord, ClinVarRecordCollection
from snpdb.enums import AlleleOriginBucket

TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")


class TestClinVarXmlParser(TestCase):
    """ Expected values are read off the ClinicalAssertion elements of the fixture XML """

    def _load_and_save(self, clinvar_variation_id: int) -> dict[str, ClinVarRecord]:
        filename = f"clinvar_variant_{clinvar_variation_id}.xml"
        xml_file = os.path.join(TEST_DATA_DIR, filename)
        with open(xml_file, "rb") as xml_input:
            records = ClinVarXmlParserViaVCV.load_from_input(xml_input)

        clinvar_record_collection = ClinVarRecordCollection.objects.create(clinvar_variation_id=clinvar_variation_id)
        for record in records:
            record.clinvar_record_collection = clinvar_record_collection
            record.save()
        saved = ClinVarRecord.objects.filter(clinvar_record_collection=clinvar_record_collection)
        return {r.record_id: r for r in saved}

    def test_parsing_17662(self):
        records = self._load_and_save(clinvar_variation_id=17662)
        self.assertEqual(len(records), 75)  # one per ClinicalAssertion (SCV)

        expert_panel = records["SCV004101425"]
        self.assertEqual(expert_panel.review_status, "reviewed by expert panel")
        self.assertEqual(expert_panel.stars, 3)
        self.assertEqual(expert_panel.clinical_significance, "P")
        self.assertEqual(expert_panel.allele_origin_bucket, AlleleOriginBucket.GERMLINE)
        self.assertEqual(expert_panel.submitter, "ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel, ClinGen")
        self.assertEqual(records["SCV000053475"].clinical_significance, "R")  # "risk factor"
        self.assertIsNone(records["SCV000086950"].clinical_significance)  # "not provided"

    def test_parsing_97006(self):
        records = self._load_and_save(clinvar_variation_id=97006)
        self.assertEqual(len(records), 7)

        self.assertEqual(records["SCV002812646"].stars, 1)  # criteria provided, single submitter
        self.assertEqual(records["SCV002812646"].clinical_significance, "P")
        tier_1 = records["SCV004697527"]  # SomaticClinicalImpact "Tier I - Strong"
        self.assertEqual(tier_1.somatic_clinical_significance, "tier_1")
        self.assertIsNone(tier_1.clinical_significance)
        self.assertEqual(tier_1.allele_origin_bucket, AlleleOriginBucket.SOMATIC)
        self.assertEqual(records["SCV004697529"].somatic_clinical_significance, "tier_3")

    # def test_doing_it_live_41650(self):
    #     records = ClinVarXmlParser.load_from_clinvar_id(14650).all_records
    #     clinvar_record_collection = ClinVarRecordCollection.objects.create(clinvar_variation_id=14650)
    #     for record in records:
    #         print(record)
    #         record.clinvar_record_collection = clinvar_record_collection
    #         record.save()
