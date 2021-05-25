"""
The switch from IVAT to VEP altered some EvidenceKeys

The key changes happened in:

classification.migrations.0042_vep_ekeys.modifying_existing_keys
"""
from collections import namedtuple
from unittest import skip

from django.contrib.auth.models import User
import unittest

from annotation.models.models_enums import VariantClass
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter
from classification.models.evidence_key import EvidenceKey, EvidenceCategory, EvidenceKeyValueType


class Test(unittest.TestCase):
    LATEST_DATA_VALUE_KEYS = {
        "cadd": {"value": "1.5"},
        "1000_genomes_af": {"value": "0.42123"},
        "esp_af": {"value": "0.42123"},
        "exac_af": {"value": "0.42123"},
        "flossies_af": {"value": "0.42123"},
        "gnomad_af": {"value": "0.42123"},
        "gnomad_popmax": {"value": 'AFR'},
        "gnomad_popmax_af": {"value": "0.42123"},
        "grantham": {"value": "5"},
        "uk10k_af": {"value": "0.42123"},
        "sample_type": {"value": "blood"},
    }

    def setUp(self):
        self.maxDiff = None

    @skip
    def test_historical_keys(self):
        # No longer required
        """ SA Path keys using old IVAT annotation as of August 2019 deployment """

        FakeVCM = namedtuple('FakeVCM', ['user'])
        existing_keys = list(EvidenceKey.objects.all())
        EvidenceKey.objects.all().delete()
        try:
            SAPATH_AUG_2019_EKEYS = [
                ("cadd", EvidenceCategory.COMPUTATIONAL_AND_PREDICTIVE_DATA, EvidenceKeyValueType.FREE_ENTRY),
                ("1000_genomes", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.FREE_ENTRY),
                ("esp", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.FREE_ENTRY),
                ("exac", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.PERCENT),
                ("flossies", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.FREE_ENTRY),
                ("gnomad", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.PERCENT),
                ("gnomad_hp", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.SELECT),
                ("gnomad_hp_maf", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.PERCENT),
                ("grantham", EvidenceCategory.COMPUTATIONAL_AND_PREDICTIVE_DATA, EvidenceKeyValueType.FREE_ENTRY),
                ("sample", EvidenceCategory.HEADER_TEST, EvidenceKeyValueType.SELECT),
                ("uk10k", EvidenceCategory.POPULATION_DATA, EvidenceKeyValueType.PERCENT),
                ("variant_type", EvidenceCategory.COMPUTATIONAL_AND_PREDICTIVE_DATA, EvidenceKeyValueType.MULTISELECT),
            ]

            evidence_keys = []
            for key, evidence_category, value_type in SAPATH_AUG_2019_EKEYS:
                ekey = EvidenceKey(key=key, evidence_category=evidence_category, value_type=value_type)
                evidence_keys.append(ekey)

            EvidenceKey.objects.bulk_create(evidence_keys)

            HISTORICAL_TEST_DATA = {
                "cadd": "1.5",
                "1000_genomes": "0.42123",
                "esp": "42.123",
                "exac": "42.123",
                "flossies": "42.123",
                "gnomad": "42.123",
                "gnomad_hp": 'AFR',
                "gnomad_hp_maf": "42.123",
                "grantham": "5",
                "sample": "blood",
                "uk10k": "0.42123",
            }
            HISTORICAL_TEST_DATA_VALUE_KEYS = {
                "cadd": {"value": "1.5"},
                "1000_genomes": {"value": "0.42123"},
                "esp": {"value": "42.123"},
                "exac": {"value": "42.123"},
                "flossies": {"value": "42.123"},
                "gnomad": {"value": "42.123"},
                "gnomad_hp": {"value": 'AFR'},
                "gnomad_hp_maf": {"value": "42.123"},
                "grantham": {"value": "5"},
                "sample": {"value": "blood"},
                "uk10k": {"value": "0.42123"},
            }

            DATA = {
                "plain dict": (HISTORICAL_TEST_DATA, self.LATEST_DATA_VALUE_KEYS),
                "value keys": (HISTORICAL_TEST_DATA_VALUE_KEYS, self.LATEST_DATA_VALUE_KEYS),
            }

            for description, (historical_test_data, shariant_test_data) in DATA.items():
                # add variant_type (will be removed)
                h_data = historical_test_data.copy()
                h_data["variant_type"] = "indel"

                historical_converter = HistoricalEKeyConverter()
                non_legacy_vcm = FakeVCM(user=User("fake_user"))
                s_data = historical_converter.to_shariant(non_legacy_vcm, h_data, add_note=False)
                self.assertDictEqual(s_data, shariant_test_data, f"{description} to_shariant conversion")
        finally:
            # Restore existing keys
            EvidenceKey.objects.all().delete()
            EvidenceKey.objects.bulk_create(existing_keys)

    def test_latest_keys(self):
        """ Should do nothing - failure probably due to historical_ekey_converter not updated to handle DB EKeys """

        FakeVCM = namedtuple('FakeVCM', ['user'])
        DATA = {
            "value keys": self.LATEST_DATA_VALUE_KEYS,
        }

        for description, latest_test_data in DATA.items():
            latest_data = latest_test_data.copy()
            latest_data["variant_class"] = VariantClass.SNV

            historical_converter = HistoricalEKeyConverter()
            non_legacy_vcm = FakeVCM(user=User("fake_user"))
            to_shariant_data = historical_converter.to_shariant(non_legacy_vcm, latest_data)
            self.assertDictEqual(to_shariant_data, latest_data, f"{description} to_shariant conversion (latest keys)")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
