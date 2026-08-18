"""
The switch from IVAT to VEP altered some EvidenceKeys

The key changes happened in:

classification.migrations.0042_vep_ekeys.modifying_existing_keys
"""
import unittest
from collections import namedtuple

from django.contrib.auth.models import User
from django.test import TestCase

from library.genomics.vcf_enums import VariantClass
from sync.shariant.historical_ekey_converter import HistoricalEKeyConverter


class Test(TestCase):
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

            latest_data_only_values = {key: {'value': value.get('value')} for key, value in latest_data.items() if isinstance(value, dict)}
            only_values = {key: {'value': value.get('value')} for key, value in to_shariant_data.items() if isinstance(value, dict)}

            self.assertDictEqual(only_values, latest_data_only_values, f"{description} to_shariant conversion (latest keys)")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
