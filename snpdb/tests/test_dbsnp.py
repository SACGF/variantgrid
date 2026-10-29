from unittest import mock

from django.test import TestCase

from snpdb.models.models_dbsnp import DbSNP

VALID_RESPONSE = {"primary_snapshot_data": {"placements_with_allele": []}}
WITHDRAWN_RESPONSE = {"withdrawn_snapshot_data": {"withdrawn_time": "2020-01-01"}}


class DbSNPTest(TestCase):
    def test_get_rs_id_above_32_bit(self):
        with mock.patch("snpdb.models.models_dbsnp.request") as mock_request:
            mock_request.return_value.json.return_value = VALID_RESPONSE
            dbsnp = DbSNP.get("rs2150111024")
        self.assertEqual(str(dbsnp), "rs2150111024")
        self.assertTrue(DbSNP.objects.filter(pk=2150111024).exists())

    def test_get_for_variant_multiple_ids_skips_withdrawn(self):
        DbSNP.objects.create(pk=1, api_response=WITHDRAWN_RESPONSE)
        DbSNP.objects.create(pk=2, api_response=VALID_RESPONSE)
        variant = mock.MagicMock()
        variant.variantannotation_set.get.return_value.dbsnp_rs_id = "rs1&rs2"
        dbsnp = DbSNP.get_for_variant(variant, None)
        self.assertEqual(dbsnp.pk, 2)

    def test_get_for_variant_all_withdrawn_raises(self):
        DbSNP.objects.create(pk=1, api_response=WITHDRAWN_RESPONSE)
        variant = mock.MagicMock()
        variant.variantannotation_set.get.return_value.dbsnp_rs_id = "rs1"
        with self.assertRaises(ValueError):
            DbSNP.get_for_variant(variant, None)
