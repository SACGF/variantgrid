from unittest import mock

from django.test import TestCase

from snpdb.models.models_dbsnp import DbSNP
from snpdb.models.models_genome import GenomeBuild

VALID_RESPONSE = {"primary_snapshot_data": {"placements_with_allele": []}}
WITHDRAWN_RESPONSE = {"withdrawn_snapshot_data": {"withdrawn_time": "2020-01-01"}}


def _placement(seq_id: str, assembly_name: str, is_chromosome: bool) -> dict:
    return {
        "placement_annot": {
            "seq_type": "refseq_chromosome",
            "seq_id_traits_by_assembly": [{"assembly_name": assembly_name, "is_chromosome": is_chromosome}],
        },
        "alleles": [{"hgvs": f"{seq_id}:g.1T>A"}],
    }


# rs2150111024 lists the GRCh38 fix-patch contig before the chromosome
PATCH_FIRST_RESPONSE = {"primary_snapshot_data": {"placements_with_allele": [
    _placement("NW_025791778.1", "GRCh38.p14", False),
    _placement("NC_000005.10", "GRCh38.p14", True),
]}}


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

    def test_alleles_for_genome_build_from_chromosome_not_patch(self):
        dbsnp = DbSNP(pk=2150111024, api_response=PATCH_FIRST_RESPONSE)
        with mock.patch.object(GenomeBuild, "get_build_with_patch", return_value="GRCh38.p14"):
            alleles = dbsnp.get_alleles_for_genome_build(GenomeBuild.grch38())
        self.assertEqual([a["hgvs"] for a in alleles], ["NC_000005.10:g.1T>A"])
