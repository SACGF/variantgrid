from unittest.mock import patch

from django.test import TestCase

from genes.cdot_data_release import (
    find_latest_release_asset,
    get_latest_release_assets,
    parse_release_asset_url,
)

DOWNLOAD_URL = "https://github.com/SACGF/cdot/releases/download/data_v0.2.33"

# The assets published by cdot data release data_v0.2.33
RELEASE_ASSET_NAMES = [
    "cdot-0.2.33.Homo_sapiens_GRCh37_Ensembl_87.gtf.json.gz",
    "cdot-0.2.33.Homo_sapiens_GRCh37_RefSeq_105.20201022.gff.json.gz",
    "cdot-0.2.33.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz",
    "cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_116.gtf.json.gz",
    "cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_110.gff.json.gz",
    "cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz",
]

COMBO_ASSET_NAMES = [
    "cdot-0.2.33.all-builds-ensembl-grch37_grch38_t2t-chm13v2.0.json.gz",
    "cdot-0.2.33.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz",
    "cdot-0.2.33.ensembl.GRCh37.json.gz",
    "cdot-0.2.33.refseq.GRCh38.json.gz",
    "cdot-0.2.33.refseq.T2T-CHM13v2.0.json.gz",
]


def _fake_data_release(asset_names):
    return {"assets": [{"browser_download_url": f"{DOWNLOAD_URL}/{name}"} for name in asset_names]}


class CdotReleaseAssetParsingTest(TestCase):
    def test_refseq_release_with_underscores(self):
        """ 'RS_2025_08' contains the separator used between name components """
        asset = parse_release_asset_url(f"{DOWNLOAD_URL}/cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz")
        self.assertIsNotNone(asset)
        self.assertEqual("0.2.33", asset.cdot_version)
        self.assertEqual("GRCh38", asset.genome_build)
        self.assertEqual("RefSeq", asset.annotation_consortium)
        self.assertEqual("RS_2025_08", asset.release)
        self.assertEqual("gff", asset.file_type)

    def test_refseq_dotted_release(self):
        asset = parse_release_asset_url(
            f"{DOWNLOAD_URL}/cdot-0.2.33.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz")
        self.assertEqual("GRCh37", asset.genome_build)
        self.assertEqual("105.20220307", asset.release)

    def test_ensembl_gtf(self):
        asset = parse_release_asset_url(f"{DOWNLOAD_URL}/cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_116.gtf.json.gz")
        self.assertEqual("Ensembl", asset.annotation_consortium)
        self.assertEqual("116", asset.release)
        self.assertEqual("gtf", asset.file_type)

    def test_combo_files_are_not_release_assets(self):
        for name in COMBO_ASSET_NAMES:
            self.assertIsNone(parse_release_asset_url(f"{DOWNLOAD_URL}/{name}"), name)


class CdotReleaseAssetLookupTest(TestCase):
    def setUp(self):
        patcher = patch("genes.cdot_data_release.get_latest_data_release",
                        return_value=_fake_data_release(RELEASE_ASSET_NAMES + COMBO_ASSET_NAMES))
        self.addCleanup(patcher.stop)
        patcher.start()

    def test_only_release_assets_returned(self):
        assets = get_latest_release_assets()
        self.assertEqual(len(RELEASE_ASSET_NAMES), len(assets))

    def test_find_refseq(self):
        asset, available = find_latest_release_asset("GRCh38", "RefSeq", "RS_2025_08")
        self.assertEqual("cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz", asset.filename)
        self.assertEqual({"110", "RS_2025_08"}, {a.release for a in available})

    def test_find_ensembl(self):
        asset, _available = find_latest_release_asset("GRCh38", "Ensembl", "116")
        self.assertEqual("cdot-0.2.33.Homo_sapiens_GRCh38_Ensembl_116.gtf.json.gz", asset.filename)

    def test_build_and_consortium_are_case_insensitive(self):
        asset, _available = find_latest_release_asset("grch38", "refseq", "rs_2025_08")
        self.assertIsNotNone(asset)

    def test_release_not_published_reports_what_is(self):
        """ Older releases get dropped from newer cdot data releases """
        asset, available = find_latest_release_asset("GRCh38", "RefSeq", "RS_2023_10")
        self.assertIsNone(asset)
        self.assertEqual({"110", "RS_2025_08"}, {a.release for a in available})

    def test_build_with_no_release_assets(self):
        """ T2T-CHM13v2.0 has combo files only """
        asset, available = find_latest_release_asset("T2T-CHM13v2.0", "Ensembl", "2022_06")
        self.assertIsNone(asset)
        self.assertEqual([], available)

    def test_wrong_consortium_does_not_match(self):
        asset, _available = find_latest_release_asset("GRCh38", "Ensembl", "RS_2025_08")
        self.assertIsNone(asset)
