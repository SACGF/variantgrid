""" A batch of classification records prefetches its ClinGen Allele Registry lookups (#2079) """
from unittest.mock import patch

from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from classification.models import ImportedAlleleInfo
from classification.models.classification_inserter import (
    BulkClassificationInserter,
    _new_record_clingen_lookup_hgvs,
)
from genes.fake_data import create_fake_transcript_version
from library.utils import md5sum_str
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI


class TestClassificationInserterClinGenPrefetch(TestCase):
    """ A batch of records looks up the c.HGVS it needs from ClinGen in one call (#2079). RUNX1's transcript only
        has GRCh38 data, so GRCh37 has to resolve it through ClinGen """
    C_HGVS = "ENST00000300305.7(RUNX1):c.352-1G>A"
    CLINGEN_HGVS = "ENST00000300305.7:c.352-1G>A"

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        get_fake_annotation_version(GenomeBuild.grch37())
        get_fake_annotation_version(GenomeBuild.grch38())
        create_fake_transcript_version(GenomeBuild.grch38())

    @staticmethod
    def _record(lab_record_id: str, genome_build: str, c_hgvs: str, operation="upsert", **kwargs) -> dict:
        return {"id": f"test_org/test_lab/{lab_record_id}",
                operation: {"genome_build": genome_build, "c.hgvs": {"value": c_hgvs}},
                **kwargs}

    def test_lookup_hgvs_only_for_new_records_resolved_through_clingen(self):
        records = [
            self._record("local", "GRCh38", self.C_HGVS),
            self._record("test_mode", "GRCh37", self.C_HGVS, test=True),
            self._record("patch", "GRCh37", self.C_HGVS, operation="patch"),
            self._record("clingen", "GRCh37", self.C_HGVS),
        ]
        self.assertEqual(_new_record_clingen_lookup_hgvs(records), [self.CLINGEN_HGVS])

        ImportedAlleleInfo.objects.create(
            imported_c_hgvs=self.C_HGVS, imported_md5_hash=md5sum_str(self.C_HGVS),
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                GenomeBuild.grch37()))
        self.assertEqual(_new_record_clingen_lookup_hgvs(records), [])

    def test_allele_info_resolves_from_prefetch(self):
        records = [self._record("clingen", "GRCh37", self.C_HGVS)]
        with patch.object(MockClinGenAlleleRegistryAPI, "get_hgvs") as mock_get_hgvs:
            with BulkClassificationInserter.clingen_prefetch(records):
                allele_info = ImportedAlleleInfo.get_or_create(
                    imported_c_hgvs=self.C_HGVS, imported_md5_hash=md5sum_str(self.C_HGVS),
                    imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                        GenomeBuild.grch37()))
        mock_get_hgvs.assert_not_called()
        self.assertTrue(allele_info.variant_coordinate.startswith("21:"), allele_info.variant_coordinate)
