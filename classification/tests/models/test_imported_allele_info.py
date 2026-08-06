from django.test import TestCase

from classification.models import ImportedAlleleInfo
from classification.models.classification_variant_info_models import ImportedAlleleInfoStatus, ResolvedVariantInfo
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion


class ImportedAlleleInfoStatusTest(TestCase):
    """ update_status is the only place ImportedAlleleInfo.status is derived - these are in memory only """

    def _allele_info(self, grch37=None, grch38=None) -> ImportedAlleleInfo:
        gbpv = GenomeBuildPatchVersion.get_unspecified_patch_version_for(GenomeBuild.grch37())
        allele_info = ImportedAlleleInfo(imported_genome_build_patch_version=gbpv)
        allele_info.grch37 = grch37
        allele_info.grch38 = grch38
        return allele_info

    def test_all_builds_resolved(self):
        allele_info = self._allele_info(grch37=ResolvedVariantInfo(), grch38=ResolvedVariantInfo())
        allele_info.update_status()
        self.assertEqual(allele_info.status, ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS)

    def test_liftover_still_running(self):
        allele_info = self._allele_info(grch37=ResolvedVariantInfo())
        allele_info.update_status()
        self.assertEqual(allele_info.status, ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD)

    def test_force_complete_after_liftover_ran(self):
        allele_info = self._allele_info(grch37=ResolvedVariantInfo())
        allele_info.update_status(force_complete=True)
        self.assertEqual(allele_info.status, ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS)

    def test_force_complete_without_imported_build(self):
        allele_info = self._allele_info()
        allele_info.update_status(force_complete=True)
        self.assertEqual(allele_info.status, ImportedAlleleInfoStatus.FAILED)
