from django.test import TestCase

from classification.models import ImportedAlleleInfo
from classification.models.classification_variant_info_models import (
    ImportedAlleleInfoStatus,
    ImportedAlleleInfoValidation,
    ResolvedVariantInfo,
)
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


class ImportedAlleleInfoValidationTest(TestCase):
    """ _calculate_validation keys off what was supplied on import - a g.HGVS only submission is accepted on its
        own terms, see https://github.com/SACGF/variantgrid/issues/1063 - these are in memory only """

    G_HGVS_37 = "NC_000009.11:g.101850184A>C"
    G_HGVS_38 = "NC_000009.12:g.99087902A>C"
    C_HGVS_37 = "NM_000059.3(BRCA2):c.1234A>G"
    C_HGVS_38 = "NM_000059.4(BRCA2):c.1234A>G"

    @staticmethod
    def _resolved(genome_build: GenomeBuild, c_hgvs=None, transcript_version_id=None, variant_id=1) -> ResolvedVariantInfo:
        return ResolvedVariantInfo(genome_build=genome_build, c_hgvs=c_hgvs,
                                   transcript_version_id=transcript_version_id, variant_id=variant_id)

    def _allele_info(self, grch37=None, grch38=None, **kwargs) -> ImportedAlleleInfo:
        """ imported against GRCh38, so 38 is the imported build and 37 the lifted over one """
        gbpv = GenomeBuildPatchVersion.get_unspecified_patch_version_for(GenomeBuild.grch38())
        allele_info = ImportedAlleleInfo(imported_genome_build_patch_version=gbpv,
                                         variant_coordinate="9:99087902 A>C", **kwargs)
        allele_info.grch37 = grch37
        allele_info.grch38 = grch38
        return allele_info

    def test_imported_as_c_hgvs(self):
        self.assertTrue(self._allele_info(imported_c_hgvs=self.C_HGVS_38).imported_as_c_hgvs)
        self.assertTrue(self._allele_info(imported_g_hgvs=self.G_HGVS_38,
                                          imported_transcript="NM_000059.4").imported_as_c_hgvs)
        self.assertFalse(self._allele_info(imported_g_hgvs=self.G_HGVS_38).imported_as_c_hgvs)

    def test_g_hgvs_resolving_to_genomic_form(self):
        allele_info = self._allele_info(
            imported_g_hgvs=self.G_HGVS_38,
            grch37=self._resolved(GenomeBuild.grch37(), c_hgvs=self.G_HGVS_37),
            grch38=self._resolved(GenomeBuild.grch38(), c_hgvs=self.G_HGVS_38))
        validation_tags = allele_info._calculate_validation()
        self.assertEqual(validation_tags, {})
        self.assertTrue(ImportedAlleleInfoValidation.should_include(validation_tags))

    def test_g_hgvs_resolving_to_transcript_reports_liftover_as_info(self):
        allele_info = self._allele_info(
            imported_g_hgvs=self.G_HGVS_38,
            grch37=self._resolved(GenomeBuild.grch37(), c_hgvs=self.C_HGVS_37, transcript_version_id=1),
            grch38=self._resolved(GenomeBuild.grch38(), c_hgvs="NM_000059.4(BRCA2):c.1240A>G",
                                  transcript_version_id=2))
        validation_tags = allele_info._calculate_validation()
        liftover = validation_tags["liftover"]
        self.assertTrue(liftover)
        self.assertEqual(set(liftover.values()), {"W"})
        self.assertTrue(ImportedAlleleInfoValidation.should_include(validation_tags))

    def test_unsupported_transcript_still_errors(self):
        allele_info = self._allele_info(
            imported_c_hgvs="NX_000059.4(BRCA2):c.1234A>G",
            grch37=self._resolved(GenomeBuild.grch37(), c_hgvs=self.C_HGVS_37, transcript_version_id=1),
            grch38=self._resolved(GenomeBuild.grch38(), c_hgvs=self.C_HGVS_38, transcript_version_id=2))
        validation_tags = allele_info._calculate_validation()
        self.assertEqual(validation_tags["general"]["transcript_type_not_supported"], "E")
        self.assertFalse(ImportedAlleleInfoValidation.should_include(validation_tags))

    def test_g_hgvs_build_coverage_tracks_variant(self):
        """ a g.HGVS submission only needs the variant coordinate per build, the c.HGVS is a bonus """
        both_builds = self._allele_info(
            imported_g_hgvs=self.G_HGVS_38,
            grch37=self._resolved(GenomeBuild.grch37()),
            grch38=self._resolved(GenomeBuild.grch38()))
        self.assertNotIn("builds", both_builds._calculate_validation())

        no_liftover = self._allele_info(
            imported_g_hgvs=self.G_HGVS_38,
            grch38=self._resolved(GenomeBuild.grch38()))
        self.assertEqual(no_liftover._calculate_validation()["builds"], {"missing_37": "W"})

    def test_c_hgvs_build_coverage_requires_c_hgvs(self):
        allele_info = self._allele_info(
            imported_c_hgvs=self.C_HGVS_38,
            grch37=self._resolved(GenomeBuild.grch37()),
            grch38=self._resolved(GenomeBuild.grch38(), c_hgvs=self.C_HGVS_38, transcript_version_id=2))
        self.assertEqual(allele_info._calculate_validation()["builds"], {"missing_37": "W"})
