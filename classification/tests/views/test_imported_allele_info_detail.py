from django.test import TestCase

from classification.models import ImportedAlleleInfo
from classification.models.classification_variant_info_models import ResolvedVariantInfo
from classification.views.imported_allele_info_view import _hgvs_resolution_steps
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion


class HGVSResolutionStepsTest(TestCase):
    """ The detail modal's three c.HGVS rows (imported, normalised, liftover) - built in memory only """

    @staticmethod
    def _allele_info(imported_build: GenomeBuild, imported_c_hgvs: str | None = None,
                     grch37: str | None = None, grch38: str | None = None) -> ImportedAlleleInfo:
        gbpv = GenomeBuildPatchVersion.get_unspecified_patch_version_for(imported_build)
        allele_info = ImportedAlleleInfo(imported_genome_build_patch_version=gbpv, imported_c_hgvs=imported_c_hgvs)
        allele_info.grch37 = ResolvedVariantInfo(genome_build=GenomeBuild.grch37(), resolved_hgvs=grch37) if grch37 else None
        allele_info.grch38 = ResolvedVariantInfo(genome_build=GenomeBuild.grch38(), resolved_hgvs=grch38) if grch38 else None
        return allele_info

    def test_imported_build_is_normalised_other_build_is_liftover(self):
        allele_info = self._allele_info(GenomeBuild.grch38(), imported_c_hgvs="NM_000059.3(BRCA2):c.1234A>G",
                                        grch38="NM_000059.4(BRCA2):c.1234A>G", grch37="NM_000059.4(BRCA2):c.1234A>G")
        imported, normalised, liftover = _hgvs_resolution_steps(allele_info)

        self.assertEqual(imported.label, "Imported (GRCh38)")
        self.assertIsNone(imported.variant_info)
        self.assertEqual(normalised.label, "Normalised (GRCh38)")
        self.assertIs(normalised.variant_info, allele_info.grch38)
        self.assertEqual(liftover.label, "Liftover (GRCh37)")
        self.assertIs(liftover.variant_info, allele_info.grch37)

    def test_diff_from_previous_step(self):
        allele_info = self._allele_info(GenomeBuild.grch37(), imported_c_hgvs="NM_000059.3(BRCA2):c.1234A>G",
                                        grch37="NM_000059.4(BRCA2):c.1234A>G", grch38="NM_000059.4(BRCA2):c.1234A>G")
        _, normalised, liftover = _hgvs_resolution_steps(allele_info)

        self.assertEqual(normalised.diff_from_previous, ["Different transcript version"])
        self.assertIsNone(liftover.diff_from_previous)

    def test_g_hgvs_import_has_no_imported_step(self):
        allele_info = self._allele_info(GenomeBuild.grch37(), grch37="NM_000059.4(BRCA2):c.1234A>G")
        imported, normalised, liftover = _hgvs_resolution_steps(allele_info)

        self.assertIsNone(imported)
        self.assertEqual(normalised.label, "Normalised (GRCh37)")
        self.assertEqual(liftover.label, "Liftover")
        self.assertIsNone(liftover.variant_info)
        self.assertIsNone(liftover.p_hgvs)
        self.assertEqual(liftover.html, "")
