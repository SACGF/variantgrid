from unittest.mock import patch

from django.test import TestCase, override_settings

from annotation.fake_annotation import get_fake_annotation_version
from classification.models import ImportedAlleleInfo
from classification.models.classification_variant_info_models import (
    ImportedAlleleInfoStatus,
    ImportedAlleleInfoValidation,
    ResolvedVariantInfo,
)
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils import sha256sum_str
from snpdb.models import (
    GenomeBuild,
    GenomeBuildPatchVersion,
    Locus,
    Sequence,
    Variant,
)


class ImportedAlleleInfoValidationTagsTest(TestCase):
    def test_no_tags_renders_as_no_issues(self):
        """ the grid reads the column raw, so a record validated before it had anything to validate hits None """
        self.assertEqual(ImportedAlleleInfoValidation.validation_tags_list_from_dict(None), [])


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
    def _resolved(genome_build: GenomeBuild, resolved_hgvs=None, transcript_version_id=None, variant_id=1) -> ResolvedVariantInfo:
        return ResolvedVariantInfo(genome_build=genome_build, resolved_hgvs=resolved_hgvs,
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
            grch37=self._resolved(GenomeBuild.grch37(), resolved_hgvs=self.G_HGVS_37),
            grch38=self._resolved(GenomeBuild.grch38(), resolved_hgvs=self.G_HGVS_38))
        validation_tags = allele_info._calculate_validation()
        self.assertEqual(validation_tags, {})
        self.assertTrue(ImportedAlleleInfoValidation.should_include(validation_tags))

    def test_g_hgvs_resolving_to_transcript_reports_liftover_as_info(self):
        allele_info = self._allele_info(
            imported_g_hgvs=self.G_HGVS_38,
            grch37=self._resolved(GenomeBuild.grch37(), resolved_hgvs=self.C_HGVS_37, transcript_version_id=1),
            grch38=self._resolved(GenomeBuild.grch38(), resolved_hgvs="NM_000059.4(BRCA2):c.1240A>G",
                                  transcript_version_id=2))
        validation_tags = allele_info._calculate_validation()
        liftover = validation_tags["liftover"]
        self.assertTrue(liftover)
        self.assertEqual(set(liftover.values()), {"W"})
        self.assertTrue(ImportedAlleleInfoValidation.should_include(validation_tags))

    def _gene_level_allele_info(self, imported_c_hgvs: str, coordinate: str, **kwargs) -> ImportedAlleleInfo:
        allele_info = self._allele_info(imported_c_hgvs=imported_c_hgvs, **kwargs)
        allele_info.variant_coordinate = coordinate
        return allele_info

    def test_gene_level_is_not_a_c_hgvs_submission(self):
        """ the value arrives in imported_c_hgvs but names genes, so nothing may expect a transcript of it """
        allele_info = self._gene_level_allele_info("BRCA2::PICALM", "GENE_LEVEL:1101-1101 <FUSION:HGNC:15514>")
        self.assertTrue(allele_info.is_gene_level)
        self.assertFalse(allele_info.imported_as_c_hgvs)

    def test_gene_level_without_a_coordinate_still_reads_as_gene_level(self):
        """ a record whose gene turned out to be a typo has no coordinate to read, so the shape of the
            imported value has to answer - or it would be reported as a broken HGVS """
        allele_info = self._allele_info(imported_c_hgvs="ARHGEF::TP53")
        allele_info.variant_coordinate = None
        self.assertTrue(allele_info.is_gene_level)
        self.assertFalse(allele_info.imported_as_c_hgvs)

    def test_gene_level_that_did_not_resolve_says_so(self):
        """ the value names genes, so the error is about the genes rather than about a transcript """
        allele_info = self._allele_info(imported_c_hgvs="ARHGEF::TP53")
        allele_info.variant_coordinate = None
        general = allele_info._calculate_validation()["general"]
        self.assertEqual({"gene_level_unresolved": "E"}, general)

    def test_gene_level_resolved_in_both_builds_is_included(self):
        """ a gene-level variant sits on no transcript, so its ResolvedVariantInfo has no c.HGVS - the build
            check has to track the variant the way a g.HGVS submission does, or it can never be included """
        allele_info = self._gene_level_allele_info(
            "ARV7", "GENE_LEVEL:644-644 <SPLICE:HGNC:644:V7>",
            grch37=self._resolved(GenomeBuild.grch37()),
            grch38=self._resolved(GenomeBuild.grch38()))
        validation_tags = allele_info._calculate_validation()
        self.assertEqual(validation_tags, {})
        self.assertTrue(ImportedAlleleInfoValidation.should_include(validation_tags))

    def test_gene_level_matched_displays_as_resolved(self):
        """ a grid finding no c.HGVS on either build must not call a matched gene-level variant unresolved """
        grch37, grch38 = GenomeBuild.grch37(), GenomeBuild.grch38()
        allele_info = self._gene_level_allele_info(
            "ARV7", "GENE_LEVEL:644-644 <SPLICE:HGNC:644:V7>", grch37=self._resolved(grch37))
        display = allele_info.matched_without_resolved_hgvs_display(grch38)
        self.assertEqual("ARV7", display.full_hgvs)
        # the imported value is shown, so it carries the imported build rather than the one that matched
        self.assertEqual(grch38, display.genome_build)
        self.assertTrue(display.is_normalised)
        self.assertTrue(display.is_desired_build)
        self.assertTrue(display.is_resolved_without_hgvs)
        self.assertFalse(allele_info.matched_without_resolved_hgvs_display(grch37).is_desired_build)

        self.assertIsNone(self._allele_info(imported_c_hgvs=self.C_HGVS_38).matched_without_resolved_hgvs_display(grch38))
        with_c_hgvs = self._allele_info(imported_c_hgvs=self.C_HGVS_38,
                                        grch38=self._resolved(grch38, resolved_hgvs=self.C_HGVS_38))
        self.assertIsNone(with_c_hgvs.matched_without_resolved_hgvs_display(grch38))

    def test_unsupported_transcript_still_errors(self):
        allele_info = self._allele_info(
            imported_c_hgvs="NX_000059.4(BRCA2):c.1234A>G",
            grch37=self._resolved(GenomeBuild.grch37(), resolved_hgvs=self.C_HGVS_37, transcript_version_id=1),
            grch38=self._resolved(GenomeBuild.grch38(), resolved_hgvs=self.C_HGVS_38, transcript_version_id=2))
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
            grch38=self._resolved(GenomeBuild.grch38(), resolved_hgvs=self.C_HGVS_38, transcript_version_id=2))
        self.assertEqual(allele_info._calculate_validation()["builds"], {"missing_37": "W"})


class ImportedAlleleInfoGeneLevelDisabledTest(TestCase):
    """ With gene-level variants off, a value naming genes is just a c.HGVS that fails to parse """

    @override_settings(VARIANT_GENE_LEVEL_ENABLED=False)
    def test_fusion_string_fails_as_an_hgvs(self):
        allele_info = ImportedAlleleInfo.get_or_create(
            imported_c_hgvs="BCR::ABL1",
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                GenomeBuild.grch38()))
        self.assertIsNone(allele_info.variant_coordinate)
        self.assertEqual(ImportedAlleleInfoStatus.FAILED, allele_info.status)
        self.assertFalse(allele_info.is_gene_level)
        general = allele_info.latest_validation.validation_tags["general"]
        self.assertIn("cant_resolve_to_variant_coordinate", general)
        self.assertNotIn("gene_level_unresolved", general)


class ResolvedVariantInfoCNVTest(TestCase):
    """ #1574 - a <CNV> is a valid variant that HGVS simply cannot write, so it is recorded as an
        error on the ResolvedVariantInfo rather than reported as a bug """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.grch37()
        get_fake_annotation_version(cls.genome_build)

        contig = cls.genome_build.contigs.get(name="3")
        ref, _ = Sequence.objects.get_or_create(seq="N", seq_sha256_hash=sha256sum_str("N"))
        alt, _ = Sequence.objects.get_or_create(seq=VCFSymbolicAllele.CNV,
                                                seq_sha256_hash=sha256sum_str(VCFSymbolicAllele.CNV))
        locus, _ = Locus.objects.get_or_create(contig=contig, position=128200000, ref=ref)
        cls.variant, _ = Variant.objects.get_or_create(locus=locus, alt=alt, svlen=1000,
                                                       defaults={"end": 128201000})

        cls.allele_info = ImportedAlleleInfo.objects.create(
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                cls.genome_build),
            imported_c_hgvs="NM_001145661.2(GATA2):c.1018-213_1304del")

    @patch("classification.models.classification_variant_info_models.report_exc_info")
    def test_cnv_records_error_without_reporting(self, mock_report_exc_info):
        variant_info = ResolvedVariantInfo(genome_build=self.genome_build, allele_info=self.allele_info,
                                           variant=self.variant)
        variant_info.set_variant_and_save(self.variant)

        self.assertIsNone(variant_info.resolved_hgvs)
        self.assertIn("has no HGVS representation", variant_info.error)
        mock_report_exc_info.assert_not_called()
