"""
ResolvedVariantInfo's canonical p.HGVS: source order (VEP, computed, ClinGen, evidence), the guards that reject a
bad value from a source, and the set_variant_and_save hook (#1383)
"""
from io import StringIO
from unittest.mock import patch

from django.core.management import call_command
from django.test import TestCase
from hgvs.exceptions import HGVSDataNotAvailableError
from hgvs.parsers import Parser

from annotation.fake_data import create_fake_variant_annotation, get_fake_annotation_version
from annotation.models import VariantTranscriptAnnotation
from classification.enums import PHGVSSource, SpecialEKeys, SubmissionSource
from classification.models import Classification, ImportedAlleleInfo
from classification.models.classification_variant_info_models import (
    HGVSResolution,
    ResolvedVariantInfo,
)
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.fake_data import create_gata2_transcript_version
from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.biocommons_hgvs.hgvs_converter_biocommons import BioCommonsHGVSConverter
from genes.models import TranscriptVersionSequenceInfo
from snpdb.models import Allele, ClinGenAllele, GenomeBuild, GenomeBuildPatchVersion
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

C_HGVS = "NM_001145661.2(GATA2):c.4G>A"
COMPUTED_P_HGVS = "NP_001139133.1:p.(Glu2Lys)"  # GAG -> AAG, on the fixture's cdot protein


class ResolvedVariantInfoPHGVSTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.transcript_version = create_gata2_transcript_version(cls.genome_build)
        # cache the sequence (from the mocked fetcher) so the computed path never needs to retrieve it
        TranscriptVersionSequenceInfo.get(cls.transcript_version.accession, retrieve=True)
        cls.variant = slowly_create_test_variant("3", 128205871, "C", "T", cls.genome_build)

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.allele = Allele.objects.create()
        self.allele_info = ImportedAlleleInfo.objects.create(
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                self.genome_build),
            imported_c_hgvs=C_HGVS, allele=self.allele)
        self.rvi = ResolvedVariantInfo.objects.create(allele_info=self.allele_info, genome_build=self.genome_build,
                                                      variant=self.variant, resolved_hgvs=C_HGVS,
                                                      transcript_version=self.transcript_version)

    def _vep(self, hgvs_p: str):
        vav = self.annotation_version.variant_annotation_version
        va = create_fake_variant_annotation(self.variant, vav)
        VariantTranscriptAnnotation.objects.create(version=vav, variant=self.variant, annotation_run=va.annotation_run,
                                                   transcript_version=self.transcript_version, hgvs_p=hgvs_p)

    def _clingen(self, p_hgvs: str):
        api_response = {"@id": "http://reg.genome.network/allele/CA1383",
                        "transcriptAlleles": [{"hgvs": ["NM_001145661.2:c.4G>A"],
                                               "proteinEffect": {"hgvs": p_hgvs}}]}
        self.allele.clingen_allele = ClinGenAllele.objects.create(id=1383, api_response=api_response)
        self.allele.save()

    def _evidence(self, p_hgvs: str):
        classification = Classification.create(user=self.user, lab=self.lab, save=True, source=SubmissionSource.API,
                                               data={SpecialEKeys.P_HGVS: {"value": p_hgvs}})
        classification.publish_latest(user=self.user)
        Classification.objects.filter(pk=classification.pk).update(allele_info=self.allele_info)

    def _computed_fails(self):
        return patch.object(BioCommonsHGVSConverter, "c_hgvs_to_p_sequence_variant",
                            side_effect=NotImplementedError("sequence length of 1378 is not divisible by 3"))

    def test_vep_wins_over_computed(self):
        self._vep("NP_001139133.1:p.Glu2Ter")
        resolution = self.rvi.recalc_p_hgvs()
        self.assertEqual((resolution.p_hgvs, resolution.source), ("NP_001139133.1:p.(Glu2Ter)", PHGVSSource.VEP))

    def test_computed_without_vep_for_exact_transcript_version(self):
        resolution = self.rvi.recalc_p_hgvs()
        self.assertEqual((resolution.p_hgvs, resolution.source), (COMPUTED_P_HGVS, PHGVSSource.COMPUTED))

    def test_clingen_used_when_computed_fails(self):
        self._clingen("NP_001139133.1:p.Glu2Lys")
        with self._computed_fails():
            resolution = self.rvi.recalc_p_hgvs()
        self.assertEqual((resolution.p_hgvs, resolution.source), ("NP_001139133.1:p.(Glu2Lys)", PHGVSSource.CLINGEN))

    def test_clingen_non_protein_rejected(self):
        """ ClinGen returns n. on the protein accession for intronic variants """
        self._clingen("NP_001139133.1:n.123+1G>A")
        with self._computed_fails():
            resolution = self.rvi.recalc_p_hgvs()
        self.assertIsNone(resolution.p_hgvs)
        self.assertIn("not divisible by 3", resolution.error)

    def test_computed_without_protein_accession_is_an_error(self):
        """ Pre-cdot RefSeq rows have no "protein" - biocommons would otherwise fall through to an unimplemented
            lookup. The other sources are still tried """
        self._clingen("NP_001139133.1:p.Glu2Lys")
        # biocommons wraps the provider's methods in its own lru_cache at init, which needs a real function
        with patch.object(DjangoTranscriptDataProvider, "get_pro_ac_for_tx_ac", lambda _self, _tx_ac: None):
            converter = BioCommonsHGVSConverter(self.genome_build)
            with self.assertRaises(HGVSDataNotAvailableError):
                converter.c_hgvs_to_p_sequence_variant(C_HGVS)
            resolution = self.rvi.recalc_p_hgvs(converter=converter)
        self.assertEqual((resolution.p_hgvs, resolution.source), ("NP_001139133.1:p.(Glu2Lys)", PHGVSSource.CLINGEN))

    def test_evidence_typo_skipped(self):
        """ The lab's accession is stored as typed - here none """
        self._evidence("p.Glu2Lys")
        self._evidence("p.Glx2Lsy")  # latest, but doesn't parse
        with self._computed_fails():
            resolution = self.rvi.recalc_p_hgvs()
        self.assertEqual((resolution.p_hgvs, resolution.source), ("p.(Glu2Lys)", PHGVSSource.EVIDENCE))

    def test_frameshift_in_disguise_falls_through(self):
        """ A delins running to the last residue is a frameshift biocommons couldn't name - and ClinGen's
            unknown stop is no better than that """
        disguised = Parser().parse("NP_001139133.1:p.Glu2_Gly480delinsIleLeuTer")
        self._clingen("NP_001139133.1:p.Glu2IlefsTer?")
        with patch.object(BioCommonsHGVSConverter, "c_hgvs_to_p_sequence_variant", return_value=disguised):
            resolution = self.rvi.recalc_p_hgvs()
        self.assertIsNone(resolution.p_hgvs)
        self.assertIn("frameshift", resolution.error)

    def test_set_variant_and_save_computed_failure_keeps_c_hgvs(self):
        hgvs_resolution = HGVSResolution(resolved_hgvs=C_HGVS, hgvs_converter_version=None,
                                         hgvs_converter_data_version="", transcript_version=self.transcript_version,
                                         gene_symbol=None)
        with patch.object(ResolvedVariantInfo, "recalc_resolved_hgvs", return_value=hgvs_resolution), \
                self._computed_fails():
            self.rvi.set_variant_and_save(self.variant)
        self.rvi.refresh_from_db()
        self.assertEqual(C_HGVS, self.rvi.resolved_hgvs)
        self.assertIsNone(self.rvi.p_hgvs)
        self.assertIn("not divisible by 3", self.rvi.p_hgvs_error)

    def test_set_variant_and_save_without_transcript_leaves_p_hgvs_null(self):
        self.rvi.p_hgvs = COMPUTED_P_HGVS
        self.rvi.p_hgvs_source = PHGVSSource.COMPUTED
        g_hgvs_resolution = HGVSResolution(resolved_hgvs="NC_000003.11:g.128205871C>T", hgvs_converter_version=None,
                                           hgvs_converter_data_version="", transcript_version=None, gene_symbol=None)
        with patch.object(ResolvedVariantInfo, "recalc_resolved_hgvs", return_value=g_hgvs_resolution):
            self.rvi.set_variant_and_save(self.variant)
        self.rvi.refresh_from_db()
        self.assertEqual((None, None, None), (self.rvi.p_hgvs, self.rvi.p_hgvs_source, self.rvi.p_hgvs_error))

    def test_backfill_command(self):
        call_command("one_off_backfill_resolved_variant_info_p_hgvs", "--genome-build", "GRCh37", "--sync",
                     stdout=StringIO())
        self.rvi.refresh_from_db()
        self.assertEqual((COMPUTED_P_HGVS, PHGVSSource.COMPUTED), (self.rvi.p_hgvs, self.rvi.p_hgvs_source))
