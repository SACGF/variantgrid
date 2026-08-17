from django.test import SimpleTestCase, TestCase

from annotation.models.models_enums import NMDEscapeStatus
from annotation.ptc import calculate_ptc_position, parse_ptc_distance_codons, predict_nmd_escape
from annotation.tests.test_data_fake_genes import create_gata2_transcript_version
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import (
    TranscriptGeometryCache,
    add_calculated_ptc,
)
from snpdb.models.models_genome import GenomeBuild


class TestParsePTCDistanceCodons(SimpleTestCase):
    """ NM_000251.3(MSH2):c.2527delT -> NP_000242.1:p.Cys843ValfsTer49 (#579) """

    def test_frameshift_ter(self):
        self.assertEqual(parse_ptc_distance_codons("NP_000242.1:p.Cys843ValfsTer49"), 49)

    def test_frameshift_ter_no_accession(self):
        self.assertEqual(parse_ptc_distance_codons("p.Cys843ValfsTer49"), 49)

    def test_changed_codon_is_the_stop(self):
        """ p.Tyr535Ter - the AR variant from #665. HGVS omits 'fs' when the changed codon
            is itself a stop, so the PTC is the changed codon. """
        self.assertEqual(parse_ptc_distance_codons("NP_000035.2:p.Tyr535Ter"), 1)

    def test_unlocatable_stop(self):
        self.assertIsNone(parse_ptc_distance_codons("NP_000242.1:p.Cys843ValfsTer?"))

    def test_no_stop_in_hgvs_p(self):
        self.assertIsNone(parse_ptc_distance_codons("NP_000242.1:p.Cys843Val"))

    def test_empty(self):
        self.assertIsNone(parse_ptc_distance_codons(""))
        self.assertIsNone(parse_ptc_distance_codons(None))


class TestCalculatePTCPosition(SimpleTestCase):

    def test_msh2_deletion(self):
        """ c.2527delT / p.Cys843ValfsTer49: PTC protein position 843 + 49 - 1 = 891,
            CDS 3 * 890 + 1 = 2671 in the mutant, + 1 deleted base = 2672 in the reference. """
        self.assertEqual(calculate_ptc_position(843, 49, indel_offset=1), 2672)

    def test_insertion_offset_shifts_the_other_way(self):
        """ A 2 base insertion (len(ref) - len(alt) = -2) puts the reference position earlier. """
        self.assertEqual(calculate_ptc_position(843, 49, indel_offset=-2), 2669)

    def test_changed_codon_is_the_stop(self):
        """ p.Tyr535Ter -> PTC at protein position 535, CDS 3 * 534 + 1 = 1603. """
        self.assertEqual(calculate_ptc_position(535, 1, indel_offset=0), 1603)


class TestPredictNMDEscape(SimpleTestCase):
    """ The four NMD.pm rules, anchored on the PTC. """

    def test_single_exon_transcript_escapes(self):
        status = predict_nmd_escape(ptc_cds=2672, ptc_last_junction_distance=500, single_exon=True)
        self.assertEqual(status, NMDEscapeStatus.ESCAPING)

    def test_ptc_in_final_exon_escapes(self):
        status = predict_nmd_escape(ptc_cds=2672, ptc_last_junction_distance=-120, single_exon=False)
        self.assertEqual(status, NMDEscapeStatus.ESCAPING)

    def test_last_junction_distance_50_escapes(self):
        status = predict_nmd_escape(ptc_cds=2672, ptc_last_junction_distance=50, single_exon=False)
        self.assertEqual(status, NMDEscapeStatus.ESCAPING)

    def test_last_junction_distance_51_predicts_nmd(self):
        status = predict_nmd_escape(ptc_cds=2672, ptc_last_junction_distance=51, single_exon=False)
        self.assertEqual(status, NMDEscapeStatus.PREDICTED_NMD)

    def test_start_proximal_ptc_cds_100_escapes(self):
        status = predict_nmd_escape(ptc_cds=100, ptc_last_junction_distance=5000, single_exon=False)
        self.assertEqual(status, NMDEscapeStatus.ESCAPING)

    def test_start_proximal_ptc_cds_101_predicts_nmd(self):
        status = predict_nmd_escape(ptc_cds=101, ptc_last_junction_distance=5000, single_exon=False)
        self.assertEqual(status, NMDEscapeStatus.PREDICTED_NMD)


class TestAddCalculatedPTC(TestCase):
    """ The maths combined with real transcript geometry - GATA2 NM_001145661.2 on GRCh37 has
        5'UTR 435, CDS 1443 and its last exon-exon junction at cDNA 1578. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.transcript_version = create_gata2_transcript_version(grch37)

    def _add_ptc(self, indel_offset=1, **transcript_data) -> dict:
        transcript_data.setdefault("transcript_version_id", self.transcript_version.pk)
        add_calculated_ptc(transcript_data, indel_offset, TranscriptGeometryCache())
        return transcript_data

    def test_frameshift_upstream_of_last_junction_predicts_nmd(self):
        """ p.Val100GlyfsTer10 -> PTC at protein 109, CDS 325 (+1 deleted base = 326),
            cDNA 761 - 817nt short of the last junction. """
        transcript_data = self._add_ptc(consequence="frameshift_variant",
                                        hgvs_p="NP_001139133.1:p.Val100GlyfsTer10",
                                        protein_position="100")
        self.assertEqual(transcript_data["ptc_distance_codons"], 10)
        self.assertEqual(transcript_data["ptc_last_junction_distance"], 817)
        self.assertEqual(transcript_data["nmd_escape_status"], NMDEscapeStatus.PREDICTED_NMD)

    def test_frameshift_in_final_exon_escapes(self):
        """ p.Val400GlyfsTer5 -> PTC cDNA 1646, past the last junction at 1578. """
        transcript_data = self._add_ptc(consequence="frameshift_variant",
                                        hgvs_p="NP_001139133.1:p.Val400GlyfsTer5",
                                        protein_position="400")
        self.assertEqual(transcript_data["ptc_last_junction_distance"], -68)
        self.assertEqual(transcript_data["nmd_escape_status"], NMDEscapeStatus.ESCAPING)

    def test_protein_position_range_uses_first_int(self):
        """ VEP writes indel protein positions as eg '100-101'. """
        transcript_data = self._add_ptc(consequence="frameshift_variant",
                                        hgvs_p="NP_001139133.1:p.Val100GlyfsTer10",
                                        protein_position="100-101")
        self.assertEqual(transcript_data["ptc_last_junction_distance"], 817)

    def test_non_frameshift_is_not_applicable(self):
        transcript_data = self._add_ptc(consequence="missense_variant",
                                        hgvs_p="NP_001139133.1:p.Val100Gly",
                                        protein_position="100")
        self.assertEqual(transcript_data["nmd_escape_status"], NMDEscapeStatus.NOT_APPLICABLE)
        self.assertNotIn("ptc_distance_codons", transcript_data)

    def test_unlocatable_stop_is_not_applicable(self):
        transcript_data = self._add_ptc(consequence="frameshift_variant",
                                        hgvs_p="NP_001139133.1:p.Val100GlyfsTer?",
                                        protein_position="100")
        self.assertEqual(transcript_data["nmd_escape_status"], NMDEscapeStatus.NOT_APPLICABLE)
        self.assertNotIn("ptc_distance_codons", transcript_data)

    def test_no_transcript_version_is_not_applicable(self):
        """ Non-coding / unlinked transcripts have no geometry to place the PTC against. """
        transcript_data = self._add_ptc(consequence="frameshift_variant",
                                        hgvs_p="NP_001139133.1:p.Val100GlyfsTer10",
                                        protein_position="100",
                                        transcript_version_id=None)
        self.assertEqual(transcript_data["nmd_escape_status"], NMDEscapeStatus.NOT_APPLICABLE)

    def test_geometry_read_once_per_transcript(self):
        """ The cdot blob is read once per transcript per run, not once per frameshift row. """
        cache = TranscriptGeometryCache()

        def _add():
            transcript_data = {"consequence": "frameshift_variant",
                               "hgvs_p": "NP_001139133.1:p.Val100GlyfsTer10",
                               "protein_position": "100",
                               "transcript_version_id": self.transcript_version.pk}
            add_calculated_ptc(transcript_data, 1, cache)

        _add()
        with self.assertNumQueries(0):
            _add()
            _add()
