"""
VariantTranscriptSelections builds the transcript table shown on the variant and "create
classification" pages. These cover how it treats symbolic alts when adding the other annotation
consortium's transcripts: DEL/DUP/INV go to the converter as coordinates (#1571) and <CNV>/<INS>
have no HGVS at all (#1574), so neither is expanded to explicit ref/alt.
"""
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import AnnotationRun
from annotation.models.models import VariantAnnotation, VariantTranscriptAnnotation
from annotation.tests.test_data_fake_genes import create_gata2_transcript_version
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils import sha256sum_str
from snpdb.models import GenomeBuild, Locus, Sequence, Variant


class VariantTranscriptSelectionsSymbolicTest(TestCase):
    GATA2_TRANSCRIPT = "NM_001145661.2"
    POSITION = 128200000
    SVLEN = 1000  # Under HGVS_MAX_SEQUENCE_LENGTH, so the length guard doesn't skip the other consortium

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.grch37()
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.vav = cls.annotation_version.variant_annotation_version  # Ensembl, so RefSeq is "other"
        create_gata2_transcript_version(cls.genome_build)
        cls.annotation_run = AnnotationRun.objects.create()

    @classmethod
    def _create_symbolic_variant(cls, alt: str, svlen: int) -> Variant:
        contig = cls.genome_build.chrom_contig_mappings["3"]
        ref_seq, _ = Sequence.objects.get_or_create(seq="N", seq_sha256_hash=sha256sum_str("N"))
        alt_seq, _ = Sequence.objects.get_or_create(seq=alt, seq_sha256_hash=sha256sum_str(alt))
        locus, _ = Locus.objects.get_or_create(contig=contig, position=cls.POSITION, ref=ref_seq)
        variant, _ = Variant.objects.get_or_create(locus=locus, alt=alt_seq, svlen=svlen,
                                                   defaults={"end": cls.POSITION + abs(svlen)})
        # VEP put it in GATA2 without a transcript - the symbol is what brings RefSeq transcripts in
        VariantAnnotation.objects.create(version=cls.vav, variant=variant, annotation_run=cls.annotation_run)
        VariantTranscriptAnnotation.objects.create(version=cls.vav, variant=variant,
                                                   annotation_run=cls.annotation_run, symbol="GATA2")
        return variant

    def _gata2_transcript_data(self, alt: str, svlen: int) -> dict:
        variant = self._create_symbolic_variant(alt, svlen)
        vts = VariantTranscriptSelections(variant, self.genome_build, self.annotation_version)
        by_refseq = {td.get(VariantTranscriptSelections.REFSEQ_TRANSCRIPT): td for td in vts.transcript_data}
        self.assertIn(self.GATA2_TRANSCRIPT, by_refseq, "RefSeq transcript added for the Ensembl-annotated variant")
        return by_refseq[self.GATA2_TRANSCRIPT]

    def test_cnv_transcript_row_without_hgvs(self):
        """ #1574 - a <CNV> shorter than HGVS_MAX_SEQUENCE_LENGTH used to raise
            "Unknown symbolic alt of '<CNV>'" out of as_external_explicit() """
        t_data = self._gata2_transcript_data(VCFSymbolicAllele.CNV, self.SVLEN)
        self.assertIsNone(t_data.get("hgvs_c"))

    def test_del_transcript_row_keeps_hgvs(self):
        t_data = self._gata2_transcript_data(VCFSymbolicAllele.DEL, -self.SVLEN)
        self.assertEqual("NM_001145661.2(GATA2):c.1018-213_1304del", t_data["hgvs_c"])
