import sys
from importlib import metadata
from typing import Tuple

import pyhgvs
from pyhgvs import get_genomic_sequence, HGVSName
from pyhgvs.utils import make_transcript

from genes.hgvs import HGVSNameExtra
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele
from genes.models import TranscriptVersion
from genes.refseq_transcripts import transcript_is_lrg
from snpdb.models import GenomeBuild, VariantCoordinate


class PyHGVSConverter(HGVSConverter):
    def __int__(self, genome_build: GenomeBuild):
        super().__init__(genome_build)

    def variant_coords_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSNameExtra:
        chrom, offset, ref, alt = vc
        hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt,
                                                self.genome_build.genome_fasta.fasta,
                                                transcript=None, max_allele_length=sys.maxsize)
        contig = self.genome_build.chrom_contig_mappings[chrom]
        hgvs_name.chrom = contig.refseq_accession
        return HGVSNameExtra(hgvs_name)

    def variant_coords_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSNameExtra:
        pyhgvs_transcript = make_transcript(transcript_version.pyhgvs_data)
        hgvs_name = pyhgvs.variant_to_hgvs_name(*vc, self.genome_build.genome_fasta.fasta,
                                                pyhgvs_transcript, max_allele_length=sys.maxsize)
        return HGVSNameExtra(hgvs_name)

    def hgvs_to_variant_coords_and_reference_match(self, hgvs_string: str, transcript_version) -> Tuple[VariantCoordinate, HgvsMatchRefAllele]:
        pyhgvs_transcript = None
        hgvs_name = HGVSName(hgvs_string)

        # Check transcript bounds
        if transcript_version:
            pyhgvs_transcript = make_transcript(transcript_version.pyhgvs_data)
            self._validate_in_transcript_range(pyhgvs_transcript, hgvs_name)

        variant_tuple = pyhgvs.parse_hgvs_name(hgvs_string, self.genome_build.genome_fasta.fasta,
                                               transcript=pyhgvs_transcript,
                                               indels_start_with_same_base=False)

        chrom, position, ref, alt = variant_tuple
        contig = self.genome_build.chrom_contig_mappings[chrom]
        chrom = contig.name

        matches_reference = self.get_hgvs_match_ref_allele(hgvs_name, pyhgvs_transcript)
        return (chrom, position, ref, alt), matches_reference

    def hgvs_clean_for_clingen(self, hgvs_string: str) -> str:
        # ClinGen Allele Registry doesn't like gene names - so strip (unless LRG_)
        hgvs_name = HGVSName(hgvs_string)
        transcript_accession = self.get_transcript_accession(hgvs_string)
        if not transcript_is_lrg(transcript_accession):
            hgvs_name.gene = None
        return hgvs_name.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        hgvs_name = HGVSName(hgvs_string)
        return hgvs_name.transcript

    def description(self) -> str:
        return f"pyhgvs v{metadata.version('pyhgvs')}"

    def get_hgvs_match_ref_allele(self, hgvs_name, pyhgvs_transcript = None) -> HgvsMatchRefAllele:
        """Return True if reference allele matches genomic sequence."""
        is_forward_strand = pyhgvs_transcript.tx_position.is_forward_strand if pyhgvs_transcript else True
        ref, _ = hgvs_name.get_ref_alt(is_forward_strand,
                                       raw_dup_alleles=True)  # get raw values so dup isn't always True
        chrom, start, end = hgvs_name.get_raw_coords(pyhgvs_transcript)
        fasta = self.genome_build.genome_fasta.fasta
        genome_ref = get_genomic_sequence(fasta, chrom, start, end)
        return HgvsMatchRefAllele(provided_ref=ref, calculated_ref=genome_ref)

    @staticmethod
    def _validate_in_transcript_range(pyhgvs_transcript, hgvs_name: HGVSName):
        # @see https://varnomen.hgvs.org/bg-material/numbering/
        # Transcript Flanking: it is not allowed to describe variants in nucleotides beyond the boundaries of the
        # reference sequence using the reference sequence
        NAMES = {
            "cdna_start": hgvs_name.cdna_start,
            "cdna_end": hgvs_name.cdna_end,
        }
        for description, cdna_coord in NAMES.items():
            genomic_coord = pyhgvs_transcript.cdna_to_genomic_coord(cdna_coord)
            within_transcript = pyhgvs_transcript.tx_position.chrom_start <= genomic_coord <= pyhgvs_transcript.tx_position.chrom_stop
            if not within_transcript:
                raise pyhgvs.InvalidHGVSName(f"'{hgvs_name.format()}' {description} {cdna_coord} resolves outside of transcript")
