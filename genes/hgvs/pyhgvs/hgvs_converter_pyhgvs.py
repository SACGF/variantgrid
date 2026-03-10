import sys
from typing import Optional

import pyhgvs
from pyhgvs import get_genomic_sequence, HGVSName
from pyhgvs.utils import make_transcript

from hgvs_shim import PyHGVSConverter as _PyHGVSConverter, TranscriptInfo
from hgvs_shim.hgvs_converter_pyhgvs import PyHGVSVariant
from genes.hgvs import HGVSException
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele, HGVSConverterType
from genes.transcripts_utils import transcript_is_lrg
from snpdb.models import VariantCoordinate


class PyHGVSConverter(HGVSConverter):
    """Django-coupled adapter for hgvs_shim.PyHGVSConverter"""

    def __init__(self, genome_build, local_resolution=True, clingen_resolution=True):
        super().__init__(genome_build, local_resolution=local_resolution,
                         clingen_resolution=clingen_resolution)
        fasta = genome_build.genome_fasta.fasta
        self._shim = _PyHGVSConverter(fasta_file=fasta, get_transcript=self._load_pyhgvs_transcript)

    def _load_pyhgvs_transcript(self, accession: str):
        """Look up pyhgvs transcript data from the Django DB"""
        from genes.models import TranscriptVersion
        tv = TranscriptVersion.get_transcript_version(self.genome_build, accession)
        return make_transcript(tv.pyhgvs_data)

    @staticmethod
    def _hgvs_name(hgvs_string) -> HGVSName:
        """Catches PyHGVS specific exceptions and converts to HGVSException"""
        HGVSConverter._hgvs_string_validation(hgvs_string)
        try:
            return HGVSName(hgvs_string)
        except pyhgvs.InvalidHGVSName as e:
            raise HGVSException(str(e)) from e

    def create_hgvs_variant(self, hgvs_string: str) -> PyHGVSVariant:
        return self._shim.create_hgvs_variant(hgvs_string)

    def normalize(self, hgvs_variant: PyHGVSVariant) -> PyHGVSVariant:
        return self._shim.normalize(hgvs_variant)

    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> PyHGVSVariant:
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt,
                                                self.genome_build.genome_fasta.fasta,
                                                transcript=None, max_allele_length=sys.maxsize)
        contig = self.genome_build.chrom_contig_mappings[chrom]
        hgvs_name.chrom = contig.refseq_accession
        return PyHGVSVariant(hgvs_name)

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> PyHGVSVariant:
        pyhgvs_transcript = make_transcript(transcript_version.pyhgvs_data)
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt, self.genome_build.genome_fasta.fasta,
                                                pyhgvs_transcript, max_allele_length=sys.maxsize)
        return PyHGVSVariant(hgvs_name)

    def hgvs_to_variant_coordinate_reference_match_and_normalized(
            self, hgvs_string: str, transcript_version=None
    ) -> tuple[VariantCoordinate, HgvsMatchRefAllele, bool]:
        pyhgvs_transcript = None
        hgvs_name = self._hgvs_name(hgvs_string)

        if transcript_version:
            pyhgvs_transcript = make_transcript(transcript_version.pyhgvs_data)
            self._validate_in_transcript_range(pyhgvs_transcript, hgvs_name)

        chrom, position, ref, alt = pyhgvs.parse_hgvs_name(hgvs_string, self.genome_build.genome_fasta.fasta,
                                                        transcript=pyhgvs_transcript,
                                                        indels_start_with_same_base=False)
        contig = self.genome_build.chrom_contig_mappings[chrom]
        chrom = contig.name
        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref, alt)
        matches_reference = self.get_hgvs_match_ref_allele(hgvs_name, pyhgvs_transcript)
        return vc, matches_reference, True

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        # ClinGen Allele Registry doesn't like gene names - so strip (unless LRG_)
        hgvs_name = self._hgvs_name(hgvs_string)
        transcript_accession = self.get_transcript_accession(hgvs_string)
        if not transcript_is_lrg(transcript_accession):
            hgvs_name.gene = None
        return hgvs_name.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        return self._shim.get_transcript_accession(hgvs_string)

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.PYHGVS

    def get_version(self) -> str:
        return self._shim.get_version()

    # --- VG-specific methods ---

    def get_hgvs_match_ref_allele(self, hgvs_name: HGVSName, pyhgvs_transcript=None) -> HgvsMatchRefAllele:
        """Return reference match info by comparing provided ref against genomic sequence."""
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
                raise HGVSException(f"{hgvs_name}: {description} {cdna_coord} resolves outside of transcript")
