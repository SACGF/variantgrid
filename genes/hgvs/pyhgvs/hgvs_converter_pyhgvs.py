import logging
import sys
from importlib import metadata
from typing import Optional

import pyhgvs
from django.conf import settings
from pyhgvs import get_genomic_sequence, HGVSName
from pyhgvs.utils import make_transcript

from genes.hgvs import HGVSVariant, HGVSException
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele, HGVSConverterType
from genes.transcripts_utils import transcript_is_lrg
from snpdb.models import VariantCoordinate


class PyHGVSVariant(HGVSVariant):
    def __init__(self, hgvs_name: Optional[HGVSName] = None):
        if isinstance(hgvs_name, HGVSVariant):
            raise ValueError("Double extra!")
        self._hgvs_name = hgvs_name

    def _get_length(self) -> int:
        return self._hgvs_name.end + 1 - self._hgvs_name.start

    def _get_contig_accession(self) -> str:
        return self._hgvs_name.chrom

    def _get_gene(self) -> str:
        return self._hgvs_name.gene

    def _set_gene(self, value):
        self._hgvs_name.gene = value

    def _get_transcript(self) -> str:
        return self._hgvs_name.transcript

    def _set_transcript(self, value):
        self._hgvs_name.transcript = value

    def _get_kind(self) -> str:
        return self._hgvs_name.kind

    def _set_kind(self, value):
        self._hgvs_name.kind = value

    def _get_mutation_type(self):
        return self._hgvs_name.mutation_type

    def _get_ref_allele(self):
        return self._hgvs_name.ref_allele

    def _set_ref_allele(self, value):
        self._hgvs_name.ref_allele = value

    def _safe(self) -> HGVSName:
        params = vars(self._hgvs_name)
        params.pop('name', None)  # don't provide name a second time as parsing of name redundantly recalculates values
        copy = HGVSName(**params)
        return copy

    def get_ref_alt(self):
        # We don't ever use the kwargs so don't implement it
        return self._hgvs_name.get_ref_alt()

    def get_cdna_coords(self) -> str:
        return self._hgvs_name.format_cdna_coords()

    def format(self, use_compat=False, max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH) -> Optional[str]:
        # would be better practise to throw an error if we couldn't generate
        # but this keeps existing behaviour
        if not self._hgvs_name:
            return None

        if self.can_shrink_long_ref(self._hgvs_name, max_ref_length=max_ref_length):
            hgvs_name = self._safe()
            self.format_hgvs_remove_long_ref(hgvs_name, max_ref_length=max_ref_length)
            return hgvs_name.format()
        return self._hgvs_name.format()

    @staticmethod
    def can_shrink_long_ref(hgvs_name, max_ref_length=10) -> bool:
        SHRINKABLE_MUTATION_TYPES = {"del", "dup", "delins"}
        return hgvs_name.mutation_type in SHRINKABLE_MUTATION_TYPES and len(hgvs_name.ref_allele) > max_ref_length

    @staticmethod
    def format_hgvs_remove_long_ref(hgvs_name, max_ref_length=10):
        """ Similar to pyhgvs.variant_to_hgvs_name but only for dels, delins and dups and we don't specify length

            From a Facebook post:
            Q: What is the correct way to describe a deletion, c.7432-2025_7536+372del2502 or c.7432-2025_7536+372del.
            While ClinVar seems to prefer the first, #HGVS seems to prefer the second format.
            A: HGVS descriptions do not contain redundant information. The size of the deletion, in the example 2502
            nucleotides, can be deduced from the variant description. HGVS thus suggests to use c.7432-2025_7536+372del.
        """

        if PyHGVSVariant.can_shrink_long_ref(hgvs_name, max_ref_length=max_ref_length):
            hgvs_name.ref_allele = ""

    def ref_lengths(self) -> int:
        if not self._hgvs_name:
            return 0

        count = 0
        h_name = self._hgvs_name
        parts = [
            h_name.ref_allele,
            h_name.ref2_allele
        ]
        for part in parts:
            if part:
                count += len(part)
        return count

    def get_gene_symbol_if_no_transcript(self) -> Optional[str]:
        # PyHGVS HGVSName works like:
        # NM_001145661.2:c.1113dup          transcript=NM_001145661.2, gene=''
        # NM_001145661.2(GATA2):c.1113dup   transcript=NM_001145661.2, gene=GATA2
        # GATA2:c.1113dup                   transcript='',             gene=GATA2
        gene_symbol = None
        if not self._get_transcript():
            gene_symbol = self.gene
        return gene_symbol

    def __str__(self):
        return self.format()


class PyHGVSConverter(HGVSConverter):
    @staticmethod
    def _hgvs_name(hgvs_string):
        """ Catches PyHGVS specific exceptions and converts to HGVSException """

        HGVSConverter._hgvs_string_validation(hgvs_string)
        try:
            return HGVSName(hgvs_string)
        except pyhgvs.InvalidHGVSName as e:
            raise HGVSException(str(e)) from e

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        return PyHGVSVariant(self._hgvs_name(hgvs_string))

    def normalize(self, hgvs_variant: PyHGVSVariant) -> HGVSVariant:
        logging.warning("PyHGVSConverter normalize is a no-op")
        return hgvs_variant

    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt,
                                                self.genome_build.genome_fasta.fasta,
                                                transcript=None, max_allele_length=sys.maxsize)
        contig = self.genome_build.chrom_contig_mappings[chrom]
        hgvs_name.chrom = contig.refseq_accession
        return PyHGVSVariant(hgvs_name)

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        pyhgvs_transcript = make_transcript(transcript_version.pyhgvs_data)
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt, self.genome_build.genome_fasta.fasta,
                                                pyhgvs_transcript, max_allele_length=sys.maxsize)
        return PyHGVSVariant(hgvs_name)

    def hgvs_to_variant_coordinate_reference_match_and_normalized(self, hgvs_string: str, transcript_version) -> tuple[VariantCoordinate, HgvsMatchRefAllele, bool]:
        pyhgvs_transcript = None
        hgvs_name = self._hgvs_name(hgvs_string)

        # Check transcript bounds
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
        """ Only returns anything for c.HGVS """
        if hgvs_string is not None:
            hgvs_name = self._hgvs_name(hgvs_string)
            transcript_accession = hgvs_name.transcript
        else:
            transcript_accession = ''
        return transcript_accession

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.PYHGVS

    def get_version(self) -> str:
        return metadata.version('pyhgvs')

    def get_hgvs_match_ref_allele(self, hgvs_name, pyhgvs_transcript=None) -> HgvsMatchRefAllele:
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
                raise HGVSException(f"{hgvs_name}: {description} {cdna_coord} resolves outside of transcript")
