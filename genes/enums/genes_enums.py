"""
Enums for genes and transcripts: AnnotationConsortium (RefSeq / Ensembl - the axis every Gene and
Transcript id carries), HGNC status, alias sources, LRG RefSeqGene category, MANE status, PanelApp
confidence and GeneCopyNumberEventKind (gain / loss of a whole-gene copy number event).
Stored codes: add values, never renumber.
"""
from typing import Optional

from django.db import models

from library.genomics.vcf_enums import GeneLevelSymbolicAlt


class AnnotationConsortium(models.TextChoices):
    REFSEQ = "R", "RefSeq"
    ENSEMBL = "E", "Ensembl"

    @staticmethod
    def get_from_transcript_accession(transcript_accession: str):
        TRANSCRIPT_PREFIXES = [
            (AnnotationConsortium.ENSEMBL, 4, {"ENST"}),
            (AnnotationConsortium.REFSEQ, 3, {'NM_', 'NR_', 'NG_', 'XM_', 'XR_'}),
        ]

        for (ac, length, prefixes) in TRANSCRIPT_PREFIXES:
            prefix = transcript_accession[:length]
            if prefix in prefixes:
                return ac
        raise ValueError(f"Couldn't determine annotation consortium for \"{transcript_accession}\"")


class HGNCStatus(models.TextChoices):
    APPROVED = 'A', 'Approved'
    SYMBOL_WITHDRAWN = 'S', 'Symbol Withdrawn'
    ENTRY_WITHDRAWN = 'E', 'Entry Withdrawn'


class GeneSymbolAliasSource(models.TextChoices):
    NCBI = "N", 'NCBI'
    UCSC = "U", 'UCSC'
    HGNC = "H", 'HGNC'
    MANUAL = "M", 'Manual'


class LRGRefSeqGeneCategory(models.TextChoices):
    ALIGNED_SELECTED = "S", 'aligned: Selected'
    ALIGNED_HISTORICAL = "H", "aligned: historical"
    REFERENCE_STANDARD = "R", "reference standard"


class MANEStatus(models.TextChoices):
    """
    MANE Select = 1 transcript per gene
    MANE Plus Clinical - additional transcripts for genes where MANE Select alone is not sufficient to report all
                         Pathogenic or Likely Pathogenic clinical variants available in public resources.
    """
    MANE_SELECT = "M", 'MANE Select'
    MANE_PLUS_CLINICAL = "C", "MANE Plus Clinical"


class PanelAppConfidence(models.TextChoices):
    """
    https://panelapp.genomicsengland.co.uk/media/files/PanelApp_User_Guide.pdf
    * Red = lowest level of confidence; 1 of the 4 sources or from other sources
    * Amber = intermediate; a gene from 2 sources.
    * Green = highest level of confidence; a gene from 3 or 4 sources.
    """
    LOW = "1", 'Low'
    INTERMEDIATE = "2", "Intermediate"
    HIGH = "3", "High"


class GeneCopyNumberEventKind(models.TextChoices):
    """ Which way a whole-gene copy number call went. GAIN rather than amplification because the
        threshold that makes a gain an amplification is the lab's, not ours - identity never encodes
        one, and the copy ratio stays per-sample on CohortGenotype. """

    GAIN = "G", "Gain"
    LOSS = "L", "Loss"

    @property
    def alt_kind(self) -> str:
        """ The GeneLevelSymbolicAlt this kind's Variant alt is written with """
        return GENE_COPY_NUMBER_ALT_KINDS[self]

    @property
    def canonical_word(self) -> str:
        """ 'amplification' is what the reports say; 'gain' is what we store """
        return GENE_COPY_NUMBER_CANONICAL_WORDS[self]

    @staticmethod
    def from_alt_kind(alt_kind: str) -> Optional['GeneCopyNumberEventKind']:
        return GENE_COPY_NUMBER_KINDS_BY_ALT.get(alt_kind)


GENE_COPY_NUMBER_ALT_KINDS = {
    GeneCopyNumberEventKind.GAIN: GeneLevelSymbolicAlt.GAIN,
    GeneCopyNumberEventKind.LOSS: GeneLevelSymbolicAlt.LOSS,
}
GENE_COPY_NUMBER_KINDS_BY_ALT = {alt_kind: kind for kind, alt_kind in GENE_COPY_NUMBER_ALT_KINDS.items()}
GENE_COPY_NUMBER_CANONICAL_WORDS = {
    GeneCopyNumberEventKind.GAIN: "amplification",
    GeneCopyNumberEventKind.LOSS: "loss",
}
