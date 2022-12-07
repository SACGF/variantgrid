from django.db import models


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
        raise ValueError(f"Couldn't determine annotation consortium for {transcript_accession}")


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
