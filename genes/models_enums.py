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
        raise ValueError(f"Couldn't determine annotation consortium for \"{transcript_accession}\"")


class HGNCStatus(models.TextChoices):
    APPROVED = 'A', 'Approved'
    SYMBOL_WITHDRAWN = 'S', 'Symbol Withdrawn'
    ENTRY_WITHDRAWN = 'E', 'Entry Withdrawn'


class HGVSKind(models.TextChoices):
    CODING = "c", "coding"
    GENOMIC = "g", "genomic"
    MITOCHONDRIA = "m", "mitochondria"
    NON_CODING = "n", "non-coding"
    CIRCULAR_GENOMIC = "o", "circular genomic"
    PROTEIN = "p", "protein"
    RNA = "r", "RNA transcript"


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

    def get_css_class(self):
        CSS_CLASSES = {
            self.HIGH: "HighEvidence",
            self.INTERMEDIATE: "ModerateEvidence",
            self.LOW: "LowEvidence",
        }
        return CSS_CLASSES[self.value]
