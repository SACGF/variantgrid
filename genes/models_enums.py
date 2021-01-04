from django.db import models


class AnnotationConsortium(models.TextChoices):
    REFSEQ = "R", "RefSeq"
    ENSEMBL = "E", "Ensembl"


class HGNCStatus(models.TextChoices):
    APPROVED = 'A', 'Approved'
    SYMBOL_WITHDRAWN = 'S', 'Symbol Withdrawn'
    ENTRY_WITHDRAWN = 'E', 'Entry Withdrawn'


class GeneSymbolAliasSource(models.TextChoices):
    NCBI = "N", 'NCBI'
    UCSC = "U", 'UCSC'
    HGNC = "H", 'HGNC'
    MANUAL = "M", 'Manual'
