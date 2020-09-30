class AnnotationConsortium:
    REFSEQ = "R"
    ENSEMBL = 'E'

    CHOICES = (
        (REFSEQ, "RefSeq"),
        (ENSEMBL, "Ensembl"),
    )


class HGNCStatus:
    APPROVED = 'A'
    SYMBOL_WITHDRAWN = 'S'
    ENTRY_WITHDRAWN = 'E'

    CHOICES = (
        (APPROVED, 'Approved'),
        (SYMBOL_WITHDRAWN, 'Symbol Withdrawn'),
        (ENTRY_WITHDRAWN, 'Entry Withdrawn'),
    )


class GeneSymbolAliasSource:
    NCBI = 'N'
    UCSC = "U"
    HGNC = 'H'
    MANUAL = 'M'

    CHOICES = (
        (NCBI, 'NCBI'),
        (UCSC, 'UCSC'),
        (HGNC, 'HGNC'),
        (MANUAL, 'Manual'),
    )
