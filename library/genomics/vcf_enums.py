class VCFColumns:
    CHROM = 0
    POS = 1
    ID = 2
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8


class VCFSymbolicAllele:
    CNV = "<CNV>"
    DEL = "<DEL>"
    DUP = "<DUP>"
    INS = "<INS>"
    INV = "<INV>"


class VCFConstant:
    FREEBAYES = "freeBayes"
    CLCAD2 = "CLCAD2"  # CLC Genomics Workbench - variant track counts for ref,alt (1 for each alt)
    DEFAULT_ALLELE_FIELD = 'AD'
    DEFAULT_ALLELE_FREQUENCY_FIELD = "AF"
    DEFAULT_READ_DEPTH_FIELD = 'DP'
    DEFAULT_GENOTYPE_FIELD = 'GT'
    DEFAULT_GENOTYPE_QUALITY_FIELD = 'GQ'
    DEFAULT_PHRED_LIKILIHOOD_FIELD = 'PL'
    DEFAULT_SAMPLE_FILTERS_FIELD = 'FT'
    GENOTYPE_LIKELIHOOD = "GL"
    ALT_DEPTH_FIELD = "AO"  # FreeBayes - Alternate allele observation count
    REF_DEPTH_FIELD = "RO"  # FreeBayes - Reference allele observation count
