from django.db import models


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


class VariantClass(models.TextChoices):
    """ https://asia.ensembl.org/info/genome/variation/prediction/classification.html#classes """

    SNV = 'SN', "SNV"
    GENETIC_MARKER = 'GM', "genetic_marker"
    SUBSTITUTION = 'SU', "substitution"
    TANDEM_REPEAT = 'TR', "tandem_repeat"
    ALU_INSERTION = 'AI', "Alu_insertion"
    COMPLEX_STRUCTURAL_ALTERATION = 'CA', "complex_structural_alteration"
    COMPLEX_SUBSTITUTION = 'CS', "complex_substitution"
    COPY_NUMBER_GAIN = 'CG', "copy_number_gain"
    COPY_NUMBER_LOSS = 'CL', "copy_number_loss"
    COPY_NUMBER_VARIATION = 'CN', "copy_number_variation"
    DUPLICATION = 'DU', "duplication"
    INTERCHROMOSOMAL_BREAKPOINT = 'IB', "interchromosomal_breakpoint"
    INTERCHROMOSOMAL_TRANSLOCATION = 'IT', "interchromosomal_translocation"
    INTRACHROMOSOMAL_BREAKPOINT = 'CB', "intrachromosomal_breakpoint"
    INTRACHROMOSOMAL_TRANSLOCATION = 'CT', "intrachromosomal_translocation"
    INVERSION = 'IN', "inversion"
    LOSS_OF_HETEROZYGOSITY = 'LO', "loss_of_heterozygosity"
    MOBILE_ELEMENT_DELETION = 'MD', "mobile_element_deletion"
    MOBILE_ELEMENT_INSERTION = 'MI', "mobile_element_insertion"
    NOVEL_SEQUENCE_INSERTION = 'NI', "novel_sequence_insertion"
    SHORT_TANDEM_REPEAT_VARIATION = 'ST', "short_tandem_repeat_variation"
    TANDEM_DUPLICATION = 'TD', "tandem_duplication"
    TRANSLOCATION = 'TL', "translocation"
    DELETION = 'DE', "deletion"
    INDEL = 'ND', "indel"
    INSERTION = 'IS', "insertion"
    SEQUENCE_ALTERATION = 'SA', "sequence_alteration"
    PROBE = 'PR', "probe"


INFO_LIFTOVER_SWAPPED_REF_ALT = "VG_LIFTOVER_SWAPPED_REF_ALT"
