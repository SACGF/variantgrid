import re
from typing import Optional

from django.db import models

from library.utils import Constant


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


class GeneIdNamespace(models.TextChoices):
    """ Whether the number in a gene-level alt means anything outside this deployment.
        HGNC is the same gene everywhere; GENE is a local id for a symbol HGNC doesn't carry -
        @see genes.models.FusionGeneId for what to send instead when a record leaves. """
    HGNC = "HGNC", "HGNC ID"
    GENE = "GENE", "Local gene ID"


class GeneLevelSymbolicAlt(models.TextChoices):
    """ Symbolic alts for gene-level events (gene fusions), which live on the shared gene-level contig
        with the anchor gene's FusionGeneId as position. The alt carries the partner's id, so
        biological identity hashes to its own Sequence and therefore its own Variant.

        Encoding identity in the alt is what lets the existing (locus, alt, svlen) unique constraint do
        the work - @see snpdb.gene_level_variants for why these are Variants at all.

        FUSION is directional - the anchor is the 5' partner, so BCR-ABL1 and ABL1-BCR are distinct.
        FUSION_UNORDERED anchors on the smaller id, because an unordered report asserts no direction.
        AMP/LOSS are for callers reporting a gene-level copy event with no coordinates at all. """

    FUSION = "FUSION", "Gene fusion"
    FUSION_UNORDERED = "FUSION_UNORDERED", "Gene fusion (direction not asserted)"
    AMP = "AMP", "Gene amplification"
    LOSS = "LOSS", "Gene loss"

    # A partner the caller left unspecified, in place of the namespace:id
    UNKNOWN_PARTNER = Constant("UNKNOWN")

    @staticmethod
    def format(kind: str, namespace: Optional[str], gene_id: Optional[int]) -> str:
        if gene_id is None:
            return f"<{kind}:{GeneLevelSymbolicAlt.UNKNOWN_PARTNER}>"
        return f"<{kind}:{namespace}:{gene_id}>"

    @staticmethod
    def parse(alt) -> Optional[tuple[str, Optional[str], Optional[int]]]:
        """ Returns (kind, namespace, gene id) - namespace/id are None for an unknown partner.
            None if this isn't a gene-level alt at all """
        if m := GENE_LEVEL_ALT_PATTERN.fullmatch(str(alt)):
            kind, namespace, gene_id = m.groups()
            return kind, namespace, int(gene_id) if gene_id else None
        return None


# Longest kind first, so FUSION_UNORDERED isn't shadowed by FUSION
GENE_LEVEL_ALT_PATTERN = re.compile(
    rf"<({'|'.join(sorted(GeneLevelSymbolicAlt.values, key=len, reverse=True))}):"
    rf"(?:({'|'.join(GeneIdNamespace.values)}):(\d+)|{GeneLevelSymbolicAlt.UNKNOWN_PARTNER})>"
)


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
    """ Ensembl variant classes (VEP `VARIANT_CLASS` / VariationFeature class_SO_term).

        Human-readable reference (note: does NOT list every value VEP can emit - e.g.
        chromosome_breakpoint and tandem_repeat are missing from it):
          https://asia.ensembl.org/info/genome/variation/prediction/classification.html#classes

        Authoritative source for the values VEP actually outputs:
          - SVs:        Bio/EnsEMBL/Variation/Utils/Config.pm  %SO_TERMS  (moved here in rel 114;
                        previously defined inline in Parser.pm::get_SO_term)
          - SVTYPE map: Bio/EnsEMBL/VEP/Parser.pm  get_SO_term()  (e.g. BND -> chromosome_breakpoint)
          - small vars: Bio/EnsEMBL/Variation/Utils/Sequence.pm  SO_variation_class()

        Keep every existing member even if a current VEP release no longer emits it - older
        annotation data may already contain it.
    """

    SNV = 'SN', "SNV"
    GENETIC_MARKER = 'GM', "genetic_marker"
    SUBSTITUTION = 'SU', "substitution"
    TANDEM_REPEAT = 'TR', "tandem_repeat"
    ALU_INSERTION = 'AI', "Alu_insertion"
    HERV_INSERTION = 'HI', "HERV_insertion"
    LINE1_INSERTION = 'LI', "LINE1_insertion"
    SVA_INSERTION = 'VI', "SVA_insertion"
    COMPLEX_STRUCTURAL_ALTERATION = 'CA', "complex_structural_alteration"
    COMPLEX_SUBSTITUTION = 'CS', "complex_substitution"
    COPY_NUMBER_GAIN = 'CG', "copy_number_gain"
    COPY_NUMBER_LOSS = 'CL', "copy_number_loss"
    COPY_NUMBER_VARIATION = 'CN', "copy_number_variation"
    DUPLICATION = 'DU', "duplication"
    CHROMOSOME_BREAKPOINT = 'CH', "chromosome_breakpoint"
    INTERCHROMOSOMAL_BREAKPOINT = 'IB', "interchromosomal_breakpoint"
    INTERCHROMOSOMAL_TRANSLOCATION = 'IT', "interchromosomal_translocation"
    INTRACHROMOSOMAL_BREAKPOINT = 'CB', "intrachromosomal_breakpoint"
    INTRACHROMOSOMAL_TRANSLOCATION = 'CT', "intrachromosomal_translocation"
    INVERSION = 'IN', "inversion"
    LOSS_OF_HETEROZYGOSITY = 'LO', "loss_of_heterozygosity"
    MOBILE_ELEMENT_DELETION = 'MD', "mobile_element_deletion"
    ALU_DELETION = 'AD', "Alu_deletion"
    HERV_DELETION = 'HD', "HERV_deletion"
    LINE1_DELETION = 'LD', "LINE1_deletion"
    SVA_DELETION = 'VD', "SVA_deletion"
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

# FILTER values the source header never declared. vcf_clean_and_filter moves them here because bcftools
# norm dies on an undeclared FILTER, and the genotype processor puts them back at insert.
# Separator is '|' as the VCF spec already bars whitespace and ';' from FILTER IDs, and ',' would read
# as a multi-value INFO
UNDECLARED_FILTERS_INFO = "VG_UNDECLARED_FILTERS"
UNDECLARED_FILTERS_SEPARATOR = "|"
