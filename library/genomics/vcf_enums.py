"""
VCF-level constants that both the importer and the writers share: column positions (VCFColumns),
the symbolic alts (VCFSymbolicAllele, and GeneLevelSymbolicAlt for gene-level events on the fake
contig), header constants, and VariantClass (Ensembl's VARIANT_CLASS terms).
"""
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
        @see genes.models.GeneLevelId for what to send instead when a record leaves. """
    HGNC = "HGNC", "HGNC ID"
    GENE = "GENE", "Local gene ID"


class GeneLevelSymbolicAlt(models.TextChoices):
    """ Symbolic alts for gene-level events, which live on the shared gene-level contig with a
        GeneLevelId as position. The alt carries the id the event is about, so biological identity
        hashes to its own Sequence and therefore its own Variant.

        Encoding identity in the alt is what lets the existing (locus, alt, svlen) unique constraint do
        the work - @see snpdb.gene_level_variants for why these are Variants at all.

        FUSION is directional - the anchor is the 5' partner, so BCR-ABL1 and ABL1-BCR are distinct.
        FUSION_UNORDERED anchors on the smaller id, because an unordered report asserts no direction.
        GAIN/LOSS are a whole-gene copy number call with no coordinates at all: they repeat the
        position's own gene, so the alt alone says what the variant is. GAIN rather than AMP because
        the threshold that makes a gain an amplification is the lab's, not ours.

        SPLICE repeats the gene too, and carries the name of the junction as a third segment
        (<SPLICE:HGNC:644:V7>), so two events in one gene are two variants. @see genes.gene_splice
        for what the labels are and where the coordinates go. """

    FUSION = "FUSION", "Gene fusion"
    FUSION_UNORDERED = "FUSION_UNORDERED", "Gene fusion (direction not asserted)"
    GAIN = "GAIN", "Gene copy number gain"
    LOSS = "LOSS", "Gene copy number loss"
    SPLICE = "SPLICE", "Splice event"

    # A partner the caller left unspecified, in place of the namespace:id
    UNKNOWN_PARTNER = Constant("UNKNOWN")

    @staticmethod
    def format(kind: str, namespace: Optional[str], gene_id: Optional[int],
               label: Optional[str] = None) -> str:
        if gene_id is None:
            return f"<{kind}:{GeneLevelSymbolicAlt.UNKNOWN_PARTNER}>"
        if label is not None:
            return f"<{kind}:{namespace}:{gene_id}:{label}>"
        return f"<{kind}:{namespace}:{gene_id}>"

    @staticmethod
    def parse(alt) -> Optional[tuple[str, Optional[str], Optional[int], Optional[str]]]:
        """ Returns (kind, namespace, gene id, label) - namespace/id are None for an unknown partner
            and the label is None for every kind but SPLICE. None if this isn't a gene-level alt """
        if m := GENE_LEVEL_ALT_PATTERN.fullmatch(str(alt)):
            kind = m.group("splice_kind") or m.group("kind")
            namespace = m.group("splice_namespace") or m.group("namespace")
            gene_id = m.group("splice_gene_id") or m.group("gene_id")
            return kind, namespace, int(gene_id) if gene_id else None, m.group("label")
        return None


# What a splice junction is named by - a seeded label (V7, vIII, ex14skip) or, for a junction we have
# no name for, its own coordinates (X_66905968_66914514). @see genes.gene_splice
GENE_LEVEL_LABEL = r"[A-Za-z0-9._-]+"
# The kinds whose alt is <KIND:NAMESPACE:id>, longest first so FUSION_UNORDERED isn't shadowed by
# FUSION. SPLICE is its own branch, as it alone carries the label segment
_UNLABELLED_KINDS = sorted((v for v in GeneLevelSymbolicAlt.values if v != GeneLevelSymbolicAlt.SPLICE),
                           key=len, reverse=True)
_NAMESPACES = "|".join(GeneIdNamespace.values)
GENE_LEVEL_ALT_PATTERN = re.compile(
    rf"<(?:"
    rf"(?P<splice_kind>{GeneLevelSymbolicAlt.SPLICE}):(?P<splice_namespace>{_NAMESPACES}):"
    rf"(?P<splice_gene_id>\d+):(?P<label>{GENE_LEVEL_LABEL})"
    rf"|(?P<kind>{'|'.join(_UNLABELLED_KINDS)}):"
    rf"(?:(?P<namespace>{_NAMESPACES}):(?P<gene_id>\d+)|{GeneLevelSymbolicAlt.UNKNOWN_PARTNER})"
    rf")>"
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
    # The FORMAT (or, for a single-sample VCF, INFO) keys a caller writes copy number under, best
    # first: CN is an integer copy number, SM a linear copy ratio (DRAGEN), FC a fold change (Pisces)
    COPY_NUMBER_FIELDS = ("CN", "SM", "FC")
    # Which of those is a ratio against the normal rather than an absolute count. They are different
    # quantities, so a classification stores each under its own evidence key (copy_number / fold_change)
    COPY_NUMBER_FIELD_IS_RATIO = {"CN": False, "SM": True, "FC": True}


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
    # Not an Ensembl class - VEP has none for a fusion, and a gene-level variant never reaches it anyway.
    # SO:0001565, the term SnpEff emits. @see snpdb.gene_level_variants
    GENE_FUSION = 'GF', "gene_fusion"
    # Nor for a splice event reported as a junction rather than a variant in a splice site -
    # SO:0001568. @see genes.gene_splice
    SPLICING_VARIANT = 'SP', "splicing_variant"


# Presentation grouping for variant type filters - every VariantClass belongs to exactly one group.
# A new VEP class is an enum member plus a line here.
VARIANT_CLASS_GROUPS = {
    "SNV": [
        VariantClass.SNV,
    ],
    "Indel": [
        VariantClass.INSERTION,
        VariantClass.DELETION,
        VariantClass.INDEL,
        VariantClass.SUBSTITUTION,
        VariantClass.COMPLEX_SUBSTITUTION,
        VariantClass.SEQUENCE_ALTERATION,
    ],
    "Copy number": [
        VariantClass.COPY_NUMBER_GAIN,
        VariantClass.COPY_NUMBER_LOSS,
        VariantClass.COPY_NUMBER_VARIATION,
        VariantClass.DUPLICATION,
        VariantClass.TANDEM_DUPLICATION,
    ],
    "Rearrangement": [
        VariantClass.INVERSION,
        VariantClass.TRANSLOCATION,
        VariantClass.INTERCHROMOSOMAL_TRANSLOCATION,
        VariantClass.INTRACHROMOSOMAL_TRANSLOCATION,
        VariantClass.CHROMOSOME_BREAKPOINT,
        VariantClass.INTERCHROMOSOMAL_BREAKPOINT,
        VariantClass.INTRACHROMOSOMAL_BREAKPOINT,
        VariantClass.COMPLEX_STRUCTURAL_ALTERATION,
        VariantClass.LOSS_OF_HETEROZYGOSITY,
    ],
    "Fusion": [
        VariantClass.GENE_FUSION,
    ],
    "Splicing": [
        VariantClass.SPLICING_VARIANT,
    ],
    "Other": [
        VariantClass.ALU_INSERTION,
        VariantClass.HERV_INSERTION,
        VariantClass.LINE1_INSERTION,
        VariantClass.SVA_INSERTION,
        VariantClass.MOBILE_ELEMENT_INSERTION,
        VariantClass.MOBILE_ELEMENT_DELETION,
        VariantClass.ALU_DELETION,
        VariantClass.HERV_DELETION,
        VariantClass.LINE1_DELETION,
        VariantClass.SVA_DELETION,
        VariantClass.NOVEL_SEQUENCE_INSERTION,
        VariantClass.TANDEM_REPEAT,
        VariantClass.SHORT_TANDEM_REPEAT_VARIATION,
        VariantClass.GENETIC_MARKER,
        VariantClass.PROBE,
    ],
}


INFO_LIFTOVER_SWAPPED_REF_ALT = "VG_LIFTOVER_SWAPPED_REF_ALT"

# FILTER values the source header never declared. vcf_clean_and_filter moves them here because bcftools
# norm dies on an undeclared FILTER, and the genotype processor puts them back at insert.
# Separator is '|' as the VCF spec already bars whitespace and ';' from FILTER IDs, and ',' would read
# as a multi-value INFO
UNDECLARED_FILTERS_INFO = "VG_UNDECLARED_FILTERS"
UNDECLARED_FILTERS_SEPARATOR = "|"
