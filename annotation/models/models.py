"""
The variant annotation rows VEP writes: VariantAnnotation (the representative transcript, one row per variant) and
VariantTranscriptAnnotation (every transcript) share AbstractVariantAnnotation's per-transcript columns, and
VariantGeneOverlap has a row per gene a variant overlaps, for gene list queries. All three are partitioned per
VariantAnnotationVersion, so query them through annotation/annotation_version_querysets.py, never
VariantAnnotation.objects directly. Plugin display parsing is in annotation/open_targets.py and
annotation/vep_plugin_display.py. Also the gene citation / PubMed count tables a CachedWebResource loads.
"""
import logging
import re
from collections import defaultdict
from collections.abc import Iterable
from functools import cached_property
from typing import Optional

from Bio.Data.IUPACData import protein_letters_1to3
from django.conf import settings
from django.db import models
from django.db.models import F, Max, Min, Q
from django.db.models.deletion import CASCADE, DO_NOTHING, SET_NULL
from django.db.models.functions import Coalesce, Greatest

from annotation import open_targets, vep_plugin_display
from annotation.enums import (
    ALoFTPrediction,
    AlphaMissensePrediction,
    ClinPredPrediction,
    FATHMMPrediction,
    MetaRNNPrediction,
    MutationAssessorPrediction,
    MutationTasterPrediction,
    NMDEscapeStatus,
    Pathogenicity,
    PathogenicityImpact,
    Polyphen2Prediction,
    PrimateAIPrediction,
    SIFTPrediction,
    VariantAnnotationPipelineType,
    VEPSkippedReason,
)
from annotation.external_search_terms import (
    get_variant_pubmed_search_terms,
    get_variant_search_terms,
)
from annotation.models.models_citations import Citation
from annotation.models.models_version import VariantAnnotationVersion
from annotation.models.repeat_masker import RepeatMaskerSummary
from genes.enums import AnnotationConsortium
from genes.hgvs.hgvs_matcher import HGVSMatcher
from genes.models import (
    Gene,
    GeneSymbol,
    Transcript,
    TranscriptVersion,
)
from library.django_utils.data_archive_mixin import DataArchivedError
from library.genomics.vcf_enums import VariantClass
from library.utils import invert_dict
from patients.enums import GnomADPopulation
from snpdb.models import Variant
from snpdb.models.models_cached_web_resource import CachedWebResource


class AbstractVariantAnnotation(models.Model):
    """ Common fields between VariantAnnotation and VariantTranscriptAnnotation
        These fields are PER-TRANSCRIPT """
    SV_HGVS_TOO_LONG_MESSAGE = "HGVS not calculated due to length"
    SV_HGVS_ERROR_MESSAGE = "Error creating HGVS"

    # DO_NOTHING as records are deleted by dropping the partition (a CASCADE locks every other version's)
    version = models.ForeignKey(VariantAnnotationVersion, on_delete=DO_NOTHING)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    annotation_run = models.ForeignKey("AnnotationRun", on_delete=CASCADE)

    # gene/transcript are set by VEP and don't have a version.
    # can be set when up/downstream of gene (in which case HGVS_C and thus transcript_version will be null)
    gene = models.ForeignKey(Gene, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=SET_NULL)
    # Linked from HGVS transcript (to get version) @see BulkVEPVCFAnnotationInserter.get_transcript_version_id
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)

    # VEP Fields
    # The best way to see how these map to VEP fields is via the annotation details page
    alphamissense_pred = models.CharField(max_length=1, choices=AlphaMissensePrediction.choices, null=True, blank=True)
    alphamissense_score = models.FloatField(null=True, blank=True)
    amino_acids = models.TextField(null=True, blank=True)
    bayesdel_noaf_score = models.FloatField(null=True, blank=True)
    cadd_phred = models.FloatField(null=True, blank=True)
    cadd_raw = models.FloatField(null=True, blank=True)
    # TODO: This doesn't need to be nullable (default=False) - but will be slow. Change with next schema change
    canonical = models.BooleanField(null=True, blank=True)
    clinpred_pred = models.CharField(max_length=1, choices=ClinPredPrediction.CHOICES, null=True, blank=True)
    clinpred_score = models.FloatField(null=True, blank=True)
    nmd_escaping_variant = models.BooleanField(null=True, blank=True)

    # #579 - where the frameshift's new premature termination codon lands, and what that
    # implies for NMD. @see annotation/ptc.py
    # Codons from the first changed residue to the new stop, read from hgvs_p "fsTer<n>".
    # 1 = the changed codon is itself the stop (e.g. p.Tyr535Ter)
    ptc_distance_codons = models.IntegerField(null=True, blank=True)
    # Nucleotides from the PTC to the last exon-exon junction.
    # Negative = PTC lies in the final exon. Drives the 50nt rule.
    ptc_last_junction_distance = models.IntegerField(null=True, blank=True)
    # PTC-aware NMD prediction - see NMDEscapeStatus
    nmd_escape_status = models.CharField(max_length=1, choices=NMDEscapeStatus.choices, null=True, blank=True)

    codons = models.TextField(null=True, blank=True)
    consequence = models.TextField(null=True, blank=True)
    distance = models.IntegerField(null=True, blank=True)
    domains = models.TextField(null=True, blank=True)
    ensembl_protein = models.TextField(null=True, blank=True)
    exon = models.TextField(null=True, blank=True)
    fathmm_pred_most_damaging = models.CharField(max_length=1, choices=FATHMMPrediction.CHOICES, null=True, blank=True)
    flags = models.TextField(null=True, blank=True)
    gerp_pp_rs = models.FloatField(null=True, blank=True)  # Genomic Evolutionary Rate Profiling (GERP++)
    grantham = models.IntegerField(null=True, blank=True)
    hgvs_c = models.TextField(null=True, blank=True)
    hgvs_p = models.TextField(null=True, blank=True)
    impact = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    interpro_domain = models.TextField(null=True, blank=True)
    intron = models.TextField(null=True, blank=True)
    maxentscan_alt = models.FloatField(null=True, blank=True)
    maxentscan_diff = models.FloatField(null=True, blank=True)
    maxentscan_ref = models.FloatField(null=True, blank=True)
    maxentscan_percent_diff_ref = models.FloatField(null=True, blank=True)
    metarnn_pred = models.CharField(max_length=1, choices=MetaRNNPrediction.CHOICES, null=True, blank=True)
    metarnn_score = models.FloatField(null=True, blank=True)
    mpc_score = models.FloatField(null=True, blank=True)
    mutation_assessor_pred_most_damaging = models.CharField(max_length=1, choices=MutationAssessorPrediction.CHOICES, null=True, blank=True)
    mutation_taster_pred_most_damaging = models.CharField(max_length=1, choices=MutationTasterPrediction.CHOICES, null=True, blank=True)
    mutpred2_score = models.FloatField(null=True, blank=True)
    mutpred2_top5_mechanisms = models.TextField(null=True, blank=True)
    polyphen2_hvar_pred_most_damaging = models.CharField(max_length=1, choices=Polyphen2Prediction.CHOICES, null=True, blank=True)
    primateai_pred = models.CharField(max_length=1, choices=PrimateAIPrediction.CHOICES, null=True, blank=True)
    primateai_score = models.FloatField(null=True, blank=True)
    # protein_position = text as it can be e.g. indel: "22-23" or splicing: "?-10" or "10-?"
    protein_position = models.TextField(null=True, blank=True)
    revel_score = models.FloatField(null=True, blank=True)
    sift = models.CharField(max_length=1, choices=SIFTPrediction.CHOICES, null=True, blank=True)
    splice_region = models.TextField(null=True, blank=True)
    symbol = models.TextField(null=True, blank=True)
    varity_er_score = models.FloatField(null=True, blank=True)
    varity_r_score = models.FloatField(null=True, blank=True)
    vest4_score = models.FloatField(null=True, blank=True)

    mavedb_score = models.FloatField(null=True, blank=True)
    mavedb_urn = models.TextField(null=True, blank=True)

    # UTRAnnotator and SpliceVault (#2426). Rows annotated before the plugins were switched on are null
    # until backfilled (VariantAnnotationVersion.backfilled_utr_avi_splicevault), and the backfill only
    # reaches VariantAnnotation, never VariantTranscriptAnnotation.

    # UTRAnnotator - 5' UTR variants only
    utr_annotator_consequence = models.TextField(null=True, blank=True)  # eg 5_prime_UTR_premature_start_codon_gain_variant
    # One record per uORF effect, '&'-joined, each 'key=value:key=value' (KozakStrength, DistanceToCDS, type, ...)
    utr_annotator_annotation = models.TextField(null=True, blank=True)
    utr_annotator_existing_uorfs = models.IntegerField(null=True, blank=True)
    utr_annotator_existing_out_of_frame_oorfs = models.IntegerField(null=True, blank=True)
    utr_annotator_existing_inframe_oorfs = models.IntegerField(null=True, blank=True)

    # SpliceVault - RNA-seq mis-splicing events at a splice site SpliceAI predicts lost. Keyed by MANE
    # transcript, so only that transcript's entry has values (annotation/splicevault.py)
    # One record per event, '&'-joined, each 'rank:type:transcript_impact:percent_of_samples:frame'
    splicevault_top_events = models.TextField(null=True, blank=True)
    splicevault_site_sample_count = models.IntegerField(null=True, blank=True)
    splicevault_site_max_depth = models.IntegerField(null=True, blank=True)
    splicevault_out_of_frame_events = models.TextField(null=True, blank=True)  # eg Frameshift:2/4
    splicevault_site_pos = models.TextField(null=True, blank=True)  # eg chr1:925801
    splicevault_site_type = models.TextField(null=True, blank=True)  # Donor_loss / Acceptor_loss
    splicevault_spliceai_delta = models.FloatField(null=True, blank=True)

    class Meta:
        abstract = True

    @property
    def has_hgvs_c(self) -> bool:
        return self.hgvs_c and self.hgvs_c != self.SV_HGVS_TOO_LONG_MESSAGE

    @property
    def transcript_accession(self):
        """ Get transcript_id (with version if possible) """
        if self.transcript_version:
            return self.transcript_version.accession
        if self.transcript:
            return self.transcript.identifier
        return None

    @staticmethod
    def get_domains_components(domains):
        domain_dict = defaultdict(list)
        for d in domains.split("&"):
            key, *values = d.split(":")
            for v in values:
                domain_dict[key].append(v)
        return domain_dict

    @staticmethod
    def amino_acid_3_to_1(protein_string: str) -> str:
        """ p.Ala2351Thr -> p.A2351T """

        aa_3_to_1 = invert_dict(protein_letters_1to3)
        aa_3_to_1["Ter"] = "*"  # BioPython omits stop codon; HGVS uses Ter → *
        aa_3_regex_or = "|".join(aa_3_to_1)

        protein_aa1 = protein_string
        for aa3 in re.findall(f"({aa_3_regex_or})", protein_string):
            protein_aa1 = protein_aa1.replace(aa3, aa_3_to_1[aa3])
        return protein_aa1

    @staticmethod
    def protein_position_to_int(protein_position: str) -> int:
        """ VEP protein position can be:
            185, 185-187, ?-187, 185-? Takes 1st integer it finds """

        for i in protein_position.split("-"):
            try:
                return int(i)
            except ValueError:
                pass  # Skip "?"
        raise ValueError(f"Unable to handle protein_position={protein_position}")

    @cached_property
    def is_gene_level_annotation(self) -> bool:
        """ Computed from the gene identity rather than by VEP - @see snpdb.gene_level_variants.
            Everything VEP produces is absent, so the detail page hides those sections """
        return self.annotation_run.pipeline_type == VariantAnnotationPipelineType.GENE_LEVEL

    def get_hgvs_c_with_symbol(self) -> str:
        hgvs_c = self.hgvs_c
        if self.is_gene_level_annotation:
            # Already VICC gene-level nomenclature, which names both genes - there is no HGVS here to
            # parse, and nothing to graft a symbol onto. @see genes.models.fusion_canonical_str
            return hgvs_c
        if self.has_hgvs_c and self.symbol:
            hgvs_matcher = HGVSMatcher.instance(self.version.genome_build)
            hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_c)
            hgvs_variant.gene = self.symbol
            hgvs_c = str(hgvs_variant)
        return hgvs_c

    def get_short_label(self) -> str:
        """ Shortest label that still tells a human which variant this is, eg "RUNX1:p.Ala547Val" -
            for tabs and other places with no room for the transcript """
        change = None
        if self.hgvs_p:
            change = self.hgvs_p.split(":", 1)[-1]
        elif self.has_hgvs_c:
            change = self.hgvs_c.split(":", 1)[-1]

        if change is None:
            return str(self.variant)
        if self.symbol:
            return f"{self.symbol}:{change}"
        return change


class VariantAnnotation(AbstractVariantAnnotation):
    """ This is the "representative transcript" chosen (1 per variant/annotation version) """
    GENE_COLUMN = "variantannotation__gene"

    # Only need this once per variant
    hgvs_g = models.TextField(null=True, blank=True)

    # Population frequency
    af_1kg = models.FloatField(null=True, blank=True)
    af_uk10k = models.FloatField(null=True, blank=True)
    topmed_af = models.FloatField(null=True, blank=True)
    gnomad_af = models.FloatField(null=True, blank=True)
    gnomad2_liftover_af = models.FloatField(null=True, blank=True)
    gnomad_ac = models.IntegerField(null=True, blank=True)
    gnomad_an = models.IntegerField(null=True, blank=True)
    gnomad_hom_alt = models.IntegerField(null=True, blank=True)
    gnomad_afr_af = models.FloatField(null=True, blank=True)
    gnomad_amr_af = models.FloatField(null=True, blank=True)
    gnomad_asj_af = models.FloatField(null=True, blank=True)
    gnomad_eas_af = models.FloatField(null=True, blank=True)
    gnomad_fin_af = models.FloatField(null=True, blank=True)
    gnomad_mid_af = models.FloatField(null=True, blank=True)  # Middle East - added in gnomADv4
    gnomad_nfe_af = models.FloatField(null=True, blank=True)
    gnomad_oth_af = models.FloatField(null=True, blank=True)
    gnomad_sas_af = models.FloatField(null=True, blank=True)
    # filtering allele frequencies (new in gnomADv4)
    gnomad_faf95 = models.FloatField(null=True, blank=True)
    gnomad_faf99 = models.FloatField(null=True, blank=True)
    gnomad_fafmax_faf95_max = models.FloatField(null=True, blank=True)
    gnomad_fafmax_faf99_max = models.FloatField(null=True, blank=True)
    gnomad_xy_af = models.FloatField(null=True, blank=True)
    gnomad_xy_ac = models.IntegerField(null=True, blank=True)
    gnomad_xy_an = models.IntegerField(null=True, blank=True)
    gnomad_hemi_count = models.IntegerField(null=True, blank=True)  # This is set from gnomad_xy_ac if gnomad_non_par
    gnomad_popmax_af = models.FloatField(null=True, blank=True)
    gnomad_popmax_ac = models.IntegerField(null=True, blank=True)
    gnomad_popmax_an = models.IntegerField(null=True, blank=True)
    gnomad_popmax_hom_alt = models.IntegerField(null=True, blank=True)
    gnomad_filtered = models.BooleanField(null=True, blank=True)
    gnomad_non_par = models.BooleanField(null=True, blank=True)  # Not pseudoautosomal regions
    gnomad_popmax = models.CharField(max_length=3, choices=GnomADPopulation.choices, null=True, blank=True)

    # These are populated for SVs using VEP plugin StructuralVariantOverlap
    # They are text as they can have multiple entries joined via '&'
    # SVOverlapProcessor picks one of those records (ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD) and copies its
    # values onto the regular gnomad_af/gnomad_ac/gnomad_popmax_af/per-population columns, so that the analysis
    # PopulationNode (which compares floats) filters SVs the same way it does small variants. These text fields
    # keep every overlapping record, for display and to show which one was used.
    gnomad_sv_overlap_af = models.TextField(null=True, blank=True)
    gnomad_sv_overlap_percent = models.TextField(null=True, blank=True)
    gnomad_sv_overlap_name = models.TextField(null=True, blank=True)
    gnomad_sv_overlap_coords = models.TextField(null=True, blank=True)
    # Can't use filters unfortunately due to VEP custom bug, @see https://github.com/Ensembl/ensembl-vep/issues/1646
    # gnomad_sv_overlap_filters = models.TextField(null=True, blank=True)

    # AnnotSV (full-line) annotations - SV-only, populated by the AnnotSV stage on
    # the STRUCTURAL_VARIANT pipeline. Per-gene split-line rows are deferred (#1533).
    annotsv_acmg_class = models.IntegerField(null=True, blank=True, choices=Pathogenicity.choices)
    annotsv_acmg_score = models.FloatField(null=True, blank=True)          # AnnotSV_ranking_score
    annotsv_re_gene = models.TextField(null=True, blank=True)              # RE_gene
    annotsv_repeat_type_left = models.TextField(null=True, blank=True)
    annotsv_repeat_type_right = models.TextField(null=True, blank=True)
    annotsv_segdup_left = models.TextField(null=True, blank=True)
    annotsv_segdup_right = models.TextField(null=True, blank=True)
    annotsv_encode_blacklist_left = models.TextField(null=True, blank=True)
    annotsv_encode_blacklist_right = models.TextField(null=True, blank=True)
    annotsv_encode_blacklist_characteristics_left = models.TextField(null=True, blank=True)
    annotsv_encode_blacklist_characteristics_right = models.TextField(null=True, blank=True)
    # Bundled population AFs from AnnotSV's BenignSV sources, per SV-type.
    annotsv_b_gain_af_max = models.FloatField(null=True, blank=True)       # B_gain_AFmax
    annotsv_b_loss_af_max = models.FloatField(null=True, blank=True)       # B_loss_AFmax
    annotsv_b_ins_af_max = models.FloatField(null=True, blank=True)        # B_ins_AFmax
    annotsv_b_inv_af_max = models.FloatField(null=True, blank=True)        # B_inv_AFmax

    # AnnotSV (cont.) - additional fields captured in columns_version 4 (#1040 / #1533).
    annotsv_acmg_criteria = models.TextField(null=True, blank=True)        # AnnotSV_ranking_criteria
    annotsv_frameshift = models.BooleanField(null=True)                    # Frameshift (yes/no)
    annotsv_exons_spanned = models.IntegerField(null=True, blank=True)     # Exons_spanned
    annotsv_dist_nearest_ss = models.IntegerField(null=True, blank=True)   # Dist_nearest_SS
    annotsv_nearest_ss_type = models.TextField(null=True, blank=True)      # Nearest_SS_type
    annotsv_omim_inheritance = models.TextField(null=True, blank=True)     # OMIM_inheritance
    annotsv_omim_morbid = models.BooleanField(null=True)                   # OMIM_morbid
    annotsv_omim_phenotype = models.TextField(null=True, blank=True)       # OMIM_phenotype
    annotsv_omim_id = models.TextField(null=True, blank=True)              # OMIM_ID
    # Reference-only summary of pathogenic-SV overlaps (ClinVar / dbVar / ClinGen / OMIM-morbid).
    # Folds 16 P_{event}_{phen,hpo,source,coord} TSV columns into one JSONB blob keyed by event-type.
    annotsv_pathogenic_overlaps = models.JSONField(null=True)

    # From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267638/
    # "optimum cutoff value identified in the ROC analysis (0.6)"
    dbscsnv_ada_score = models.FloatField(null=True, blank=True)
    dbscsnv_rf_score = models.FloatField(null=True, blank=True)
    dbsnp_rs_id = models.TextField(null=True, blank=True)
    cosmic_id = models.TextField(null=True, blank=True)  # COSV - Genomic Mutation ID
    cosmic_legacy_id = models.TextField(null=True, blank=True)  # COSM
    cosmic_count = models.IntegerField(null=True, blank=True)
    # denovo-db - https://denovo-db.gs.washington.edu/. '&'-separated parallel arrays
    # (one element per contributing record), counts split case vs. PrimaryPhenotype=control
    denovo_db_studies = models.TextField(null=True, blank=True)
    denovo_db_pubmed_ids = models.TextField(null=True, blank=True)
    denovo_db_primary_phenotypes = models.TextField(null=True, blank=True)
    denovo_db_case_count = models.IntegerField(null=True, blank=True)
    denovo_db_control_count = models.IntegerField(null=True, blank=True)
    # ProtVar (columns_version >= 5) - https://www.ebi.ac.uk/ProtVar/
    # ddG protein stability (>2 likely destabilising); pocket/int are '&'-joined compound blobs
    protvar_stability = models.FloatField(null=True, blank=True)  # ProtVar_stability (ddG)
    protvar_pocket = models.TextField(null=True, blank=True)  # id&score&MpLDDT&energy&buriedness&RoG&residues
    protvar_int = models.TextField(null=True, blank=True)  # protein_id&pDockQ
    # Open Targets (columns_version >= 5, GRCh38) - https://platform.opentargets.org/
    open_targets_gwas_l2g_score = models.FloatField(null=True, blank=True)  # Locus-to-gene score (max)
    # Raw '&'-joined per-record L2G scores, parallel to the other open_targets_* arrays
    open_targets_gwas_l2g_scores = models.TextField(null=True, blank=True)
    open_targets_gwas_gene_id = models.TextField(null=True, blank=True)  # Ensembl gene id
    open_targets_gwas_diseases = models.TextField(null=True, blank=True)
    open_targets_study_type = models.TextField(null=True, blank=True)
    open_targets_study_id = models.TextField(null=True, blank=True)  # external lookup key
    # 'true'/'false' per record - whether this variant is the credible set's lead, or just a member
    open_targets_is_lead = models.TextField(null=True, blank=True)
    open_targets_variant_id = models.TextField(null=True, blank=True)  # external lookup key
    open_targets_qtl_gene_id = models.TextField(null=True, blank=True)  # Ensembl gene id
    open_targets_qtl_biosample = models.TextField(null=True, blank=True)
    # EVE / popEVE (columns_version >= 5, GRCh38, VEP >= 116) - https://evemodel.org/
    eve_score = models.FloatField(null=True, blank=True)
    eve_class = models.TextField(null=True, blank=True)  # Benign / Uncertain / Pathogenic
    popeve_score = models.FloatField(null=True, blank=True)
    # PromoterAI (columns_version >= 5, GRCh38, VEP >= 116)
    promoter_ai_score = models.FloatField(null=True, blank=True)
    promoter_ai_tss_pos = models.IntegerField(null=True, blank=True)
    # AlphaGenome Variant Impact (#2426) - GRCh38 SNVs, VEP >= 116. Null before the backfill as above
    avi_phred = models.FloatField(null=True, blank=True)
    avi_raw = models.FloatField(null=True, blank=True)
    # AlphaGenome merged splicing score (#2433) - GRCh38 SNVs, non-commercial data
    # (ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED). Null before the backfill as above
    alphagenome_splicing = models.FloatField(null=True, blank=True)
    pubmed = models.TextField(null=True, blank=True)
    # Mastermind Cited Variants Reference. @see https://www.genomenon.com/cvr/
    mastermind_count_1_cdna = models.IntegerField(null=True, blank=True)
    mastermind_count_2_cdna_prot = models.IntegerField(null=True, blank=True)
    mastermind_count_3_aa_change = models.IntegerField(null=True, blank=True)
    mastermind_mmid3 = models.TextField(null=True, blank=True)  # gene:key for mastermind_count_3_aa_change

    # dbNSFP pathogenicity rank scores (variant only)
    cadd_raw_rankscore = models.FloatField(null=True, blank=True)
    revel_rankscore = models.FloatField(null=True, blank=True)
    bayesdel_noaf_rankscore = models.FloatField(null=True, blank=True)
    clinpred_rankscore = models.FloatField(null=True, blank=True)
    vest4_rankscore = models.FloatField(null=True, blank=True)
    metalr_rankscore = models.FloatField(null=True, blank=True)
    alphamissense_rankscore = models.FloatField(null=True, blank=True)

    # ALoFT (from dbNSFP)
    aloft_prob_tolerant = models.FloatField(null=True, blank=True)
    aloft_prob_recessive = models.FloatField(null=True, blank=True)
    aloft_prob_dominant = models.FloatField(null=True, blank=True)
    aloft_pred = models.CharField(max_length=1, choices=ALoFTPrediction.choices, null=True, blank=True)
    aloft_high_confidence = models.BooleanField(null=True, blank=True)
    aloft_ensembl_transcript = models.TextField(null=True, blank=True)  # Transcript of most damaging prediction chosen

    # Not all builds have all phylop/phastcons
    phylop_30_way_mammalian = models.FloatField(null=True, blank=True)
    phylop_46_way_mammalian = models.FloatField(null=True, blank=True)
    phylop_100_way_vertebrate = models.FloatField(null=True, blank=True)
    phastcons_30_way_mammalian = models.FloatField(null=True, blank=True)
    phastcons_46_way_mammalian = models.FloatField(null=True, blank=True)
    phastcons_100_way_vertebrate = models.FloatField(null=True, blank=True)

    # SpliceAI - @see https://pubmed.ncbi.nlm.nih.gov/30661751/
    # This has so many columns, perhaps a highest score and summary col of "AG: -5, AL: 35, DG: -5, DL: 17"
    # 0.2 = high recal, 0.5 = recommended, 0.8 = high precision
    spliceai_pred_dp_ag = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_al = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_dg = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_dl = models.IntegerField(null=True, blank=True)
    spliceai_pred_ds_ag = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_al = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_dg = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_dl = models.FloatField(null=True, blank=True)
    spliceai_max_ds = models.FloatField(null=True, blank=True, db_index=True)
    spliceai_gene_symbol = models.TextField(null=True, blank=True)

    repeat_masker = models.TextField(null=True, blank=True)
    overlapping_symbols = models.TextField(null=True, blank=True)
    # Summary of most_damaging fields for faster DamageNode queries
    predictions_num_pathogenic = models.IntegerField(default=0)
    predictions_num_benign = models.IntegerField(default=0)

    somatic = models.BooleanField(null=True, blank=True)
    variant_class = models.CharField(max_length=2, choices=VariantClass.choices, null=True, blank=True)
    vep_skipped_reason = models.CharField(max_length=1, choices=VEPSkippedReason.choices, null=True, blank=True)

    GNOMAD_FIELDS = {
        GnomADPopulation.AFRICAN_AFRICAN_AMERICAN: 'gnomad_afr_af',
        GnomADPopulation.LATINO: 'gnomad_amr_af',
        GnomADPopulation.ASHKENAZI_JEWISH: 'gnomad_asj_af',
        GnomADPopulation.EAST_ASIAN: 'gnomad_eas_af',
        GnomADPopulation.FINNISH: 'gnomad_fin_af',
        GnomADPopulation.MIDDLE_EASTERN: 'gnomad_mid_af',
        GnomADPopulation.NON_FINNISH_EUROPEAN: 'gnomad_nfe_af',
        GnomADPopulation.OTHER: 'gnomad_oth_af',
        GnomADPopulation.SOUTH_ASIAN: 'gnomad_sas_af',
    }

    ALOFT_FIELDS = {
        "aloft_pred": "Pred",
        "aloft_high_confidence": "High Confidence",
        "aloft_prob_tolerant": "Tol",
        "aloft_prob_recessive": "Rec",
        "aloft_prob_dominant": "Dom",
        "aloft_ensembl_transcript": "Transcript",
    }

    MASTERMIND_FIELDS = {
        "mastermind_count_1_cdna": "cDNA",
        "mastermind_count_2_cdna_prot": "cDNA/Prot",
        "mastermind_count_3_aa_change": "AA change",
        "mastermind_mmid3": "MMID3",
    }

    SPLICEAI_DS_DP = {
        "AG": ("spliceai_pred_ds_ag", "spliceai_pred_dp_ag"),
        "AL": ("spliceai_pred_ds_al", "spliceai_pred_dp_al"),
        "DG": ("spliceai_pred_ds_dg", "spliceai_pred_dp_dg"),
        "DL": ("spliceai_pred_ds_dl", "spliceai_pred_dp_dl"),
    }

    # 'conserved' is where a score starts to read as conserved - the variant grid conservation cell
    # fills a dot at or above it. GERP++ and PhyloP 46 way are the deleterious thresholds from
    # https://academic.oup.com/hmg/article/24/8/2125/651446, PhyloP 100 way and PhastCons the ones the
    # conservation filter node was built around (@see ConservationNode). PhyloP 30 way has no published
    # cutoff - 1.0 is ~75% of the way to its 1.312 ceiling, in line with the PhyloP 100 way cutoff
    CONSERVATION_SCORES = {
        "gerp_pp_rs": {
            # UCSC says RS scores range from a maximum of 6.18 down to a below-zero minimum, which we cap at -12.36
            "min": -12.36,
            "max": 6.18,
            "conserved": 4.4,
        },
        # BigWig stats obtained via kent-335 bigWigInfo
        "phylop_30_way_mammalian": {
            "mean": 0.097833,
            "min": -20.0,
            "max": 1.312,
            "std": 0.727453,
            "conserved": 1.0,
        },
        "phylop_46_way_mammalian": {
            "mean": 0.035934,
            "min": -13.796,
            "max": 2.941,
            "std": 0.779426,
            "conserved": 1.6,
        },
        "phylop_100_way_vertebrate": {
            "mean": 0.093059,
            "min": -20.0,
            "max": 10.003,
            "std": 1.036944,
            "conserved": 1.4,
        },
        "phastcons_30_way_mammalian": {
            "mean": 0.128025,
            "min": 0.0,
            "max": 1.0,
            "std": 0.247422,
            "conserved": 0.85,
        },
        "phastcons_46_way_mammalian": {
            "mean": 0.088576,
            "min": 0.0,
            "max": 1.0,
            "std": 0.210242,
            "conserved": 0.85,
        },
        "phastcons_100_way_vertebrate": {
            "mean": 0.101765,
            "min": 0.0,
            "max": 1.0,
            "std": 0.237072,
            "conserved": 0.85,
        }
    }

    @staticmethod
    def get_variant_annotation_q_list() -> list[Q]:
        """ Filters to describe variants that can be annotated - a method rather than a class constant,
            as the reference Q looks up a Sequence pk.
            Gene-level variants belong here - they get a VariantAnnotation row like anything else, just
            written by the GENE_LEVEL pipeline rather than VEP. Which pipeline claims them is
            pipeline_type_variant_q's business, and it subtracts them from both VEP types. """
        return [
            Variant.get_no_reference_q(),
            ~Q(alt__seq__in=['.', '*']),  # Exclude non-standard variants
        ]

    @cached_property
    def is_standard_annotation(self) -> bool:
        return self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STANDARD

    @cached_property
    def is_structural_variant_annotation(self) -> bool:
        return self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT

    @cached_property
    def repeat_masker_summary(self) -> RepeatMaskerSummary:
        """ Group the (possibly '&'-joined) RepeatMasker value by repeat class - see #1580 """
        return RepeatMaskerSummary.from_value(self.repeat_masker)

    @property
    def has_pathogenicity(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_dbnsfp(self) -> bool:
        """ dbNSFP only annotates non-synonymous SNVs in coding regions, so for indels,
            synonymous/intronic/intergenic/UTR variants it contributes nothing. Single check
            on cadd_raw_rankscore would catch ~all of them, but use OR across several fields
            to cover the rare case where dbNSFP omits CADD but emits another score. #1045 """
        return any(v is not None for v in (
            self.cadd_raw_rankscore,
            self.revel_rankscore,
            self.bayesdel_noaf_rankscore,
            self.alphamissense_rankscore,
            self.clinpred_rankscore,
            self.vest4_rankscore,
            self.metalr_rankscore,
            self.aloft_pred,
            self.gerp_pp_rs,
            self.cadd_phred,
        ))

    @cached_property
    def show_rankscores(self) -> bool:
        """ dbNSFP rankscores arrived in columns_version 2. Once raw scores exist (columns_version 4+)
            they're legacy - shown only when ANNOTATION_SHOW_LEGACY_RANKSCORES is set. columns_version
            2-3 carry only rankscores (no raw scores), so always show them. """
        cv = self.version.columns_version
        if cv < 2:
            return False
        return settings.ANNOTATION_SHOW_LEGACY_RANKSCORES or cv < 4

    @property
    def has_conservation(self) -> bool:
        """ Thanks to summary stats we can now do this in VEP112 """
        if self.is_gene_level_annotation:
            return False
        return self.is_standard_annotation or self.version.vep >= 112

    @property
    def has_splicing(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_gnomad(self) -> bool:
        """ A gnomAD site with AC=0 has a real af of 0.0, so test for a value not a truthy one -
            'absent from gnomAD' and 'present, never observed' are different answers """
        return self.gnomad_af is not None or self.gnomad2_liftover_af is not None

    @property
    def annotsv_acmg_clinical_significance(self) -> Optional[str]:
        """ AnnotSV's 5-tier ACMG_class as a clinical_significance ekey value, so SVs can be shown
            with the same pill as classifications """
        if self.annotsv_acmg_class is None:
            return None
        return Pathogenicity(self.annotsv_acmg_class).short_label

    @property
    def has_annotsv(self) -> bool:
        """ AnnotSV is opt-in and its own pipeline run (#720), so an SV can have VEP annotation with
            none of these columns filled in - annotated before it was enabled, or its run is yet to
            finish. Checks every annotsv_ column so it can't drift as fields are added. """
        return any(getattr(self, f.name) is not None
                   for f in self._meta.fields if f.name.startswith("annotsv_"))

    @cached_property
    def visible_columns(self) -> frozenset[str]:
        """ All VariantGrid columns populated for this annotation's build / pipeline / version / data files.
            Drives variant detail per-row show/hide (via `labelled visible_fields=`) off the same VEP_COLUMNS
            table that controls what annotation writes, so the two can't drift (#1148). Passing `vep_config`
            drops columns whose data file isn't configured - matching the `VariantAnnotationVersion.has_*`
            flags (e.g. PhastCons/PhyloP mammalian tracks). """
        return self.version.get_visible_columns(self.annotation_run.pipeline_type)

    @property
    def has_non_gnomad_population_frequency(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_cosmic(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_mavedb(self) -> bool:
        return self.is_standard_annotation and self.version.columns_version >= 3

    @property
    def protvar_int_parsed(self) -> Optional[dict]:
        return vep_plugin_display.protvar_int_parsed(self)

    @property
    def protvar_pocket_parsed(self) -> Optional[dict]:
        return vep_plugin_display.protvar_pocket_parsed(self)

    @property
    def utr_annotator_effects(self) -> list[dict[str, str]]:
        return vep_plugin_display.utr_annotator_effects(self.utr_annotator_annotation)

    @property
    def splicevault_top_events_list(self) -> list[str]:
        return [e for e in (self.splicevault_top_events or "").split("&") if e]

    @property
    def gnomad4_or_later(self) -> bool:
        return self.version.gnomad_major_version >= 4

    @property
    def has_hemi(self):
        return self.gnomad4_or_later and self.variant.locus.contig.name == 'X'

    @property
    def gnomad_url(self):
        GNOMAD2 = "gnomad_r2_1"
        GNOMAD3 = "gnomad_r3"
        GNOMAD4 = "gnomad_r4"

        gnomad_dataset = None
        if self.gnomad_af and self.version.gnomad.startswith("2.1"):
            gnomad_dataset = GNOMAD2
        elif self.version.gnomad.startswith("3"):
            if self.gnomad_af is not None:
                gnomad_dataset = GNOMAD3
        elif self.gnomad_af and self.version.gnomad.startswith("4"):
            gnomad_dataset = GNOMAD4

        url = None
        if gnomad_dataset:
            v = self.variant
            gnomad_variant = f"{v.locus.chrom}-{v.locus.position}-{v.locus.ref}-{v.alt}"
            url = vep_plugin_display.get_gnomad_url(gnomad_variant, gnomad_dataset)
        return url

    def get_mave_urls(self) -> dict[str, str]:
        """ dict of label: url
            MAVE team intend on creating a variant landing page, so we'll replace this with that when ready
        """
        return vep_plugin_display.mavedb_urn_to_urls(self.mavedb_urn)

    @property
    def mmid3_mastermind_urls(self) -> dict:
        return vep_plugin_display.get_mmid3_mastermind_urls(self.mastermind_mmid3)

    @property
    def denovo_db_records(self) -> list[dict]:
        return vep_plugin_display.denovo_db_records(self)

    @cached_property
    def open_targets_records(self) -> list[dict]:
        return open_targets.open_targets_records(self)

    @property
    def open_targets_gwas_genes(self) -> list[dict]:
        return open_targets.open_targets_gwas_genes(self)

    @property
    def open_targets_qtl_genes(self) -> list[dict]:
        return open_targets.open_targets_qtl_genes(self)

    @property
    def open_targets_variant(self) -> Optional[dict]:
        return open_targets.open_targets_variant(self)

    def has_spliceai(self):
        """ SpliceAI scoring a variant 0 for all 4 deltas is a prediction (no splicing impact),
            not missing data - only a null means VEP found no SpliceAI record """
        return any(ds is not None for ds in (self.spliceai_pred_ds_ag, self.spliceai_pred_ds_al,
                                             self.spliceai_pred_ds_dg, self.spliceai_pred_ds_dl))

    def highest_spliceai(self) -> int|None:
        values = []
        for (ds, _) in self.SPLICEAI_DS_DP.values():
            if v := getattr(self, ds):
                values.append(v)
        if values:
            return max(values)
        return None

    @classmethod
    def backfill_spliceai_max_ds(cls, version: 'VariantAnnotationVersion',
                                 chunk_size: int = 10_000) -> int:
        any_ds_set = (
            Q(spliceai_pred_ds_ag__isnull=False)
            | Q(spliceai_pred_ds_al__isnull=False)
            | Q(spliceai_pred_ds_dg__isnull=False)
            | Q(spliceai_pred_ds_dl__isnull=False)
        )
        base_qs = cls.objects.filter(version=version, spliceai_max_ds__isnull=True).filter(any_ds_set)
        bounds = base_qs.aggregate(lo=Min("pk"), hi=Max("pk"))
        lo, hi = bounds["lo"], bounds["hi"]
        if lo is None:
            return 0

        total = 0
        cursor = lo
        while cursor <= hi:
            nxt = cursor + chunk_size
            updated = (cls.objects
                       .filter(version=version, spliceai_max_ds__isnull=True,
                               pk__gte=cursor, pk__lt=nxt)
                       .filter(any_ds_set)
                       .update(spliceai_max_ds=Greatest(
                           Coalesce(F("spliceai_pred_ds_ag"), 0.0),
                           Coalesce(F("spliceai_pred_ds_al"), 0.0),
                           Coalesce(F("spliceai_pred_ds_dg"), 0.0),
                           Coalesce(F("spliceai_pred_ds_dl"), 0.0),
                       )))
            total += updated
            logging.info("backfill_spliceai_max_ds vav=%s chunk pk=[%s,%s) updated=%s",
                         version.pk, cursor, nxt, updated)
            cursor = nxt
        return total

    @staticmethod
    def get_gnomad_population_field(population):
        return VariantAnnotation.GNOMAD_FIELDS.get(population)

    @cached_property
    def transcript_annotation(self) -> list['VariantTranscriptAnnotation']:
        return self.variant.varianttranscriptannotation_set.filter(version=self.version)

    def get_search_terms(self):
        extra_terms = []
        if self.dbsnp_rs_id:
            extra_terms.append(self.dbsnp_rs_id)
        return get_variant_search_terms(self.transcript_annotation, extra_terms=extra_terms)

    def get_pubmed_search_terms(self):
        return get_variant_pubmed_search_terms(self.transcript_annotation)

    @property
    def gnomad_sv_overlap(self) -> list[dict]:
        return vep_plugin_display.get_gnomad_sv_overlap(self.__dict__, self.version.gnomad_sv)

    def get_short_label(self) -> str:
        """ Falls back to g.HGVS, which only this (per-variant) annotation has - calculated for a row
            VEP skipped (eg an SV over its size cap), which has none of its own """
        if self.hgvs_p or self.has_hgvs_c:
            return super().get_short_label()
        return self.hgvs_g or VariantAnnotation.get_hgvs_g(self.variant) or str(self.variant)

    @staticmethod
    def get_hgvs_g(variant: Variant) -> Optional[str]:
        """ This is very slow - only use on a page with just 1 variant, not in a loop """
        # They should all be the same so any will be ok - but take latest in case we fixed HGVS
        qs = variant.variantannotation_set.filter(hgvs_g__isnull=False).order_by("-version")
        data = qs.values_list("hgvs_g", flat=True)[:1]
        hgvs_g = None
        if data:
            hgvs_g = data[0]
        if hgvs_g is None:
            hgvs_g = VariantAnnotation._generate_hgvs_g(variant)
        return hgvs_g

    @staticmethod
    def get_hgvs_g_for_variant_ids(variant_ids: Iterable[int]) -> dict[int, Optional[str]]:
        """ get_hgvs_g for a page of variants in one query, rather than one per variant """
        variant_ids = set(variant_ids)
        qs = VariantAnnotation.objects.filter(variant__in=variant_ids, hgvs_g__isnull=False)
        qs = qs.order_by("variant", "-version").distinct("variant")
        hgvs_g_by_variant_id = dict(qs.values_list("variant", "hgvs_g"))
        if unannotated_ids := variant_ids - hgvs_g_by_variant_id.keys():
            unannotated_qs = Variant.objects.filter(pk__in=unannotated_ids).select_related("locus__contig", "locus__ref", "alt")
            for variant in unannotated_qs:
                hgvs_g_by_variant_id[variant.pk] = VariantAnnotation._generate_hgvs_g(variant)
        return hgvs_g_by_variant_id

    @staticmethod
    def _generate_hgvs_g(variant: Variant) -> Optional[str]:
        """ Reference variants have no annotation - so we'll have to fall back to generating it.
            None for what has no g.HGVS to write - <CNV>/<INS> and gene-level """
        if not variant.can_make_g_hgvs:
            return None
        matcher = HGVSMatcher.instance(variant.any_genome_build)
        return matcher.variant_to_g_hgvs(variant)

    def __str__(self):
        return f"{self.variant}: {self.version}"

    class Meta:
        unique_together = ("version", "variant")


class VariantTranscriptAnnotation(AbstractVariantAnnotation):
    """ This is ALL transcript predictions for a variant (>=1 per variant) """

    class Meta:
        unique_together = ("version", "variant", "transcript_version")

    @staticmethod
    def get_overlapping_genes_q(variant_annotation_version: VariantAnnotationVersion, gene_ids_qs) -> Q:
        """ Returns Q object where any transcript for a variant matches gene in genes
            Variants can overlap with multiple genes, and the VariantAnnotation (ie "pick" or representative annotation)
            may not be the one in the gene list. Thus we have to check everything that was in transcript annotation too
        """
        if variant_annotation_version.data_archived:
            raise DataArchivedError(variant_annotation_version)

        vto = VariantGeneOverlap.objects.filter(version=variant_annotation_version, gene__in=gene_ids_qs)
        # Use pk__in so we don't return multiple records per variant
        return Q(pk__in=vto.values_list("variant_id", flat=True))


class VariantGeneOverlap(models.Model):
    """ Created for every gene that overlaps a variant (per annotation version)
        Allows efficient gene list queries (10% larger than variant annotation) while handling multiple transcripts """
    # DO_NOTHING as records are deleted by dropping the partition (a CASCADE locks every other version's)
    version = models.ForeignKey(VariantAnnotationVersion, on_delete=DO_NOTHING)
    annotation_run = models.ForeignKey("AnnotationRun", on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)

    class Meta:
        unique_together = ('version', 'variant', 'gene')


class GeneSymbolCitation(models.Model):
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    citation = models.ForeignKey(Citation, null=True, on_delete=CASCADE)

    class Meta:
        unique_together = ('gene_symbol', 'citation')


class GenePubMedCount(models.Model):
    """ From https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz """
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    count = models.IntegerField()
    cached_web_resource = models.ForeignKey(CachedWebResource, on_delete=CASCADE)

    @staticmethod
    def get_count_and_note_for_gene_symbol(gene_symbol_id) -> tuple[int, str]:
        gene_symbol = GeneSymbol.objects.get(pk=gene_symbol_id)
        genes_qs = gene_symbol.get_genes().filter(annotation_consortium=AnnotationConsortium.REFSEQ)
        notes = []
        count = 0
        if gene_counts := list(GenePubMedCount.objects.filter(gene__in=genes_qs)):
            for gc in gene_counts:
                count = max(count, gc.count)
                notes.append(f"{gc.gene}={gc.count}")
            modified = gene_counts[0].retrieved_date
        else:
            if first_gpmc := GenePubMedCount.objects.first():
                count = 0
                notes.append(f"No results (could be no RefSeq gene for {gene_symbol=})")
                modified = first_gpmc.retrieved_date
            else:
                raise ValueError("No GenePubMedCount entries")

        retrieved = modified.date().isoformat()
        notes.insert(0, f"NCBI curated gene2pubmed retrieved {retrieved}. Homo Sapien only. ")
        return count, " ".join(notes)

    @property
    def retrieved_date(self):
        return self.cached_web_resource.created


class MutationalSignatureInfo(models.Model):
    signature_id = models.IntegerField(primary_key=True)
    cancer_types = models.TextField()
    proposed_aetiology = models.TextField()
    additional_mutational_features = models.TextField()
    comments = models.TextField()
