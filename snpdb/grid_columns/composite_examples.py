"""
Fictional grid rows the annotation descriptions page draws its example composite cells from.

A row holds what the *client* receives - keyed by each member's `variant_column`, carrying the value
after any server side rendering (choice fields hold the label the grid shows, not the code). Allele
frequencies are the exception: they are stored here as unit values and the page runs them through the
same formatter the grid does, so the example follows the deployment's percent setting.
@see AF_UNIT_COLUMNS in snpdb/grids.py

Nothing here is real data - gene EXAMPL1, transcript NM_012345.6, sample NA00001, lab "Example Lab".
Numbers are plausible and internally consistent, so a reader can follow one cell's headline into its
hover detail. snpdb/tests/test_composite_examples.py fails if a composite gains or loses a member
without its example keeping up.
"""
from typing import Any


def _va(**values) -> dict[str, Any]:
    """ variantannotation__ paths, keyed by the field name each ends in """
    return {f"variantannotation__{k}": v for k, v in values.items()}


def _clinvar(**values) -> dict[str, Any]:
    return {f"clinvar__{k}": v for k, v in values.items()}


def _dbnsfp_gene(**values) -> dict[str, Any]:
    return {f"variantannotation__gene__geneannotation__dbnsfp_gene__{k}": v for k, v in values.items()}


# Composite pk -> the grid row its cell draws
COMPOSITE_EXAMPLE_ROWS: dict[str, dict[str, Any]] = {
    "variant": {
        "id": 1,
        "locus__contig__name": "12",
        "locus__position": 12345678,
        "locus__ref__seq": "C",
        "alt__seq": "T",
        "svlen": None,  # not a symbolic alt, so there's no length to draw
        **_va(
            hgvs_c="NM_012345.6(EXAMPL1):c.1234C>T",
            hgvs_p="NP_012345.1:p.(Arg412Trp)",
            hgvs_g="NC_000012.12:g.12345678C>T",
            symbol="EXAMPL1",
        ),
    },
    # ClinVar has a germline LP with 2 stars and a somatic Tier III; internally one lab has called it
    # a VUS and nobody has classified it somatically - so that chip draws empty
    "classifications": {
        **_clinvar(
            highest_pathogenicity=4,
            review_status="Criteria provided - multiple submitters w/no conflicts",  # 2 stars
            clinical_significance="Likely pathogenic",
            conflicting_clinical_significance=None,
            clinvar_variation_id=123456,
            clinvar_allele_id=654321,
            origin=1,  # germline
            preferred_disease_name="Example syndrome",
            disease_database_name="MONDO:0000001",
            clinical_sources="Example Lab",
            drug_response=False,
            highest_oncogenicity=None,
            oncogenic_classification=None,
            oncogenic_review_status=None,
            somatic_tier="Tier III",
            somatic_review_status="Criteria provided - single submitter",  # 1 star
        ),
        "max_internal_classification": "3",
        "internally_classified": "3",
        "internally_classified_labs": "Example Lab",
        "max_internal_somatic_classification": None,
        "internally_classified_somatic": None,
    },
    "consequence_impact": _va(consequence="missense_variant", impact="MODERATE"),
    "gnomad": _va(
        gnomad_af=0.000123,
        gnomad_popmax="Non-Finnish European",
        gnomad_filtered=False,
        gnomad_popmax_af=0.00051,
        gnomad_ac=31,
        gnomad_an=251432,
        gnomad_hom_alt=0,
        gnomad_hemi_count=0,
        gnomad_popmax_ac=29,
        gnomad_popmax_an=56780,
        gnomad_popmax_hom_alt=0,
        gnomad_afr_af=0.00002,
        gnomad_amr_af=0.00009,
        gnomad_asj_af=0.00003,
        gnomad_eas_af=0.00001,
        gnomad_fin_af=0.00004,
        gnomad_mid_af=0.00006,
        gnomad_nfe_af=0.00051,
        gnomad_oth_af=0.0001,
        gnomad_sas_af=0.00001,
        gnomad_xy_af=0.00012,
        gnomad_xy_ac=15,
        gnomad_xy_an=124800,
        gnomad_non_par=True,
        gnomad_faf95=0.00009,
        gnomad_faf99=0.00008,
        gnomad_fafmax_faf95_max=0.00038,
        gnomad_fafmax_faf99_max=0.00034,
        gnomad2_liftover_af=0.00011,
    ),
    # TOPMed's larger cohort has seen it more often than the other two
    "pop_freq_other": _va(af_1kg=0.0002, af_uk10k=0.00013, topmed_af=0.00031),
    # An acceptor loss 3 bases into the intron carries the max delta score
    "spliceai": _va(
        spliceai_max_ds=0.62,
        spliceai_pred_ds_ag=0.02,
        spliceai_pred_ds_al=0.62,
        spliceai_pred_ds_dg=0.0,
        spliceai_pred_ds_dl=0.01,
        spliceai_pred_dp_ag=-12,
        spliceai_pred_dp_al=3,
        spliceai_pred_dp_dg=41,
        spliceai_pred_dp_dl=-7,
        spliceai_gene_symbol="EXAMPL1",
    ),
    # % diff/ref is 100 * diff / ref
    "maxentscan": _va(maxentscan_percent_diff_ref=-38.4, maxentscan_ref=9.21,
                      maxentscan_alt=5.67, maxentscan_diff=-3.54),
    "mastermind": _va(mastermind_count_1_cdna=7, mastermind_count_2_cdna_prot=9,
                      mastermind_count_3_aa_change=14, mastermind_mmid3="EXAMPL1:R412W"),
    "aloft": _va(aloft_pred="Recessive", aloft_prob_dominant=0.12, aloft_prob_recessive=0.81,
                 aloft_prob_tolerant=0.07, aloft_high_confidence=True,
                 aloft_ensembl_transcript="ENST00000123456"),
    # Four tools call it damaging (SIFT, Polyphen2, MutationTaster, MetaLR), two benign
    "predictions": _va(
        predictions_num_pathogenic=4,
        predictions_num_benign=2,
        sift="Damaging",
        polyphen2_hvar_pred_most_damaging="Probably Damaging",
        mutation_taster_pred_most_damaging="Disease Causing",
        mutation_assessor_pred_most_damaging="Low",
        fathmm_pred_most_damaging="Tolerated",
        metalr_rankscore=0.91,
    ),
    "db_zygosity": {
        "global_variant_zygosity__het_count": 3,
        "global_variant_zygosity__hom_count": 1,
        "global_variant_zygosity__ref_count": 42,
        "global_variant_zygosity__unk_count": 0,
    },

    "alphamissense": _va(alphamissense_pred="likely_pathogenic", alphamissense_score=0.87,
                         alphamissense_rankscore=0.92),
    "bayesdel": _va(bayesdel_noaf_score=0.31, bayesdel_noaf_rankscore=0.88),
    "cadd": _va(cadd_phred=24.3, cadd_raw=3.91, cadd_raw_rankscore=0.93),
    "clinpred": _va(clinpred_pred="Damaging", clinpred_score=0.96, clinpred_rankscore=0.94),
    "eve": _va(eve_class="Pathogenic", eve_score=0.81),
    "metarnn": _va(metarnn_pred="Damaging", metarnn_score=0.89),
    "primateai": _va(primateai_pred="Damaging", primateai_score=0.83),
    "revel": _va(revel_score=0.78, revel_rankscore=0.90),
    "vest4": _va(vest4_score=0.86, vest4_rankscore=0.91),
    "varity": _va(varity_r_score=0.74, varity_er_score=0.69),
    "mutpred2": _va(mutpred2_score=0.71,
                    mutpred2_top5_mechanisms="Loss of helix (P = 0.02); Gain of loop (P = 0.04); "
                                             "Altered stability (P = 0.08); Loss of catalytic site (P = 0.11); "
                                             "Altered ordered interface (P = 0.15)"),
    "hipred": _dbnsfp_gene(hipred_prediction=True, hipred_score=0.72),
    "gene_damage_index": _dbnsfp_gene(gene_damage_index_phred=3.2, gene_damage_index_score=0.41),
    "gene_indispensability": _dbnsfp_gene(gene_indispensability_pred=True,
                                          gene_indispensability_score=0.88),
    "gnomad_constraint": {
        **_dbnsfp_gene(gnomad_pli=0.97, gnomad_pnull=0.01, gnomad_prec=0.02),
        "variantannotation__gene__geneannotation__gnomad_oe_lof": 0.12,
    },
    "uniprot": {
        "variantannotation__transcript_version__gene_version__hgnc__uniprot__accession": "P00001",
        "variantannotation__transcript_version__gene_version__hgnc__uniprot__function":
            "Catalyses the example reaction in the cytosol. Required for normal development of the "
            "example tissue",
        "variantannotation__transcript_version__gene_version__hgnc__uniprot__pathway":
            "Example biosynthesis; example-1 from example-2: step 1/3",
        "variantannotation__transcript_version__gene_version__hgnc__uniprot__tissue_specificity":
            "Expressed in liver and kidney",
        "variantannotation__transcript_version__gene_version__hgnc__uniprot__reactome":
            "R-HSA-000001;Example pathway",
    },
    # The overlap columns are text - AnnotSV reports them per overlapping SV
    "gnomad_sv": _va(gnomad_sv_overlap_af="0.0021", gnomad_sv_overlap_percent="94.5",
                     gnomad_sv_overlap_name="gnomAD-SV_v3_DEL_12_00001",
                     gnomad_sv_overlap_coords="12:12300000-12400000"),
    # Conserved on every score bar PhyloP 30 way, so the cell draws one muted dot among the filled ones
    "conservation": _va(
        phylop_100_way_vertebrate=7.85,
        phylop_30_way_mammalian=0.62,
        phylop_46_way_mammalian=2.03,
        phastcons_100_way_vertebrate=1.0,
        phastcons_30_way_mammalian=0.99,
        phastcons_46_way_mammalian=0.98,
        gerp_pp_rs=5.63,
    ),
    "cosmic": _va(cosmic_id="COSV100000001", cosmic_count=3, cosmic_legacy_id="COSM1000001"),
    "denovo_db": _va(denovo_db_case_count=2, denovo_db_control_count=0,
                     denovo_db_studies="Example2021,Example2023",
                     denovo_db_primary_phenotypes="developmental disorder",
                     denovo_db_pubmed_ids="10000001,10000002"),
    "open_targets": _va(
        open_targets_gwas_l2g_score=0.63,
        open_targets_gwas_diseases="Example trait",
        open_targets_gwas_gene_id="ENSG00000123456",
        open_targets_qtl_biosample="whole blood",
        open_targets_qtl_gene_id="ENSG00000123456",
        open_targets_study_id="GCST000001",
        open_targets_study_type="gwas",
        open_targets_is_lead="true",
        open_targets_variant_id="12_12345678_C_T",
        open_targets_gwas_l2g_scores="ENSG00000123456:0.63",
    ),
    "promoter_ai": _va(promoter_ai_score=0.42, promoter_ai_tss_pos=-85),
    "mavedb": _va(mavedb_score=-1.24, mavedb_urn="urn:mavedb:00000001-a-1"),
    "dbscsnv": _va(dbscsnv_ada_score=0.97, dbscsnv_rf_score=0.88),
    "nmd_ptc": _va(nmd_escape_status="Escapes NMD", nmd_escaping_variant=True,
                   ptc_distance_codons=12, ptc_last_junction_distance=-34),
    "protvar": _va(protvar_stability=-1.8, protvar_int="yes", protvar_pocket="no"),
    "essential_gene": _dbnsfp_gene(essential_gene_crispr="Essential",
                                   essential_gene_crispr2="Essential",
                                   essential_gene_gene_trap="Non-essential phenotype-changing"),
    "annotsv_acmg": _va(annotsv_acmg_class=4, annotsv_acmg_score=0.9,
                        annotsv_acmg_criteria="1A,2C,5F",
                        annotsv_pathogenic_overlaps="loss: ClinVar / Example syndrome / "
                                                    "HP:0000001 / chr12:12000000-12500000"),
    "annotsv_gene_effect": _va(annotsv_exons_spanned=3, annotsv_frameshift=True,
                               annotsv_re_gene="EXAMPL1"),
    "annotsv_benign_af": _va(annotsv_b_loss_af_max=0.012, annotsv_b_gain_af_max=0.004,
                             annotsv_b_ins_af_max=0.0, annotsv_b_inv_af_max=0.0),
    "annotsv_omim": _va(annotsv_omim_morbid=True, annotsv_omim_id="600001",
                        annotsv_omim_inheritance="AD", annotsv_omim_phenotype="Example syndrome"),
    "annotsv_region": _va(
        annotsv_repeat_type_left="AluY",
        annotsv_repeat_type_right="L1",
        annotsv_segdup_left="chr12:12000000-12010000",
        annotsv_segdup_right=None,
        annotsv_encode_blacklist_left=None,
        annotsv_encode_blacklist_right="chr12:12400000-12410000",
        annotsv_encode_blacklist_characteristics_left=None,
        annotsv_encode_blacklist_characteristics_right="High Signal Region",
    ),
    "annotsv_splice": _va(annotsv_dist_nearest_ss=42, annotsv_nearest_ss_type="5'"),
}

# The sample call cell is built per sample by the analysis grid rather than catalogued, so its example
# carries its own column - the row is keyed by this prefix plus the CohortGenotype format column
SAMPLE_EXAMPLE_PREFIX = "sample_1_"
SAMPLE_EXAMPLE_SAMPLE_NAME = "NA00001"
SAMPLE_EXAMPLE_ROW: dict[str, Any] = {
    f"{SAMPLE_EXAMPLE_PREFIX}{column}": value for column, value in {
        "samples_zygosity": "HET",
        "samples_allele_frequency": 0.48,
        "samples_allele_depth": 23,
        "samples_read_depth": 48,
        "samples_genotype_quality": 99,
        "samples_phred_likelihood": 0,
        "samples_filters": "PASS",
    }.items()
}
