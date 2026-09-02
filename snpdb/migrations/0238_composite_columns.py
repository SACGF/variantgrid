"""
Composite grid cells become catalogue entries with their members recorded in the database.

Reverse is lossy: the composite rows and every member row go, and a collection keeps only whatever
members it still had - the columns collapsed into a composite are not put back.
"""
from django.db import migrations

from library.django_utils import bulk_insert_class_data
from library.django_utils.composite_columns import collapse_into_composite

_VARIANT = 'V'
_TRANSCRIPT = 'T'
_GENE = 'G'
_DATABASE = 'D'

# What each anchor column was before #686 drew a composite cell on it - it is a plain member again,
# and the composite row below takes over the #686 label, width and wording
_RESTORED_MEMBERS = {
    'consequence': (None, 'Consequence',
                    'Sequence Ontology consequence type of input variant allele on transcript'),
    'gnomad_af': (None, 'gnomAD AF',
                  'Allele Frequency (0-1) among all gnomAD genotypes (exome+genome)'),
    'gnomad_popmax_af': (None, 'gnomAD PopMax AF',
                         'Allele Frequency (0-1) in pop with highest AF (excluding Ashkenazi, Finnish, '
                         'and indeterminate ancestry) (exome+genome)'),
    'spliceai_max_ds': (None, 'SpliceAI max delta score',
                        'Highest of the four SpliceAI delta scores'),
    'maxentscan_percent_diff_ref': (None, 'MaxEntScan % diff/ref',
                                    '<a href="http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html">'
                                    'MaxEntScan Scores Human 5 splice sites.</a> 100 * maxentscan_diff / maxentscan_ref'),
    'mastermind_count_1_cdna': (None, 'MM variant article count',
                                "Count of <a href='https://www.genomenon.com/cvr/'>Mastermind</a> articles "
                                "with cDNA matches for this specific variant"),
    'aloft_pred': (None, 'ALoFT Pred',
                   '<a href="https://www.nature.com/articles/s41467-017-00443-5">ALoFT</a> final '
                   'classification, values can be Tolerant, Recessive or Dominant'),
    'predictions_num_pathogenic': (None, 'No. Path Pred pathogenic',
                                   'No. Pathogenic prediction tools that are pathogenic'),
    'total_db_het': (35, 'Het DB count', 'Database het count'),
}

# 'A' was never a ColumnAnnotationLevel - the arrange page raises on it. Its AlphaMissense siblings
# are transcript level
_FIXED_ANNOTATION_LEVELS = {
    'alphamissense_rankscore': _TRANSCRIPT,
}

# The catalogue of composite cells. 'members' is in sort_order, the first being the headline - the
# value the cell reads as, its default sort key, and what the generic renderer prints. 'detail_only'
# members are drawn on hover but are not offered in the header's sort menu.
# 'exists' rows are already in the catalogue and only gain members here.
_COMPOSITES = [
    {'pk': 'variant', 'exists': True,
     'members': ['chrom', 'position', 'ref', 'alt', 'svlen', 'hgvs_c', 'hgvs_p', 'hgvs_g', 'symbol'],
     'detail_only': ['chrom', 'position', 'ref', 'alt', 'svlen', 'hgvs_c', 'hgvs_p', 'hgvs_g', 'symbol']},
    {'pk': 'classifications', 'exists': True,
     'members': ['clinvar_highest_pathogenicity', 'clinvar_review_status', 'clinical_significance',
                 'conflicting_clinical_significance', 'clinvar_variation_id', 'clinvar_allele_id',
                 'clinvar_origin', 'clinvar_preferred_disease_name', 'clinvar_disease_database_name',
                 'clinvar_clinical_sources', 'drug_response', 'clinvar_highest_oncogenicity',
                 'clinvar_oncogenic_classification', 'clinvar_oncogenic_review_status',
                 'clinvar_somatic_tier', 'clinvar_somatic_review_status', 'max_internal_classification',
                 'internally_classified', 'internally_classified_labs',
                 'max_internal_somatic_classification', 'internally_classified_somatic'],
     'detail_only': ['conflicting_clinical_significance', 'clinvar_variation_id', 'clinvar_allele_id',
                     'clinvar_preferred_disease_name', 'clinvar_disease_database_name',
                     'clinvar_clinical_sources', 'internally_classified_labs']},
    {'pk': 'consequence_impact', 'level': _VARIANT, 'width': 160, 'label': 'Consequence',
     'description': 'Sequence Ontology consequence, with a coloured dot for the impact '
                    '(MODIFIER/LOW/MODERATE/HIGH). Sorts by consequence',
     'members': ['consequence', 'impact']},
    {'pk': 'gnomad', 'level': _VARIANT, 'width': 130, 'label': 'gnomAD',
     'description': 'Allele Frequency (0-1) among all gnomAD genotypes (exome+genome), the popmax '
                    'frequency and its population, and the gnomAD Pass/Fail link. Hover for the '
                    'per-population frequencies, allele counts, homozygotes and FAFs',
     'members': ['gnomad_af', 'gnomad_popmax', 'gnomad_filtered', 'gnomad_popmax_af', 'gnomad_ac',
                 'gnomad_an', 'gnomad_hom_alt', 'gnomad_hemi_count', 'gnomad_popmax_ac',
                 'gnomad_popmax_an', 'gnomad_popmax_hom_alt', 'gnomad_afr_af', 'gnomad_amr_af',
                 'gnomad_asj_af', 'gnomad_eas_af', 'gnomad_fin_af', 'gnomad_mid_af', 'gnomad_nfe_af',
                 'gnomad_oth_af', 'gnomad_sas_af', 'gnomad_xy_af', 'gnomad_xy_ac', 'gnomad_xy_an',
                 'gnomad_non_par', 'gnomad_faf95', 'gnomad_faf99', 'gnomad_fafmax_faf95_max',
                 'gnomad_fafmax_faf99_max', 'gnomad2_liftover_af'],
     'detail_only': ['gnomad_popmax', 'gnomad_filtered', 'gnomad_non_par']},
    {'pk': 'spliceai', 'level': _VARIANT, 'width': 80, 'label': 'SpliceAI',
     'description': "Highest of the four <a href='https://pubmed.ncbi.nlm.nih.gov/30661751/'>SpliceAI</a> "
                    "delta scores, with a dot at the 0.2 / 0.5 / 0.8 cutoffs. Hover for each prediction "
                    "and its position offset",
     'members': ['spliceai_max_ds', 'spliceai_pred_ds_ag', 'spliceai_pred_ds_al', 'spliceai_pred_ds_dg',
                 'spliceai_pred_ds_dl', 'spliceai_pred_dp_ag', 'spliceai_pred_dp_al',
                 'spliceai_pred_dp_dg', 'spliceai_pred_dp_dl', 'spliceai_gene_symbol'],
     'detail_only': ['spliceai_pred_dp_ag', 'spliceai_pred_dp_al', 'spliceai_pred_dp_dg',
                     'spliceai_pred_dp_dl', 'spliceai_gene_symbol']},
    {'pk': 'maxentscan', 'level': _VARIANT, 'width': 90, 'label': 'MaxEntScan',
     'description': "<a href='http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html'>MaxEntScan</a> "
                    "change in the 5' splice site score as a percentage of the reference. Negative means a "
                    "weakened site. Hover for the reference and alternate scores",
     'members': ['maxentscan_percent_diff_ref', 'maxentscan_ref', 'maxentscan_alt', 'maxentscan_diff']},
    {'pk': 'mastermind', 'level': _VARIANT, 'width': 80, 'label': 'Mastermind',
     'description': "<a href='https://mastermind.genomenon.com/'>Mastermind</a> article counts at "
                    "increasing granularity - cDNA, cDNA + protein, amino acid change. Click to open the "
                    "article list. Sorts by the cDNA count",
     'members': ['mastermind_count_1_cdna', 'mastermind_count_2_cdna_prot',
                 'mastermind_count_3_aa_change', 'mastermind_mmid3'],
     'detail_only': ['mastermind_mmid3']},
    {'pk': 'aloft', 'level': _VARIANT, 'width': 125, 'label': 'ALoFT',
     'description': "<a href='https://www.nature.com/articles/s41467-017-00443-5'>ALoFT</a> classification "
                    "of a loss of function variant - tolerant, recessive or dominant - with the probability "
                    "behind it. A call ALoFT is not confident in (p &gt;= 0.05) is greyed. Hover for all "
                    "three probabilities and the transcript it chose",
     'members': ['aloft_pred', 'aloft_prob_dominant', 'aloft_prob_recessive', 'aloft_prob_tolerant',
                 'aloft_high_confidence', 'aloft_ensembl_transcript'],
     'detail_only': ['aloft_ensembl_transcript']},
    {'pk': 'predictions', 'level': _VARIANT, 'width': 80, 'label': 'Predictions',
     'description': 'Pathogenicity prediction tools that made a call, damaging first - one segment '
                    'each. Hover for the call each tool made. Sorts by the number of damaging calls',
     'members': ['predictions_num_pathogenic', 'predictions_num_benign', 'sift',
                 'polyphen2_hvar_pred_most_damaging', 'mutation_taster_pred_most_damaging',
                 'mutation_assessor_pred_most_damaging', 'fathmm_pred_most_damaging',
                 'metalr_rankscore'],
     'detail_only': ['sift', 'polyphen2_hvar_pred_most_damaging', 'mutation_taster_pred_most_damaging',
                     'mutation_assessor_pred_most_damaging', 'fathmm_pred_most_damaging',
                     'metalr_rankscore']},
    {'pk': 'db_zygosity', 'level': _DATABASE, 'width': 70, 'label': 'Database counts',
     'description': 'Times this variant has been seen in this database - hom &#146; het. Hover for the '
                    'reference and unknown counts. Sorts by het count',
     'members': ['total_db_het', 'total_db_hom', 'total_db_ref', 'total_db_unk']},

    {'pk': 'alphamissense', 'level': _VARIANT, 'width': None, 'label': 'AlphaMissense',
     'description': 'AlphaMissense pathogenicity call, with the score and rankscore behind it on hover',
     'members': ['alphamissense_pred', 'alphamissense_score', 'alphamissense_rankscore']},
    {'pk': 'bayesdel', 'level': _VARIANT, 'width': None, 'label': 'BayesDel',
     'description': 'BayesDel (no allele frequency) deleteriousness score, with its rankscore on hover',
     'members': ['bayesdel_noaf_score', 'bayesdel_noaf_rankscore']},
    {'pk': 'cadd', 'level': _VARIANT, 'width': None, 'label': 'CADD',
     'description': 'CADD phred scaled score, with the raw score and rankscore on hover',
     'members': ['cadd_phred', 'cadd_raw', 'cadd_raw_rankscore']},
    {'pk': 'clinpred', 'level': _VARIANT, 'width': None, 'label': 'ClinPred',
     'description': 'ClinPred call, with the score and rankscore behind it on hover',
     'members': ['clinpred_pred', 'clinpred_score', 'clinpred_rankscore']},
    {'pk': 'eve', 'level': _VARIANT, 'width': None, 'label': 'EVE',
     'description': 'EVE class, with the score behind it on hover',
     'members': ['eve_class', 'eve_score']},
    {'pk': 'metarnn', 'level': _VARIANT, 'width': None, 'label': 'MetaRNN',
     'description': 'MetaRNN call, with the score behind it on hover',
     'members': ['metarnn_pred', 'metarnn_score']},
    {'pk': 'primateai', 'level': _VARIANT, 'width': None, 'label': 'PrimateAI',
     'description': 'PrimateAI call, with the score behind it on hover',
     'members': ['primateai_pred', 'primateai_score']},
    {'pk': 'revel', 'level': _VARIANT, 'width': None, 'label': 'REVEL',
     'description': 'REVEL score, with its rankscore on hover',
     'members': ['revel_score', 'revel_rankscore']},
    {'pk': 'vest4', 'level': _VARIANT, 'width': None, 'label': 'VEST4',
     'description': 'VEST4 score, with its rankscore on hover',
     'members': ['vest4_score', 'vest4_rankscore']},
    {'pk': 'varity', 'level': _VARIANT, 'width': None, 'label': 'VARITY',
     'description': 'VARITY R score, with the ER score on hover',
     'members': ['varity_r_score', 'varity_er_score']},
    {'pk': 'mutpred2', 'level': _VARIANT, 'width': None, 'label': 'MutPred2',
     'description': 'MutPred2 score, with the top 5 molecular mechanisms on hover',
     'members': ['mutpred2_score', 'mutpred2_top5_mechanisms'],
     'detail_only': ['mutpred2_top5_mechanisms']},
    {'pk': 'hipred', 'level': _GENE, 'width': None, 'label': 'HIPred',
     'description': 'HIPred haploinsufficiency call for the gene, with the probability on hover',
     'members': ['hipred_prediction', 'hipred_score']},
    {'pk': 'gene_damage_index', 'level': _GENE, 'width': None, 'label': 'Gene damage index',
     'description': 'Gene Damage Index phred score, with the raw score on hover',
     'members': ['gene_damage_index_phred', 'gene_damage_index_score']},
    {'pk': 'gene_indispensability', 'level': _GENE, 'width': None, 'label': 'Gene indispensability',
     'description': 'Gene indispensability call, with the score behind it on hover',
     'members': ['gene_indispensability_pred', 'gene_indispensability_score']},
    {'pk': 'gnomad_constraint', 'level': _GENE, 'width': None, 'label': 'gnomAD constraint',
     'description': "gnomAD probability the gene is loss of function intolerant (pLI). Hover for "
                    "observed/expected LoF and the tolerant / recessive probabilities",
     'members': ['gnomad_pli', 'gnomad_oe_lof', 'gnomad_pnull', 'gnomad_prec']},
    {'pk': 'gnomad_sv', 'level': _VARIANT, 'width': None, 'label': 'gnomAD SV overlap',
     'description': 'Allele frequency of the overlapping gnomAD structural variant. Hover for how much '
                    'it overlaps, its name and its coordinates',
     'members': ['gnomad_sv_overlap_af', 'gnomad_sv_overlap_percent', 'gnomad_sv_overlap_name',
                 'gnomad_sv_overlap_coords'],
     'detail_only': ['gnomad_sv_overlap_name', 'gnomad_sv_overlap_coords']},
    {'pk': 'conservation', 'level': _VARIANT, 'width': None, 'label': 'Conservation',
     'description': 'PhyloP 100 way vertebrate conservation, with the other PhyloP / PhastCons '
                    'alignments and GERP++ RS on hover',
     'members': ['phylop_100_way_vertebrate', 'phylop_30_way_mammalian', 'phylop_46_way_mammalian',
                 'phastcons_100_way_vertebrate', 'phastcons_30_way_mammalian',
                 'phastcons_46_way_mammalian', 'gerp_pp_rs']},
    {'pk': 'cosmic', 'level': _VARIANT, 'width': None, 'label': 'COSMIC',
     'description': 'COSMIC identifier, linked to the COSMIC entry. Hover for the sample count and the '
                    'legacy identifier',
     'members': ['cosmic_id', 'cosmic_count', 'cosmic_legacy_id'],
     'detail_only': ['cosmic_id', 'cosmic_legacy_id']},
    {'pk': 'denovo_db', 'level': _VARIANT, 'width': None, 'label': 'denovo-db',
     'description': 'denovo-db case count, with the control count, the studies, their primary '
                    'phenotypes and PubMed IDs on hover',
     'members': ['denovo_db_case_count', 'denovo_db_control_count', 'denovo_db_studies',
                 'denovo_db_primary_phenotypes', 'denovo_db_pubmed_ids'],
     'detail_only': ['denovo_db_studies', 'denovo_db_primary_phenotypes', 'denovo_db_pubmed_ids']},
    {'pk': 'open_targets', 'level': _VARIANT, 'width': None, 'label': 'Open Targets',
     'description': "<a href='https://platform.opentargets.org/'>Open Targets</a> GWAS locus-to-gene "
                    "score. Hover for the diseases, the genes, the QTL biosample and the study",
     'members': ['open_targets_gwas_l2g_score', 'open_targets_gwas_diseases',
                 'open_targets_gwas_gene_id', 'open_targets_qtl_biosample', 'open_targets_qtl_gene_id',
                 'open_targets_study_id', 'open_targets_study_type', 'open_targets_variant_id',
                 'open_targets_gwas_l2g_scores'],
     'detail_only': ['open_targets_gwas_diseases', 'open_targets_gwas_gene_id',
                     'open_targets_qtl_biosample', 'open_targets_qtl_gene_id', 'open_targets_study_id',
                     'open_targets_study_type', 'open_targets_variant_id',
                     'open_targets_gwas_l2g_scores']},
    {'pk': 'promoter_ai', 'level': _VARIANT, 'width': None, 'label': 'PromoterAI',
     'description': 'PromoterAI score, with the transcription start site position on hover',
     'members': ['promoter_ai_score', 'promoter_ai_tss_pos']},
    {'pk': 'mavedb', 'level': _VARIANT, 'width': None, 'label': 'MaveDB',
     'description': 'MaveDB functional assay score, with its URN (linked to the experiment set) on hover',
     'members': ['mavedb_score', 'mavedb_urn'],
     'detail_only': ['mavedb_urn']},
    {'pk': 'dbscsnv', 'level': _VARIANT, 'width': None, 'label': 'dbscSNV',
     'description': 'dbscSNV ada boost splicing score, with the random forest score on hover',
     'members': ['dbscsnv_ada_score', 'dbscsnv_rf_score']},
    {'pk': 'nmd_ptc', 'level': _TRANSCRIPT, 'width': None, 'label': 'NMD / PTC',
     'description': 'Nonsense mediated decay escape status. Hover for the escape call and the premature '
                    'termination codon distances',
     'members': ['nmd_escape_status', 'nmd_escaping_variant', 'ptc_distance_codons',
                 'ptc_last_junction_distance']},
    {'pk': 'protvar', 'level': _VARIANT, 'width': None, 'label': 'ProtVar',
     'description': 'ProtVar ddG stability, with the interface and pocket calls on hover',
     'members': ['protvar_stability', 'protvar_int', 'protvar_pocket']},
    {'pk': 'essential_gene', 'level': _GENE, 'width': None, 'label': 'Essential gene',
     'description': 'Essential gene call from CRISPR screens, with CRISPR2 and gene trap on hover',
     'members': ['essential_gene_crispr', 'essential_gene_crispr2', 'essential_gene_gene_trap']},
    {'pk': 'annotsv_acmg', 'level': _VARIANT, 'width': None, 'label': 'AnnotSV ACMG',
     'description': 'AnnotSV ACMG class for the structural variant, with the score and the criteria it '
                    'applied on hover',
     'members': ['annotsv_acmg_class', 'annotsv_acmg_score', 'annotsv_acmg_criteria'],
     'detail_only': ['annotsv_acmg_criteria']},
    {'pk': 'annotsv_benign_af', 'level': _VARIANT, 'width': None, 'label': 'AnnotSV benign AF',
     'description': 'Highest allele frequency of an overlapping benign loss, with gain, insertion and '
                    'inversion frequencies on hover',
     'members': ['annotsv_b_loss_af_max', 'annotsv_b_gain_af_max', 'annotsv_b_ins_af_max',
                 'annotsv_b_inv_af_max']},
    {'pk': 'annotsv_splice', 'level': _VARIANT, 'width': None, 'label': 'AnnotSV splice site',
     'description': 'Distance to the nearest splice site, with the site type on hover',
     'members': ['annotsv_dist_nearest_ss', 'annotsv_nearest_ss_type'],
     'detail_only': ['annotsv_nearest_ss_type']},
    {'pk': 'annotsv_omim', 'level': _VARIANT, 'width': None, 'label': 'AnnotSV OMIM',
     'description': 'OMIM morbid gene flag, with the OMIM ID, inheritance and phenotype on hover',
     'members': ['annotsv_omim_morbid', 'annotsv_omim_id', 'annotsv_omim_inheritance',
                 'annotsv_omim_phenotype'],
     'detail_only': ['annotsv_omim_id', 'annotsv_omim_inheritance', 'annotsv_omim_phenotype']},
    {'pk': 'annotsv_region', 'level': _VARIANT, 'width': None, 'label': 'AnnotSV region',
     'description': 'Repeat type at the left breakpoint. Hover for the right breakpoint, the segmental '
                    'duplications and the ENCODE blacklist regions on either side',
     'members': ['annotsv_repeat_type_left', 'annotsv_repeat_type_right', 'annotsv_segdup_left',
                 'annotsv_segdup_right', 'annotsv_encode_blacklist_left',
                 'annotsv_encode_blacklist_right', 'annotsv_encode_blacklist_characteristics_left',
                 'annotsv_encode_blacklist_characteristics_right'],
     'detail_only': ['annotsv_repeat_type_left', 'annotsv_repeat_type_right', 'annotsv_segdup_left',
                     'annotsv_segdup_right', 'annotsv_encode_blacklist_left',
                     'annotsv_encode_blacklist_right', 'annotsv_encode_blacklist_characteristics_left',
                     'annotsv_encode_blacklist_characteristics_right']},
]

_NEW_COMPOSITE_COLUMNS = [
    {'grid_column_name': c['pk'], 'variant_column': c['pk'], 'annotation_level': c['level'],
     'width': c['width'], 'label': c['label'], 'description': c['description'],
     'model_field': False, 'queryset_field': False}
    for c in _COMPOSITES if not c.get('exists')
]


def _add_composite_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", _NEW_COMPOSITE_COLUMNS)])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    for column_id, (width, label, description) in _RESTORED_MEMBERS.items():
        VariantGridColumn.objects.filter(pk=column_id).update(width=width, label=label,
                                                              description=description)
    for column_id, annotation_level in _FIXED_ANNOTATION_LEVELS.items():
        VariantGridColumn.objects.filter(pk=column_id).update(annotation_level=annotation_level)


def _add_members(apps, _schema_editor):
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    members = []
    for composite in _COMPOSITES:
        detail_only = set(composite.get('detail_only', []))
        for sort_order, column_id in enumerate(composite['members']):
            members.append(CompositeColumnMember(composite_id=composite['pk'], column_id=column_id,
                                                 sort_order=sort_order,
                                                 in_sort_menu=column_id not in detail_only))
    CompositeColumnMember.objects.bulk_create(members)


def _collapse_collections(apps, _schema_editor):
    for composite in _COMPOSITES:
        collapse_into_composite(apps, composite['pk'])


def _remove_composite_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    CompositeColumnMember.objects.all().delete()
    grid_names = [c['grid_column_name'] for c in _NEW_COMPOSITE_COLUMNS]
    VariantGridColumn.objects.filter(pk__in=grid_names).delete()  # cascades to CustomColumn


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0237_compositecolumnmember'),
    ]

    operations = [
        migrations.RunPython(_add_composite_columns, reverse_code=_remove_composite_columns),
        migrations.RunPython(_add_members, reverse_code=migrations.RunPython.noop),
        migrations.RunPython(_collapse_collections, reverse_code=migrations.RunPython.noop),
    ]
