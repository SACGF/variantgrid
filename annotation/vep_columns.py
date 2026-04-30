from dataclasses import dataclass
from typing import Optional, Union

from annotation.models.models_enums import (
    ColumnAnnotationCategory,
    VariantAnnotationPipelineType,
    VEPCustom,
    VEPPlugin,
)
from annotation.vep_config import VEPConfig


@dataclass(frozen=True)
class VEPColumnDef:
    """ Maps a VEP CSQ / plugin / custom output field -> one or more VariantGridColumn destinations.
        Empty `genome_builds` / `pipeline_types` mean "applies to all". """

    source_field: Optional[str]
    variant_grid_columns: tuple[str, ...]
    category: ColumnAnnotationCategory
    vep_plugin: Optional[VEPPlugin] = None
    vep_custom: Optional[VEPCustom] = None
    source_field_has_custom_prefix: bool = False
    genome_builds: frozenset[str] = frozenset()
    pipeline_types: frozenset[VariantAnnotationPipelineType] = frozenset()
    min_columns_version: Optional[int] = None
    max_columns_version: Optional[int] = None
    min_vep_version: Optional[int] = None
    max_vep_version: Optional[int] = None
    gnomad4_minor_version: Optional[str] = None
    summary_stats: Optional[str] = None
    source_field_processing_description: Optional[str] = None

    @property
    def vep_info_field(self) -> Optional[str]:
        """ Replicates ColumnVEPField.vep_info_field. """
        vif = self.source_field
        if self.vep_custom and self.source_field_has_custom_prefix:
            prefix = self.vep_custom.label
            vif = f"{prefix}_{vif}" if vif else prefix
        return vif

    @property
    def columns_version_description(self) -> str:
        limits = []
        if self.min_columns_version:
            limits.append(f"column version >= {self.min_columns_version}")
        if self.max_columns_version:
            limits.append(f"column version <= {self.max_columns_version}")
        return " and ".join(limits)

    def applies_to(
        self,
        *,
        genome_build_name: Optional[str] = None,
        pipeline_type: Optional[VariantAnnotationPipelineType] = None,
        columns_version: Optional[int] = None,
        vep_version: Optional[int] = None,
        gnomad4_minor_version: Optional[str] = None,
    ) -> bool:
        if genome_build_name is not None and self.genome_builds and genome_build_name not in self.genome_builds:
            return False
        if pipeline_type is not None and self.pipeline_types and pipeline_type not in self.pipeline_types:
            return False
        if columns_version is not None:
            if self.min_columns_version is not None and columns_version < self.min_columns_version:
                return False
            if self.max_columns_version is not None and columns_version > self.max_columns_version:
                return False
        if vep_version is not None:
            if self.min_vep_version is not None and vep_version < self.min_vep_version:
                return False
            if self.max_vep_version is not None and vep_version > self.max_vep_version:
                return False
        if gnomad4_minor_version is not None and self.gnomad4_minor_version is not None \
                and gnomad4_minor_version != self.gnomad4_minor_version:
            return False
        return True


# --- Family helpers -----------------------------------------------------------
# Each helper is a thin wrapper around VEPColumnDef(...) that pre-fills the
# defaults shared by a family of rows. Callers can override any default by
# passing it as a keyword argument (the override wins via dict.update).

GRCH37 = frozenset({'GRCh37'})
GRCH38 = frozenset({'GRCh38'})
GRCH37_38 = frozenset({'GRCh37', 'GRCh38'})
GRCH38_T2T = frozenset({'GRCh38', 'T2T-CHM13v2.0'})
GRCH37_38_T2T = frozenset({'GRCh37', 'GRCh38', 'T2T-CHM13v2.0'})
STANDARD = frozenset({VariantAnnotationPipelineType.STANDARD})
STRUCTURAL = frozenset({VariantAnnotationPipelineType.STRUCTURAL_VARIANT})

_ALOFT_DESC = 'Most damaging transcript prediction chosen, and Ensembl transcript stored.'
_GNOMAD2_AF_DESC = '(exome_AC+genome_AC)/(exome_AN+genome_AN)'
_INDEL_MAX_DESC = 'max() for indels'


def _to_vgc_tuple(vgc: Union[str, tuple[str, ...]]) -> tuple[str, ...]:
    return (vgc,) if isinstance(vgc, str) else tuple(vgc)


def _dbnsfp_v2(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ DBNSFP plugin, GRCh37/38, columns_version >= 2, pathogenicity bucket. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        "vep_plugin": VEPPlugin.DBNSFP,
        "genome_builds": GRCH37_38,
        "pipeline_types": STANDARD,
        "min_columns_version": 2,
        **overrides,
    })


def _dbnsfp_v1_most_damaging(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ DBNSFP plugin, columns_version == 1 (max), 'Most damaging' processing. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        "vep_plugin": VEPPlugin.DBNSFP,
        "pipeline_types": STANDARD,
        "max_columns_version": 1,
        "source_field_processing_description": 'Most damaging',
        **overrides,
    })


def _aloft(source_field: str, vgc, **overrides) -> VEPColumnDef:
    return _dbnsfp_v2(source_field, vgc, source_field_processing_description=_ALOFT_DESC, **overrides)


def _dbnsfp_v4(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ DBNSFP plugin, GRCh37/38/T2T, columns_version >= 4 (dbNSFP 5.3.1a),
        pathogenicity bucket. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        "vep_plugin": VEPPlugin.DBNSFP,
        "genome_builds": GRCH37_38_T2T,
        "pipeline_types": STANDARD,
        "min_columns_version": 4,
        **overrides,
    })


def _aloft_v4(source_field: str, vgc) -> VEPColumnDef:
    return _dbnsfp_v4(source_field, vgc, source_field_processing_description=_ALOFT_DESC)


def _gnomad4(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ gnomAD v4, custom prefix, GRCh38 + T2T, columns_version >= 3. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.FREQUENCY_DATA,
        "vep_custom": VEPCustom.GNOMAD_4,
        "source_field_has_custom_prefix": True,
        "genome_builds": GRCH38_T2T,
        "pipeline_types": STANDARD,
        "min_columns_version": 3,
        **overrides,
    })


def _gnomad3(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ gnomAD v3, custom prefix, GRCh38, columns_version <= 2. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.FREQUENCY_DATA,
        "vep_custom": VEPCustom.GNOMAD_3,
        "source_field_has_custom_prefix": True,
        "genome_builds": GRCH38,
        "pipeline_types": STANDARD,
        "max_columns_version": 2,
        **overrides,
    })


def _gnomad2_grch37(source_field: str, vgc, *, desc: Optional[str] = None,
                    **overrides) -> VEPColumnDef:
    """ gnomAD v2, custom prefix, GRCh37 only. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.FREQUENCY_DATA,
        "vep_custom": VEPCustom.GNOMAD_2,
        "source_field_has_custom_prefix": True,
        "genome_builds": GRCH37,
        "pipeline_types": STANDARD,
        "source_field_processing_description": desc,
        **overrides,
    })


def _gnomad_sv(source_field: str, vgc, **overrides) -> VEPColumnDef:
    """ gnomAD SV, custom prefix, structural-variant pipeline. genome_builds varies. """
    return VEPColumnDef(**{
        "source_field": source_field,
        "variant_grid_columns": _to_vgc_tuple(vgc),
        "category": ColumnAnnotationCategory.FREQUENCY_DATA,
        "vep_custom": VEPCustom.GNOMAD_SV,
        "source_field_has_custom_prefix": True,
        "pipeline_types": STRUCTURAL,
        **overrides,
    })


VEP_COLUMNS: tuple[VEPColumnDef, ...] = (
    # ---------- Aloft (DBNSFP, GRCh37/38, v2-3) ----------------------------
    _aloft('Ensembl_transcriptid', 'aloft_ensembl_transcript', max_columns_version=3),
    _aloft('Aloft_Confidence',     'aloft_high_confidence',    max_columns_version=3),
    _aloft('Aloft_pred',           'aloft_pred',               max_columns_version=3),
    _aloft('Aloft_prob_Dominant',  'aloft_prob_dominant',      max_columns_version=3),
    _aloft('Aloft_prob_Recessive', 'aloft_prob_recessive',     max_columns_version=3),
    _aloft('Aloft_prob_Tolerant',  'aloft_prob_tolerant',      max_columns_version=3),

    # ---------- DBNSFP rankscore family (v2-3, GRCh37/38) ------------------
    _dbnsfp_v2('AlphaMissense_rankscore', 'alphamissense_rankscore', min_columns_version=3, max_columns_version=3),
    _dbnsfp_v2('BayesDel_noAF_rankscore', 'bayesdel_noaf_rankscore', max_columns_version=3),
    _dbnsfp_v2('CADD_raw_rankscore',      'cadd_raw_rankscore',      max_columns_version=3),
    _dbnsfp_v2('ClinPred_rankscore',      'clinpred_rankscore',      max_columns_version=3),
    _dbnsfp_v2('MetaLR_rankscore',        'metalr_rankscore',        max_columns_version=3),
    _dbnsfp_v2('REVEL_rankscore',         'revel_rankscore',         max_columns_version=3),
    _dbnsfp_v2('VEST4_rankscore',         'vest4_rankscore',         max_columns_version=3),

    # ---------- DBNSFP "Most damaging" (columns_version == 1) --------------
    _dbnsfp_v1_most_damaging('FATHMM_pred',           'fathmm_pred_most_damaging'),
    _dbnsfp_v1_most_damaging('MutationAssessor_pred', 'mutation_assessor_pred_most_damaging'),
    _dbnsfp_v1_most_damaging('MutationTaster_pred',   'mutation_taster_pred_most_damaging'),
    _dbnsfp_v1_most_damaging('Polyphen2_HVAR_pred',   'polyphen2_hvar_pred_most_damaging'),

    # ---------- DBNSFP misc (no helper fits) -------------------------------
    VEPColumnDef(
        source_field='CADD_phred',
        variant_grid_columns=('cadd_phred',),
        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        vep_plugin=VEPPlugin.DBNSFP,
        pipeline_types=STANDARD,
        max_columns_version=1,
    ),
    VEPColumnDef(
        source_field='REVEL_score',
        variant_grid_columns=('revel_score',),
        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        vep_plugin=VEPPlugin.DBNSFP,
        pipeline_types=STANDARD,
        max_columns_version=1,
    ),
    VEPColumnDef(
        source_field='GERP++_RS',
        variant_grid_columns=('gerp_pp_rs',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_plugin=VEPPlugin.DBNSFP,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
        max_columns_version=3,
    ),
    VEPColumnDef(
        source_field='Interpro_domain',
        variant_grid_columns=('interpro_domain',),
        category=ColumnAnnotationCategory.PROTEIN_DOMAINS,
        vep_plugin=VEPPlugin.DBNSFP,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
        max_columns_version=3,
    ),

    # ---------- DBNSFP v4 (dbNSFP 5.3.1a, GRCh37/38/T2T, columns_version >= 4) ----------
    # ALoFT
    _aloft_v4('Ensembl_transcriptid', 'aloft_ensembl_transcript'),
    _aloft_v4('Aloft_Confidence',     'aloft_high_confidence'),
    _aloft_v4('Aloft_pred',           'aloft_pred'),
    _aloft_v4('Aloft_prob_Dominant',  'aloft_prob_dominant'),
    _aloft_v4('Aloft_prob_Recessive', 'aloft_prob_recessive'),
    _aloft_v4('Aloft_prob_Tolerant',  'aloft_prob_tolerant'),

    # Rankscores carried forward
    _dbnsfp_v4('AlphaMissense_rankscore', 'alphamissense_rankscore'),
    _dbnsfp_v4('BayesDel_noAF_rankscore', 'bayesdel_noaf_rankscore'),
    _dbnsfp_v4('CADD_raw_rankscore',      'cadd_raw_rankscore'),
    _dbnsfp_v4('ClinPred_rankscore',      'clinpred_rankscore'),
    _dbnsfp_v4('MetaLR_rankscore',        'metalr_rankscore'),
    _dbnsfp_v4('REVEL_rankscore',         'revel_rankscore'),
    _dbnsfp_v4('VEST4_rankscore',         'vest4_rankscore'),

    # Raw scores + predictions newly extracted at v4
    _dbnsfp_v4('AlphaMissense_pred',         'alphamissense_pred'),
    _dbnsfp_v4('AlphaMissense_score',        'alphamissense_score'),
    _dbnsfp_v4('BayesDel_noAF_score',        'bayesdel_noaf_score'),
    _dbnsfp_v4('CADD_phred',                 'cadd_phred'),
    _dbnsfp_v4('CADD_raw',                   'cadd_raw'),
    _dbnsfp_v4('ClinPred_pred',              'clinpred_pred'),
    _dbnsfp_v4('ClinPred_score',             'clinpred_score'),
    _dbnsfp_v4('MPC_score',                  'mpc_score'),
    _dbnsfp_v4('MetaRNN_pred',               'metarnn_pred'),
    _dbnsfp_v4('MetaRNN_score',              'metarnn_score'),
    _dbnsfp_v4('MutPred2_score',             'mutpred2_score'),
    _dbnsfp_v4('MutPred2_top5_mechanisms',   'mutpred2_top5_mechanisms'),
    _dbnsfp_v4('PrimateAI_pred',             'primateai_pred'),
    _dbnsfp_v4('PrimateAI_score',            'primateai_score'),
    _dbnsfp_v4('REVEL_score',                'revel_score'),
    _dbnsfp_v4('VARITY_ER_score',            'varity_er_score'),
    _dbnsfp_v4('VARITY_R_score',             'varity_r_score'),
    _dbnsfp_v4('VEST4_score',                'vest4_score'),

    # Conservation + protein domains carried forward across all 3 builds
    VEPColumnDef(
        source_field='GERP++_RS',
        variant_grid_columns=('gerp_pp_rs',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_plugin=VEPPlugin.DBNSFP,
        genome_builds=GRCH37_38_T2T,
        pipeline_types=STANDARD,
        min_columns_version=4,
    ),
    VEPColumnDef(
        source_field='Interpro_domain',
        variant_grid_columns=('interpro_domain',),
        category=ColumnAnnotationCategory.PROTEIN_DOMAINS,
        vep_plugin=VEPPlugin.DBNSFP,
        genome_builds=GRCH37_38_T2T,
        pipeline_types=STANDARD,
        min_columns_version=4,
    ),

    # ---------- gnomAD v4 (GRCh38 + T2T, columns_version >= 3) -------------
    _gnomad4('AC',                  'gnomad_ac'),
    _gnomad4('AF',                  'gnomad_af'),
    _gnomad4('AN',                  'gnomad_an'),
    _gnomad4('AF_afr',              'gnomad_afr_af'),
    _gnomad4('AF_amr',              'gnomad_amr_af'),
    _gnomad4('AF_asj',              'gnomad_asj_af'),
    _gnomad4('AF_eas',              'gnomad_eas_af'),
    _gnomad4('faf95',               'gnomad_faf95',           genome_builds=GRCH38),
    _gnomad4('faf99',               'gnomad_faf99',           genome_builds=GRCH38),
    _gnomad4('fafmax_faf95_max',    'gnomad_fafmax_faf95_max'),
    _gnomad4('fafmax_faf99_max',    'gnomad_fafmax_faf99_max'),
    _gnomad4('gnomad_filtered',     'gnomad_filtered',        gnomad4_minor_version='4.0'),
    _gnomad4('FILTER',              'gnomad_filtered',        gnomad4_minor_version='4.1'),
    _gnomad4('AF_fin',              'gnomad_fin_af'),
    _gnomad4('nhomalt',             'gnomad_hom_alt'),
    _gnomad4('AF_mid',              'gnomad_mid_af'),
    _gnomad4('AF_nfe',              'gnomad_nfe_af'),
    _gnomad4('non_par',             'gnomad_non_par',
             source_field_processing_description='nonpar from genomes'),
    _gnomad4('AF_remaining',        'gnomad_oth_af'),
    _gnomad4('grpmax',              'gnomad_popmax'),
    _gnomad4('AC_grpmax',           'gnomad_popmax_ac'),
    _gnomad4('AF_grpmax',           'gnomad_popmax_af'),
    _gnomad4('AN_grpmax',           'gnomad_popmax_an'),
    _gnomad4('AF_sas',              'gnomad_sas_af'),
    _gnomad4('AC_XY',               'gnomad_xy_ac'),
    _gnomad4('AF_XY',               'gnomad_xy_af'),
    _gnomad4('AN_XY',               'gnomad_xy_an'),

    # ---------- gnomAD v3 (GRCh38, columns_version <= 2) -------------------
    _gnomad3('AC',                  'gnomad_ac'),
    _gnomad3('AF',                  'gnomad_af'),
    _gnomad3('AN',                  'gnomad_an'),
    _gnomad3('AF-afr',              'gnomad_afr_af'),
    _gnomad3('AF-amr',              'gnomad_amr_af'),
    _gnomad3('AF-asj',              'gnomad_asj_af'),
    _gnomad3('AF-eas',              'gnomad_eas_af'),
    _gnomad3('FILTER',              'gnomad_filtered'),
    _gnomad3('AF-fin',              'gnomad_fin_af'),
    _gnomad3('nhomalt',             'gnomad_hom_alt'),
    _gnomad3('AF-nfe',              'gnomad_nfe_af'),
    _gnomad3('AF-oth',              'gnomad_oth_af'),
    _gnomad3('popmax',              'gnomad_popmax'),
    _gnomad3('AC_popmax',           'gnomad_popmax_ac'),
    _gnomad3('AF_popmax',           'gnomad_popmax_af'),
    _gnomad3('AN_popmax',           'gnomad_popmax_an'),
    _gnomad3('nhomalt_popmax',      'gnomad_popmax_hom_alt'),
    _gnomad3('AF-sas',              'gnomad_sas_af'),

    # ---------- gnomAD v2 (GRCh37) -----------------------------------------
    _gnomad2_grch37('AC',           'gnomad_ac',     desc='Sum of exome AC + genome AC'),
    _gnomad2_grch37('AF',           'gnomad_af',     desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('AN',           'gnomad_an',     desc='Sum of exome AN + genome AN'),
    _gnomad2_grch37('AF_afr',       'gnomad_afr_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('AF_amr',       'gnomad_amr_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('AF_asj',       'gnomad_asj_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('AF_eas',       'gnomad_eas_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('gnomad_filtered', 'gnomad_filtered',
                    desc='Filtered in either exome or genomes VCF'),
    _gnomad2_grch37('AF_fin',       'gnomad_fin_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('nhomalt',      'gnomad_hom_alt', desc='exome_nhomalt + genome_nhomalt'),
    _gnomad2_grch37('AF_nfe',       'gnomad_nfe_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('nonpar',       'gnomad_non_par', desc='nonpar from genomes'),
    _gnomad2_grch37('AF_oth',       'gnomad_oth_af'),
    _gnomad2_grch37('popmax',       'gnomad_popmax'),
    _gnomad2_grch37('AC_popmax',    'gnomad_popmax_ac',
                    desc='Sum of exome AC_popmax + genome AC_popmax'),
    _gnomad2_grch37('AF_popmax',    'gnomad_popmax_af',
                    desc='Re-calculated from max (exome_AC+genome_AC)/(exome_AN+genome_AN)'),
    _gnomad2_grch37('AN_popmax',    'gnomad_popmax_an',
                    desc='Sum of exome AN_popmax + genome AN_popmax'),
    _gnomad2_grch37('AF_sas',       'gnomad_sas_af', desc=_GNOMAD2_AF_DESC),
    _gnomad2_grch37('AC_male',      'gnomad_xy_ac',  min_columns_version=3),
    _gnomad2_grch37('AF_male',      'gnomad_xy_af',  min_columns_version=3),
    _gnomad2_grch37('AN_male',      'gnomad_xy_an',  min_columns_version=3),

    # ---------- gnomAD v2 liftover (GRCh38) --------------------------------
    VEPColumnDef(
        source_field='AF',
        variant_grid_columns=('gnomad2_liftover_af',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.GNOMAD_2,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH38,
        pipeline_types=STANDARD,
    ),

    # ---------- gnomAD SV (structural-variant pipeline) --------------------
    _gnomad_sv('AC',        'gnomad_ac'),
    _gnomad_sv('AF',        ('gnomad_sv_overlap_af', 'gnomad_af')),
    _gnomad_sv('AN',        'gnomad_an'),
    _gnomad_sv('N_HOMALT',  'gnomad_hom_alt'),
    _gnomad_sv('POPMAX_AF', 'gnomad_popmax_af'),
    _gnomad_sv('',          ('gnomad_sv_overlap_coords', 'gnomad_sv_overlap_percent')),
    _gnomad_sv('afr_AF',    'gnomad_afr_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('amr_AF',    'gnomad_amr_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('asj_AF',    'gnomad_asj_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('eas_AF',    'gnomad_eas_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('fin_AF',    'gnomad_fin_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('mid_AF',    'gnomad_mid_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('nfe_AF',    'gnomad_nfe_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('oth_AF',    'gnomad_oth_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('sas_AF',    'gnomad_sas_af', genome_builds=GRCH38_T2T),
    _gnomad_sv('AFR_AF',    'gnomad_afr_af', genome_builds=GRCH37),
    _gnomad_sv('AMR_AF',    'gnomad_amr_af', genome_builds=GRCH37),
    _gnomad_sv('EAS_AF',    'gnomad_eas_af', genome_builds=GRCH37),
    _gnomad_sv('EUR_AF',    'gnomad_nfe_af', genome_builds=GRCH37),
    _gnomad_sv('OTH_AF',    'gnomad_oth_af', genome_builds=GRCH37),

    # ---------- gnomAD SV name (different custom track) --------------------
    VEPColumnDef(
        source_field='gnomAD_SV_name',
        variant_grid_columns=('gnomad_sv_overlap_name',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.GNOMAD_SV_NAME,
        pipeline_types=STRUCTURAL,
    ),

    # ---------- SpliceAI ---------------------------------------------------
    VEPColumnDef(
        source_field='SpliceAI_pred_SYMBOL',
        variant_grid_columns=('spliceai_gene_symbol',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DP_AG',
        variant_grid_columns=('spliceai_pred_dp_ag',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DP_AL',
        variant_grid_columns=('spliceai_pred_dp_al',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DP_DG',
        variant_grid_columns=('spliceai_pred_dp_dg',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DP_DL',
        variant_grid_columns=('spliceai_pred_dp_dl',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DS_AG',
        variant_grid_columns=('spliceai_pred_ds_ag',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DS_AL',
        variant_grid_columns=('spliceai_pred_ds_al',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DS_DG',
        variant_grid_columns=('spliceai_pred_ds_dg',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='SpliceAI_pred_DS_DL',
        variant_grid_columns=('spliceai_pred_ds_dl',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEAI,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),

    # ---------- MaxEntScan -------------------------------------------------
    VEPColumnDef(
        source_field='MaxEntScan_alt',
        variant_grid_columns=('maxentscan_alt',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.MAXENTSCAN,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='MaxEntScan_diff',
        variant_grid_columns=('maxentscan_diff',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.MAXENTSCAN,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='MaxEntScan_ref',
        variant_grid_columns=('maxentscan_ref',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.MAXENTSCAN,
        pipeline_types=STANDARD,
    ),

    # ---------- Mastermind -------------------------------------------------
    VEPColumnDef(
        source_field='Mastermind_counts',
        variant_grid_columns=('mastermind_count_1_cdna',),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_plugin=VEPPlugin.MASTERMIND,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
        source_field_processing_description='MMCNT1',
    ),
    VEPColumnDef(
        source_field='Mastermind_counts',
        variant_grid_columns=('mastermind_count_2_cdna_prot',),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_plugin=VEPPlugin.MASTERMIND,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
        source_field_processing_description='MMCNT2',
    ),
    VEPColumnDef(
        source_field='Mastermind_counts',
        variant_grid_columns=('mastermind_count_3_aa_change',),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_plugin=VEPPlugin.MASTERMIND,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
        source_field_processing_description='MMCNT3',
    ),
    VEPColumnDef(
        source_field='Mastermind_MMID3',
        variant_grid_columns=('mastermind_mmid3',),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_plugin=VEPPlugin.MASTERMIND,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),

    # ---------- PhastCons / PhyloP (VEP-version split, max() for indels) ---
    VEPColumnDef(
        source_field='phastCons100way_vertebrate',
        variant_grid_columns=('phastcons_100_way_vertebrate',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHASTCONS_100_WAY,
        pipeline_types=STANDARD,
        max_vep_version=111,
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='max',
        variant_grid_columns=('phastcons_100_way_vertebrate',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHASTCONS_100_WAY,
        source_field_has_custom_prefix=True,
        min_vep_version=112,
        summary_stats='max',
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='phastCons30way_mammalian',
        variant_grid_columns=('phastcons_30_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHASTCONS_30_WAY,
        genome_builds=GRCH38,
        pipeline_types=STANDARD,
        max_vep_version=111,
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='max',
        variant_grid_columns=('phastcons_30_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHASTCONS_30_WAY,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH38,
        min_vep_version=112,
        summary_stats='max',
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='phastCons46way_mammalian',
        variant_grid_columns=('phastcons_46_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHASTCONS_46_WAY,
        genome_builds=GRCH37,
        pipeline_types=STANDARD,
        max_vep_version=111,
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='max',
        variant_grid_columns=('phastcons_46_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHASTCONS_46_WAY,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37,
        min_vep_version=112,
        summary_stats='max',
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='phyloP100way_vertebrate',
        variant_grid_columns=('phylop_100_way_vertebrate',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHYLOP_100_WAY,
        pipeline_types=STANDARD,
        max_vep_version=111,
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='max',
        variant_grid_columns=('phylop_100_way_vertebrate',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHYLOP_100_WAY,
        source_field_has_custom_prefix=True,
        min_vep_version=112,
        summary_stats='max',
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='phyloP30way_mammalian',
        variant_grid_columns=('phylop_30_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHYLOP_30_WAY,
        genome_builds=GRCH38,
        pipeline_types=STANDARD,
        max_vep_version=111,
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='max',
        variant_grid_columns=('phylop_30_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHYLOP_30_WAY,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH38,
        min_vep_version=112,
        summary_stats='max',
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='phyloP46way_mammalian',
        variant_grid_columns=('phylop_46_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHYLOP_46_WAY,
        genome_builds=GRCH37,
        pipeline_types=STANDARD,
        max_vep_version=111,
        source_field_processing_description=_INDEL_MAX_DESC,
    ),
    VEPColumnDef(
        source_field='max',
        variant_grid_columns=('phylop_46_way_mammalian',),
        category=ColumnAnnotationCategory.CONSERVATION,
        vep_custom=VEPCustom.PHYLOP_46_WAY,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37,
        min_vep_version=112,
        summary_stats='max',
        source_field_processing_description=_INDEL_MAX_DESC,
    ),

    # ---------- Plain VEP CSQ (no plugin / no custom) ----------------------
    VEPColumnDef(source_field='AF',              variant_grid_columns=('af_1kg',),        category=ColumnAnnotationCategory.FREQUENCY_DATA),
    VEPColumnDef(source_field='Amino_acids',     variant_grid_columns=('amino_acids',),   category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='CANONICAL',       variant_grid_columns=('canonical',),     category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='Codons',          variant_grid_columns=('codons',),        category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='Consequence',     variant_grid_columns=('consequence',),   category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(
        source_field='Existing_variation',
        variant_grid_columns=('cosmic_id',),
        category=ColumnAnnotationCategory.EXTERNAL_ID,
        source_field_processing_description='Extract COSMIC IDs',
    ),
    VEPColumnDef(
        source_field='Existing_variation',
        variant_grid_columns=('dbsnp_rs_id',),
        category=ColumnAnnotationCategory.EXTERNAL_ID,
        source_field_processing_description='Extract rsIDs',
    ),
    VEPColumnDef(source_field='DISTANCE',        variant_grid_columns=('distance',),      category=ColumnAnnotationCategory.NEARBY_FEATURES),
    VEPColumnDef(source_field='DOMAINS',         variant_grid_columns=('domains',),       category=ColumnAnnotationCategory.PROTEIN_DOMAINS),
    VEPColumnDef(source_field='ENSP',            variant_grid_columns=('ensembl_protein',), category=ColumnAnnotationCategory.EXTERNAL_ID),
    VEPColumnDef(source_field='EXON',            variant_grid_columns=('exon',),          category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='FLAGS',           variant_grid_columns=('flags',),         category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='HGNC_ID',         variant_grid_columns=('hgnc_id',),       category=ColumnAnnotationCategory.EXTERNAL_ID),
    VEPColumnDef(source_field='HGVSc',           variant_grid_columns=('hgvs_c',),        category=ColumnAnnotationCategory.HGVS),
    VEPColumnDef(source_field='HGVSp',           variant_grid_columns=('hgvs_p',),        category=ColumnAnnotationCategory.HGVS),
    VEPColumnDef(source_field='IMPACT',          variant_grid_columns=('impact',),        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS),
    VEPColumnDef(source_field='INTRON',          variant_grid_columns=('intron',),        category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='Protein_position', variant_grid_columns=('protein_position',), category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='PUBMED',          variant_grid_columns=('pubmed',),        category=ColumnAnnotationCategory.LITERATURE),
    VEPColumnDef(
        source_field='SIFT',
        variant_grid_columns=('sift',),
        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        genome_builds=GRCH37_38,
    ),
    VEPColumnDef(
        source_field='SOMATIC',
        variant_grid_columns=('somatic',),
        category=ColumnAnnotationCategory.EXTERNAL_ID,
        source_field_processing_description='Any values present',
    ),
    VEPColumnDef(source_field='SYMBOL',          variant_grid_columns=('symbol',),        category=ColumnAnnotationCategory.GENE_ANNOTATIONS),
    VEPColumnDef(source_field='VARIANT_CLASS',   variant_grid_columns=('variant_class',), category=ColumnAnnotationCategory.GENE_ANNOTATIONS),

    # ---------- Other custom tracks ----------------------------------------
    VEPColumnDef(
        source_field='AF',
        variant_grid_columns=('af_uk10k',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.UK10K,
        source_field_has_custom_prefix=True,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='CNT',
        variant_grid_columns=('cosmic_count',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.COSMIC,
        source_field_has_custom_prefix=True,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='LEGACY_ID',
        variant_grid_columns=('cosmic_legacy_id',),
        category=ColumnAnnotationCategory.EXTERNAL_ID,
        vep_custom=VEPCustom.COSMIC,
        source_field_has_custom_prefix=True,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='REPEAT_MASKER',
        variant_grid_columns=('repeat_masker',),
        category=ColumnAnnotationCategory.SEQUENCE,
        vep_custom=VEPCustom.REPEAT_MASKER,
    ),
    VEPColumnDef(
        source_field='TOPMED',
        variant_grid_columns=('topmed_af',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.TOPMED,
        source_field_has_custom_prefix=True,
        pipeline_types=STANDARD,
    ),

    # ---------- denovo-db (custom VCF, GRCh37/38) --------------------------
    VEPColumnDef(
        source_field='StudyName',
        variant_grid_columns=('denovo_db_studies',),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_custom=VEPCustom.DENOVO_DB,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='PubmedID',
        variant_grid_columns=('denovo_db_pubmed_ids',),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_custom=VEPCustom.DENOVO_DB,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='PrimaryPhenotype',
        variant_grid_columns=('denovo_db_primary_phenotypes',),
        category=ColumnAnnotationCategory.PHENOTYPE,
        vep_custom=VEPCustom.DENOVO_DB,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='CASE_CT',
        variant_grid_columns=('denovo_db_case_count',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.DENOVO_DB,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='CONTROL_CT',
        variant_grid_columns=('denovo_db_control_count',),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.DENOVO_DB,
        source_field_has_custom_prefix=True,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),

    # ---------- Other plugins ----------------------------------------------
    VEPColumnDef(
        source_field='ada_score',
        variant_grid_columns=('dbscsnv_ada_score',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.DBSCSNV,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='rf_score',
        variant_grid_columns=('dbscsnv_rf_score',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.DBSCSNV,
        genome_builds=GRCH37_38,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='Grantham',
        variant_grid_columns=('grantham',),
        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        vep_plugin=VEPPlugin.GRANTHAM,
        pipeline_types=STANDARD,
    ),
    VEPColumnDef(
        source_field='MaveDB_score',
        variant_grid_columns=('mavedb_score',),
        category=ColumnAnnotationCategory.FUNCTIONAL_EFFECT,
        vep_plugin=VEPPlugin.MAVEDB,
        genome_builds=GRCH38,
        pipeline_types=STANDARD,
        min_columns_version=3,
    ),
    VEPColumnDef(
        source_field='MaveDB_urn',
        variant_grid_columns=('mavedb_urn',),
        category=ColumnAnnotationCategory.FUNCTIONAL_EFFECT,
        vep_plugin=VEPPlugin.MAVEDB,
        genome_builds=GRCH38,
        pipeline_types=STANDARD,
        min_columns_version=3,
    ),
    VEPColumnDef(
        source_field='NMD',
        variant_grid_columns=('nmd_escaping_variant',),
        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        vep_plugin=VEPPlugin.NMD,
        pipeline_types=STANDARD,
        min_columns_version=2,
    ),
    VEPColumnDef(
        source_field='SpliceRegion',
        variant_grid_columns=('splice_region',),
        category=ColumnAnnotationCategory.SPLICING_PREDICTIONS,
        vep_plugin=VEPPlugin.SPLICEREGION,
        pipeline_types=STANDARD,
    ),
)


# vep_config setting keys (under settings.ANNOTATION[<build>]["vep_config"]) that must be
# present (non-None) for a VEPColumnDef to actually be populated. Empty tuple means the
# def needs no data file (plain VEP CSQ, or a plugin with no data — Grantham/SpliceRegion/
# NMD/LoFtool). VEPCustom keys are derived from `vep_custom.label.lower()`.
_VEP_PLUGIN_SETTING_KEYS: dict[VEPPlugin, tuple[str, ...]] = {
    VEPPlugin.DBNSFP: ("dbnsfp",),
    VEPPlugin.DBSCSNV: ("dbscsnv",),
    VEPPlugin.MASTERMIND: ("mastermind",),
    VEPPlugin.MAVEDB: ("mave",),
    VEPPlugin.MAXENTSCAN: ("maxentscan",),
    VEPPlugin.SPLICEAI: ("spliceai_snv", "spliceai_indel"),
    VEPPlugin.GRANTHAM: (),
    VEPPlugin.LOFTOOL: (),
    VEPPlugin.NMD: (),
    VEPPlugin.SPLICEREGION: (),
}


def required_data_keys(c: VEPColumnDef) -> tuple[str, ...]:
    """ vep_config setting keys this def needs in order to actually be populated.
        Empty tuple = no data file required. """
    if c.vep_custom:
        return (c.vep_custom.label.lower(),)
    if c.vep_plugin:
        return _VEP_PLUGIN_SETTING_KEYS.get(c.vep_plugin, ())
    return ()


def has_data_files(c: VEPColumnDef, vep_config: VEPConfig) -> bool:
    """ True iff every required vep_config key is set (non-None) for this def's build.
        VEPConfig.__getitem__ raises KeyError when the value is missing or None. """
    for key in required_data_keys(c):
        try:
            _ = vep_config[key]
        except KeyError:
            return False
    return True


def filter_for(
    *,
    vep_config: Optional[VEPConfig] = None,
    genome_build_name: Optional[str] = None,
    pipeline_type: Optional[VariantAnnotationPipelineType] = None,
    columns_version: Optional[int] = None,
    vep_version: Optional[int] = None,
    gnomad4_minor_version: Optional[str] = None,
    vep_plugin: Optional[VEPPlugin] = None,
    vep_custom: Optional[VEPCustom] = None,
) -> tuple[VEPColumnDef, ...]:
    """ When `vep_config` is provided, fills in build/columns/vep/gnomad versions from it
        and drops any def whose data file isn't configured for the build. """
    if vep_config is not None:
        if genome_build_name is None:
            genome_build_name = vep_config.genome_build.name
        if columns_version is None:
            columns_version = vep_config.columns_version
        if vep_version is None:
            vep_version = vep_config.vep_version
        if gnomad4_minor_version is None:
            gnomad4_minor_version = vep_config.gnomad4_minor_version

    return tuple(
        c for c in VEP_COLUMNS
        if c.applies_to(
            genome_build_name=genome_build_name,
            pipeline_type=pipeline_type,
            columns_version=columns_version,
            vep_version=vep_version,
            gnomad4_minor_version=gnomad4_minor_version,
        )
        and (vep_plugin is None or c.vep_plugin == vep_plugin)
        and (vep_custom is None or c.vep_custom == vep_custom)
        and (vep_config is None or has_data_files(c, vep_config))
    )


def source_fields_for(**kwargs) -> list[str]:
    """ Distinct source_fields for the given filter, sorted case-insensitively
        (matches postgres `ORDER BY source_field` under default en_US collation). """
    return sorted({c.source_field for c in filter_for(**kwargs) if c.source_field}, key=str.lower)


def for_variant_grid_column(vgc_id: str, *, vep_config: Optional[VEPConfig] = None) -> tuple[VEPColumnDef, ...]:
    """ When `vep_config` is provided, drops defs whose data file isn't configured. """
    return tuple(
        c for c in VEP_COLUMNS
        if vgc_id in c.variant_grid_columns
        and (vep_config is None or has_data_files(c, vep_config))
    )


def all_variant_grid_column_ids() -> frozenset[str]:
    return frozenset(vgc for c in VEP_COLUMNS for vgc in c.variant_grid_columns)
