from django.db import migrations

# sacgf/variantgrid#1130 - put back raw pathogenicity scores
# sacgf/variantgrid#951 - rankscore evidence keys (deferred - not adding rankscore ekeys for now)
# sacgf/variantgrid_shariant#228 - CNV support evidence keys (svlen, copy_number)


def _update_pathogenicity_ekeys(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    # 1) Restore raw scores - point existing ekeys back at raw score columns
    raw_score_mapping = {
        "cadd": "cadd_phred",
        "revel": "revel_score",
        "bayesdel": "bayesdel_noaf_score",
        "clinpred": "clinpred_score",
        "vest": "vest4_score",
    }
    for ek_pk, variantgrid_column_id in raw_score_mapping.items():
        EvidenceKey.objects.filter(pk=ek_pk).update(variantgrid_column_id=variantgrid_column_id)

    # 2) MetaLR replaced by MetaRNN - hide the old key
    EvidenceKey.objects.filter(pk="metalr").update(hide=True)

    # 3) New evidence keys
    FREE_ENTRY = 'F'
    FLOAT = 'L'
    INTEGER = 'I'

    CP = 'CP'  # COMPUTATIONAL_AND_PREDICTIVE_DATA
    V = 'V'    # VARIANT

    order_near_cadd = 5
    order_metarnn = 6  # group MetaRNN score+pred together

    new_ekeys = [
        # MetaRNN replaces MetaLR
        EvidenceKey(
            key='metarnn',
            label='MetaRNN',
            sub_label='Meta-predictor',
            description='MetaRNN raw score for missense pathogenicity. Higher = more likely pathogenic.',
            examples=[0.1, 0.89],
            options=[],
            see='https://github.com/Chang-Li-UCF/MetaRNN',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_metarnn,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='metarnn_score',
        ),
        EvidenceKey(
            key='metarnn_pred',
            label='MetaRNN pred',
            sub_label='Meta-predictor',
            description='MetaRNN categorical prediction.',
            examples=['Damaging', 'Tolerated'],
            options=[],
            see='https://github.com/Chang-Li-UCF/MetaRNN',
            evidence_category=CP,
            value_type=FREE_ENTRY,
            order=order_metarnn,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='metarnn_pred',
        ),

        # AlphaMissense raw score - complements existing alphamissense_rankscore
        EvidenceKey(
            key='alphamissense',
            label='AlphaMissense',
            description='AlphaMissense raw pathogenicity score. Higher = more likely pathogenic.',
            examples=[0.1, 0.89],
            options=[],
            see='https://www.science.org/doi/10.1126/science.adg7492',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_near_cadd,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='alphamissense_score',
        ),

        # Other dbNSFPv4 raw score predictors
        EvidenceKey(
            key='mpc',
            label='MPC',
            description='Missense badness, PolyPhen-2, and Constraint score. Higher = more deleterious.',
            examples=[0.5, 2.5],
            options=[],
            see='https://www.biorxiv.org/content/10.1101/148353v1',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_near_cadd,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='mpc_score',
        ),
        EvidenceKey(
            key='mutpred2',
            label='MutPred2',
            description='MutPred2 general pathogenicity score.',
            examples=[0.1, 0.89],
            options=[],
            see='http://mutpred.mutdb.org/',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_near_cadd,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='mutpred2_score',
        ),
        EvidenceKey(
            key='primateai',
            label='PrimateAI',
            description='PrimateAI raw score. Higher = more likely pathogenic.',
            examples=[0.1, 0.89],
            options=[],
            see='https://github.com/Illumina/PrimateAI',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_near_cadd,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='primateai_score',
        ),
        EvidenceKey(
            key='varity_er',
            label='VARITY_ER',
            description='VARITY_ER score (rare exomes training set) for missense pathogenicity.',
            examples=[0.1, 0.89],
            options=[],
            see='http://varity.varianteffect.org/',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_near_cadd,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='varity_er_score',
        ),
        EvidenceKey(
            key='varity_r',
            label='VARITY_R',
            description='VARITY_R score (rare training set) for missense pathogenicity.',
            examples=[0.1, 0.89],
            options=[],
            see='http://varity.varianteffect.org/',
            evidence_category=CP,
            value_type=FLOAT,
            order=order_near_cadd,
            mandatory=False,
            max_share_level='logged_in_users',
            copy_consensus=True,
            variantgrid_column_id='varity_r_score',
        ),

        # CNV support - shariant#228
        EvidenceKey(
            key='svlen',
            label='SVLEN',
            description='Length of structural variant (matches VCF SVLEN field).',
            examples=[1000, -500],
            options=[],
            evidence_category=V,
            value_type=INTEGER,
            order=12,
            mandatory=False,
            max_share_level='public',
            copy_consensus=False,
            variantgrid_column_id='svlen',
        ),
        EvidenceKey(
            key='copy_number',
            label='CN',
            description='Copy number (matches VCF CN field). Used to distinguish duplications from higher-copy gains.',
            examples=[0, 1, 3, 4],
            options=[],
            evidence_category=V,
            value_type=INTEGER,
            order=13,
            mandatory=False,
            max_share_level='public',
            copy_consensus=False,
            variantgrid_column_id=None,
        ),
    ]
    EvidenceKey.objects.bulk_create(new_ekeys)


def _reverse_update_pathogenicity_ekeys(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    # Reverse: re-point existing ekeys back at rankscore columns
    rankscore_mapping = {
        "cadd": "cadd_raw_rankscore",
        "revel": "revel_rankscore",
        "bayesdel": "bayesdel_noaf_rankscore",
        "clinpred": "clinpred_rankscore",
        "vest": "vest4_rankscore",
    }
    for ek_pk, variantgrid_column_id in rankscore_mapping.items():
        EvidenceKey.objects.filter(pk=ek_pk).update(variantgrid_column_id=variantgrid_column_id)

    EvidenceKey.objects.filter(pk="metalr").update(hide=False)

    EvidenceKey.objects.filter(pk__in=[
        'metarnn', 'metarnn_pred', 'alphamissense', 'mpc', 'mutpred2',
        'primateai', 'varity_er', 'varity_r', 'svlen', 'copy_number',
    ]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0164_somatic_hrd_msi_tmb_ekeys'),
        ('snpdb', '0176_dbnsfp_v4_pathogenicity_columns'),
    ]

    operations = [
        migrations.RunPython(_update_pathogenicity_ekeys, reverse_code=_reverse_update_pathogenicity_ekeys),
    ]
