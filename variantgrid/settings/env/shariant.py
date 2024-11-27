"""
Shariant - https://shariant.org.au

See https://github.com/sacgf/variantgrid/wiki/Annotation%20Setup

"""

# IMPORTANT : THE BELOW IMPORTS ARE USED TO APPLY THEIR RESPECTIVE SETTINGS VALUES
from variantgrid.settings.env.shariantcommon import *  # pylint: disable=wildcard-import, unused-wildcard-import

# import all the base settings #
SITE_ID = 5  # shariant.org.au
SITE_NAME = "Shariant"

# HEARTBEAT_URL = 'https://heartbeat.uptimerobot.com/m788641874-4c58c98a716180f36670e551a0bd03fff47abfea'
SEND_EMAILS = True
OIDC_RP_CLIENT_ID = 'shariant'
OIDC_REQUIRED_GROUP = '/variantgrid/shariant_production'
LOGIN_URL = '/oidc_login/'
LOGOUT_REDIRECT_URL = KEY_CLOAK_PROTOCOL_BASE + '/logout?redirect_uri=https%3A%2F%2Fshariant.org.au'

CLASSIFICATION_DOWNLOADABLE_NOTES_AND_EXPLAINS = False
CLASSIFICATION_DOWNLOADABLE_FIELDS = set([
    "acmg:ba1",
    "acmg:bp1",
    "acmg:bp2",
    "acmg:bp3",
    "acmg:bp4",
    "acmg:bp5",
    "acmg:bp6",
    "acmg:bp7",
    "acmg:bs1",
    "acmg:bs2",
    "acmg:bs3",
    "acmg:bs4",
    "acmg:pm1",
    "acmg:pm2",
    "acmg:pm3",
    "acmg:pm4",
    "acmg:pm5",
    "acmg:pm6",
    "acmg:pp1",
    "acmg:pp2",
    "acmg:pp3",
    "acmg:pp4",
    "acmg:pp5",
    "acmg:ps1",
    "acmg:ps2",
    "acmg:ps3",
    "acmg:ps4",
    "acmg:pvs1",
    "affected_status",
    "allele_origin",
    "allele_origin_confirmation",
    "amp:level_a",
    "amp:level_b",
    "amp:level_c",
    "amp:level_d",
    "assertion_method",
    "c_hgvs",
    "clinical_significance",
    "condition",
    "curation_context",
    "curation_date",
    "curation_verified_date",
    "gene_penetrance",
    "gene_symbol",
    "genome_build",
    "horak:om1",
    "horak:om2",
    "horak:om3",
    "horak:om4",
    "horak:op1",
    "horak:op2",
    "horak:op3",
    "horak:op4",
    "horak:os1",
    "horak:os2",
    "horak:os3",
    "horak:ovs1",
    "horak:sbp1",
    "horak:sbp2",
    "horak:sbs1",
    "horak:sbs2",
    "horak:sbvs1",
    "interpretation_summary",
    "mechanism_of_disease",
    "mode_of_inheritance",
    "molecular_consequence",
    "p_hgvs",
    "sherloc:same_amino_acid_change_benign",
    "somatic:clinical_significance",
    "somatic:testing_context",
    "somatic:tumor_cellularity",
    "variant_class",
    "variant_inheritance",
    "variant_penetrance",
    "vcgs:0203",
    "vcgs:0207",
    "vcgs:0210",
    "vcgs:0212",
    "vcgs:0309",
    "vcgs:alternative_change",
    "vcgs:comparable_benign",
    "vcgs:constraint",
    "vcgs:disease_mechanism",
    "vcgs:nmd_benign",
    "vcgs:noncoding",
    "vcgs:proband",
    "zygosity"
])