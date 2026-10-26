from django.db import migrations

from library.django_utils.seed_data import insert_seed_data

GENE_LIST_CATEGORIES = [
    {'name': 'CustomText', 'icon_css_class': '', 'description': ''},
    {'name': 'PathologyTest', 'icon_css_class': '', 'description': ''},
    {'name': 'NodeCustomText', 'icon_css_class': '', 'hidden': True, 'description': ''},
    {'name': 'SampleGeneList', 'icon_css_class': '', 'description': ''},
    {'name': 'Uploaded', 'icon_css_class': '', 'description': ''},
    {'name': 'GeneInfo', 'icon_css_class': '', 'description': 'Special category to hold gene info (ie tag genes)'},
]

GENE_SYMBOLS = [
    {'symbol': 'AR'},
    {'symbol': 'EGFR'},
    {'symbol': 'MET'},
]

PANEL_APP_SERVERS = [
    {
        'name': 'Genomics England PanelApp',
        'url': 'https://panelapp.genomicsengland.co.uk',
        'icon_css_class': 'panel-app-england-icon',
    },
    {'name': 'PanelApp Australia', 'url': 'https://panelapp-aus.org', 'icon_css_class': 'panel-app-australia-icon'},
]

SPLICE_EVENTS = [
    {
        'gene_symbol_id': 'AR',
        'label': 'v_7',
        'display': 'AR-V7 splice variant',
        'genome_build_id': 'GRCh37',
        'contig_id': 23,
        'donor': 66905968,
        'acceptor': 66914514,
        'civic_variant_id': 362,
    },
    {
        'gene_symbol_id': 'AR',
        'label': 'v_7',
        'display': 'AR-V7 splice variant',
        'genome_build_id': 'GRCh38',
        'contig_id': 322,
        'donor': 67686126,
        'acceptor': 67694672,
        'civic_variant_id': 362,
    },
    {
        'gene_symbol_id': 'EGFR',
        'label': 'v_iii',
        'display': 'EGFRvIII splice variant',
        'genome_build_id': 'GRCh37',
        'contig_id': 7,
        'donor': 55087058,
        'acceptor': 55223522,
        'civic_variant_id': 312,
    },
    {
        'gene_symbol_id': 'EGFR',
        'label': 'v_iii',
        'display': 'EGFRvIII splice variant',
        'genome_build_id': 'GRCh38',
        'contig_id': 306,
        'donor': 55019365,
        'acceptor': 55155829,
        'civic_variant_id': 312,
    },
    {
        'gene_symbol_id': 'MET',
        'label': 'exon_14_skipping',
        'display': 'MET exon 14 skipping',
        'genome_build_id': 'GRCh37',
        'contig_id': 7,
        'donor': 116411708,
        'acceptor': 116414934,
        'civic_variant_id': 324,
    },
    {
        'gene_symbol_id': 'MET',
        'label': 'exon_14_skipping',
        'display': 'MET exon 14 skipping',
        'genome_build_id': 'GRCh38',
        'contig_id': 306,
        'donor': 116771654,
        'acceptor': 116774880,
        'civic_variant_id': 324,
    },
    {
        'gene_symbol_id': 'EGFR',
        'label': 'v_ii',
        'display': 'EGFRvII splice variant',
        'genome_build_id': 'GRCh37',
        'contig_id': 7,
        'donor': 55229324,
        'acceptor': 55238867,
    },
    {
        'gene_symbol_id': 'EGFR',
        'label': 'v_ii',
        'display': 'EGFRvII splice variant',
        'genome_build_id': 'GRCh38',
        'contig_id': 306,
        'donor': 55161631,
        'acceptor': 55171174,
    },
]


def _insert_seed_data(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [
        ("genes.GeneListCategory", GENE_LIST_CATEGORIES),
        ("genes.GeneSymbol", GENE_SYMBOLS),
        ("genes.PanelAppServer", PANEL_APP_SERVERS),
        ("genes.SpliceEvent", SPLICE_EVENTS),
    ])


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('genes', '0002_initial'),
        ('snpdb', '0001_initial'),
        ('snpdb', '0002_seed_data'),
    ]

    operations = [
        migrations.RunPython(_insert_seed_data),
    ]
