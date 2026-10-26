from django.db import migrations

from library.django_utils.seed_data import insert_seed_data
from manual.operations.manual_operations import ManualOperation

NODE_GRAPH_TYPES = [
    {
        'name': 'Homozygosity percent',
        'required_column': 'observedvariant__zygosity',
        'graph_class': 'analysis.graphs.node_homozygosity_percent_graph.NodeHomozygosityPercentGraph',
        'description': 'Homozygosity percent',
    },
    {
        'name': 'Chromosome SNP Density',
        'graph_class': 'analysis.graphs.chromosome_density_graph.ChromosomeDensityGraph',
        'description': 'Plot density of SNPs in 1Mb bins',
    },
    {
        'name': 'Allele Frequency Histogram',
        'required_column': 'observedvariant__allele_frequency',
        'graph_class': 'analysis.graphs.node_allele_frequency_histogram_graph.NodeAlleleFrequencyHistogramGraph',
        'description': 'AF Histogram (single input sample)',
    },
]


def _test_no_templates(apps):
    AnalysisTemplate = apps.get_model("analysis", "AnalysisTemplate")
    return not AnalysisTemplate.objects.exists()


def _insert_seed_data(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [
        ("analysis.NodeGraphType", NODE_GRAPH_TYPES),
    ])


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('analysis', '0002_initial'),
        ('manual', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(_insert_seed_data),
        ManualOperation.operation_manage(["analysis_create_default_templates"], test=_test_no_templates,
                                         requires=["annotation-version"]),
    ]
