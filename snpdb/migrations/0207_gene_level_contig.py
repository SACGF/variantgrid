from django.db import migrations

from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_REFSEQ_ACCESSION
from snpdb.genome.reference_contigs import create_gene_level_contig


def _add_gene_level_contig(apps, _schema_editor):
    def get_model(name):
        return apps.get_model("snpdb", name)

    GenomeBuild = get_model("GenomeBuild")
    GenomeBuildContig = get_model("GenomeBuildContig")

    for genome_build in GenomeBuild.objects.all():
        order = GenomeBuildContig.objects.filter(genome_build=genome_build).count()
        create_gene_level_contig(get_model, genome_build, order)


def _remove_gene_level_contig(apps, _schema_editor):
    Contig = apps.get_model("snpdb", "Contig")
    Contig.objects.filter(refseq_accession=GENE_LEVEL_CONTIG_REFSEQ_ACCESSION).delete()


class Migration(migrations.Migration):
    dependencies = [
        ('snpdb', '0206_sequence_role_gene_level'),
    ]

    operations = [
        migrations.RunPython(_add_gene_level_contig, _remove_gene_level_contig),
    ]
