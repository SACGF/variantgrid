"""
Cytoband rows per contig, read from each build's settings.ANNOTATION["cytoband"] file, so a band such as 11p13
resolves to coordinates with a query (#2685). The files are 0-based half-open; rows are stored 1-based inclusive.
"""
import os

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models

from library.graphs.chromosomes_graph import read_cytoband


def _load_cytobands(apps, _schema_editor):
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    GenomeBuildContig = apps.get_model("snpdb", "GenomeBuildContig")
    Cytoband = apps.get_model("snpdb", "Cytoband")

    for genome_build in GenomeBuild.objects.all():
        cytoband_filename = settings.ANNOTATION.get(genome_build.name, {}).get("cytoband")
        if not (cytoband_filename and os.path.exists(cytoband_filename)):
            continue
        contig_by_ucsc_name = {gbc.contig.ucsc_name: gbc.contig
                               for gbc in GenomeBuildContig.objects.filter(genome_build=genome_build)
                               .select_related("contig")}
        cytobands = []
        for chrom, bands in read_cytoband(cytoband_filename).items():
            if contig := contig_by_ucsc_name.get(chrom):
                for start, end, name, stain in bands:
                    cytobands.append(Cytoband(contig=contig, name=name, start=start + 1, end=end, stain=stain))
        Cytoband.objects.bulk_create(cytobands, ignore_conflicts=True)


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0010_sequence_length_index"),
    ]

    operations = [
        migrations.CreateModel(
            name="Cytoband",
            fields=[
                (
                    "id",
                    models.AutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("name", models.TextField()),
                ("start", models.IntegerField()),
                ("end", models.IntegerField()),
                ("stain", models.TextField()),
                (
                    "contig",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE, to="snpdb.contig"
                    ),
                ),
            ],
            options={
                "unique_together": {("contig", "name")},
            },
        ),
        migrations.RunPython(_load_cytobands, reverse_code=migrations.RunPython.noop),
    ]
