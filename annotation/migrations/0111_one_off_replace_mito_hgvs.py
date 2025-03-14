# Generated by Django 4.2.11 on 2024-07-09 11:13

from django.db import migrations
from django.db.models import Value
from django.db.models.functions import Replace


def _one_off_replace_mito_hgvs(apps, _schema_editor):
    VariantAnnotationVersion = apps.get_model('annotation', 'VariantAnnotationVersion')
    VariantAnnotation = apps.get_model('annotation', 'VariantAnnotation')
    GenomeBuild = apps.get_model('snpdb', 'GenomeBuild')
    Contig = apps.get_model('snpdb', 'Contig')

    for genome_build_name in ["GRCh37", "GRCh38"]:
        genome_build = GenomeBuild.objects.get(name=genome_build_name)
        if mito_contig := Contig.objects.filter(molecule_type='M', genomebuildcontig__genome_build=genome_build).first():
            for vav in VariantAnnotationVersion.objects.filter(genome_build=genome_build):
                qs = VariantAnnotation.objects.filter(version=vav, variant__locus__contig=mito_contig)
                mito_accession = mito_contig.refseq_accession
                print(f"Updating VariantAnnotationVersion={vav.pk} MT accession {mito_accession} from '.g' -> '.m'")
                qs = qs.filter(hgvs_g__startswith=f'{mito_accession}:g.')
                num = qs.update(hgvs_g=Replace('hgvs_g', Value(f'{mito_accession}:g.'), Value(f'{mito_accession}:m.')))
                print(f"Updated {num} records.")


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0110_merge_20240709_1411"),
    ]

    operations = [
        migrations.RunPython(_one_off_replace_mito_hgvs)
    ]
