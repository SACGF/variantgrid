# Generated by Django 3.1.3 on 2021-04-09 07:54

import django.db.models.deletion
import django_extensions.db.fields
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0030_one_off_fix_cohort_sample_order'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('analysis', '0033_one_off_variant_tags'),
    ]

    operations = [
        migrations.AlterField(
            model_name='varianttag',
            name='genome_build',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.genomebuild'),
        ),
        migrations.CreateModel(
            name='VariantTagsImport',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('genome_build', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.genomebuild')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ImportedVariantTag',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('variant_string', models.TextField()),
                ('genome_build_string', models.TextField()),
                ('gene_symbol_string', models.TextField()),
                ('tag_string', models.TextField()),
                ('variant_id', models.IntegerField()),
                ('analysis_id', models.IntegerField(null=True)),
                ('node_id', models.IntegerField(null=True)),
                ('analysis_name', models.TextField()),
                ('user_name', models.TextField()),
                ('created', models.DateTimeField()),
                ('variant_tags_import', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='analysis.varianttagsimport')),
            ],
        ),
        migrations.AddField(
            model_name='varianttag',
            name='imported_variant_tag',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='analysis.importedvarianttag'),
        ),
    ]