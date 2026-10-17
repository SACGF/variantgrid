# State-only: takes over ClinVarKeyExcludePattern and its existing table from snpdb (#2276)

import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0231_conditiontext_monarch_search_failed'),
        ('snpdb', '0289_move_clinvarkeyexcludepattern_to_classification'),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.CreateModel(
                    name='ClinVarKeyExcludePattern',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                        ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                        ('evidence_key', models.TextField()),
                        ('pattern', models.TextField()),
                        ('name', models.TextField(blank=True)),
                        ('case_insensitive', models.BooleanField(blank=True, default=True)),
                        ('mode', models.TextField(choices=[('X', 'Exclude If Match'), ('N', 'Exclude If No Match')], default='X')),
                        ('clinvar_key', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='snpdb.clinvarkey')),
                    ],
                    options={
                        'db_table': 'snpdb_clinvarkeyexcludepattern',
                    },
                ),
            ],
            database_operations=[],
        ),
    ]
