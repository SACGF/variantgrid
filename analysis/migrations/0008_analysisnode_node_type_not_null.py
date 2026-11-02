from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0007_analysisnode_node_type'),
    ]

    operations = [
        migrations.AlterField(
            model_name='analysisnode',
            name='node_type',
            field=models.TextField(editable=False),
        ),
    ]
