# Generated by Django 4.1.4 on 2023-09-29 04:31

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0103_alter_lab_options_alter_organization_options"),
    ]

    operations = [
        migrations.AlterField(
            model_name="settingsoverride",
            name="columns",
            field=models.ForeignKey(
                blank=True,
                help_text="Initial custom columns when creating analysis",
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to="snpdb.customcolumnscollection",
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="default_genome_build",
            field=models.ForeignKey(
                blank=True,
                help_text="Used for search (jump if 1 result for this build) and populating defaults everywhere",
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to="snpdb.genomebuild",
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="default_sort_by_column",
            field=models.ForeignKey(
                blank=True,
                help_text="Initial sort by column for analysis variant grid",
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to="snpdb.customcolumn",
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="email_discordance_updates",
            field=models.BooleanField(
                blank=True, help_text="Opt into discordance email updates", null=True
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="email_weekly_updates",
            field=models.BooleanField(
                blank=True, help_text="Opt in to email list", null=True
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="igv_port",
            field=models.IntegerField(
                blank=True,
                help_text="Port to connect to IGV on your machine",
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="import_messages",
            field=models.BooleanField(
                blank=True,
                help_text="Get internal notification (message icon top right) when imports are done (eg VCF finished processing and annotating)",
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="node_debug_tab",
            field=models.BooleanField(
                blank=True,
                help_text="If true, an extra tab appears in analysis node editor, with details about node settings + SQL query.",
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="tag_colors",
            field=models.ForeignKey(
                blank=True,
                help_text="Custom tag colors used in analyis and variant grids",
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to="snpdb.tagcolorscollection",
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="timezone",
            field=models.TextField(
                blank=True,
                help_text="Time/date used in classification download",
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="tool_tips",
            field=models.BooleanField(
                blank=True, help_text="Show/hide help popups on mouse hover", null=True
            ),
        ),
        migrations.AlterField(
            model_name="settingsoverride",
            name="variant_link_in_analysis_opens_new_tab",
            field=models.BooleanField(
                help_text="Whether left click by default opens up variant details in new window. Default is open where node editor is. You can always open in new window via right click then new window",
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name="usersettingsoverride",
            name="default_lab",
            field=models.ForeignKey(
                blank=True,
                help_text="Lab used for creating classifications (you can belong to more than 1 lab)",
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to="snpdb.lab",
            ),
        ),
    ]
