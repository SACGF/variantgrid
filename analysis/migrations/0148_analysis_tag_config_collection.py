from django.db import migrations, models
from django.db.models.deletion import SET_NULL


def resolve_tag_config_collection_id(apps, user_id):
    """ UserSettings.get_settings_overrides' order, against the historical models:
        Global -> the user's default lab's organization -> that lab -> the user, last non-null wins """
    global_settings = apps.get_model("snpdb", "GlobalSettings")
    org_override_model = apps.get_model("snpdb", "OrganizationUserSettingsOverride")
    lab_override_model = apps.get_model("snpdb", "LabUserSettingsOverride")
    user_override_model = apps.get_model("snpdb", "UserSettingsOverride")

    collection_id = None
    if gs := global_settings.objects.first():
        collection_id = gs.tag_config_id or collection_id

    user_override = user_override_model.objects.filter(user_id=user_id).first()
    lab = user_override.default_lab if user_override else None
    if lab:
        if org_override := org_override_model.objects.filter(organization=lab.organization).first():
            collection_id = org_override.tag_config_id or collection_id
        if lab_override := lab_override_model.objects.filter(lab=lab).first():
            collection_id = lab_override.tag_config_id or collection_id
    if user_override:
        collection_id = user_override.tag_config_id or collection_id
    return collection_id


def set_analysis_tag_config_collection(apps, _schema_editor=None):
    """ Existing analyses keep showing what their owner sees today """
    analysis_model = apps.get_model("analysis", "Analysis")
    collection_id_by_user_id = {}
    for analysis in analysis_model.objects.all().only("pk", "user_id").iterator():
        user_id = analysis.user_id
        if user_id not in collection_id_by_user_id:
            collection_id_by_user_id[user_id] = resolve_tag_config_collection_id(apps, user_id)
        if collection_id := collection_id_by_user_id[user_id]:
            analysis_model.objects.filter(pk=analysis.pk).update(tag_config_collection_id=collection_id)


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0147_allvariantsnode_gene_level_types'),
        ('snpdb', '0273_rename_tagcolor_tagconfig'),
    ]

    operations = [
        migrations.AddField(
            model_name='analysis',
            name='tag_config_collection',
            field=models.ForeignKey(blank=True, help_text='Tag colours, sort order, 1-click tags and staleness everyone opening this analysis sees', null=True, on_delete=SET_NULL, to='snpdb.tagconfigcollection'),
        ),
        migrations.RunPython(set_analysis_tag_config_collection, migrations.RunPython.noop),
        migrations.RemoveField(
            model_name='analysis',
            name='variant_tag_stale_days',
        ),
    ]
