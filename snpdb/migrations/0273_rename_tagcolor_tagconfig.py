from django.db import migrations, models
from django.db.models.deletion import SET_NULL

# The override levels in the order UserSettings merges them - later (more specific) wins
_OVERRIDE_MODELS = ["GlobalSettings", "OrganizationUserSettingsOverride",
                    "LabUserSettingsOverride", "UserSettingsOverride"]


def copy_stale_days_to_collections(apps, _schema_editor=None):
    """ Staleness used to be a per-user setting layered Global -> Org -> Lab -> User; it now lives on
        the collection those levels point at. Walk the levels in that order so the most specific
        override's value is what lands on a collection, and let the first user-level override that
        writes a collection keep it (two users sharing one collection can't both win) """
    written_by_user_override = set()
    for model_name in _OVERRIDE_MODELS:
        model = apps.get_model("snpdb", model_name)
        user_level = model_name == "UserSettingsOverride"
        for override in model.objects.filter(tag_config__isnull=False,
                                             variant_tag_stale_days__isnull=False):
            collection = override.tag_config
            if user_level:
                if collection.pk in written_by_user_override:
                    continue
                written_by_user_override.add(collection.pk)
            collection.variant_tag_stale_days = override.variant_tag_stale_days
            collection.save()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0272_tagcolor_quick_tag'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='TagColorsCollection',
            new_name='TagConfigCollection',
        ),
        migrations.RenameModel(
            old_name='TagColor',
            new_name='TagConfig',
        ),
        migrations.RenameField(
            model_name='settingsoverride',
            old_name='tag_colors',
            new_name='tag_config',
        ),
        migrations.AlterField(
            model_name='settingsoverride',
            name='tag_config',
            field=models.ForeignKey(blank=True, help_text="Tag colours, sort order, 1-click tags and staleness (modify/create these in 'Tag settings'). Initial tag config when creating an analysis", null=True, on_delete=SET_NULL, to='snpdb.tagconfigcollection'),
        ),
        migrations.AddField(
            model_name='tagconfigcollection',
            name='variant_tag_stale_days',
            field=models.IntegerField(blank=True, choices=[(180, '6 months'), (365, '1 year'), (545, '18 months'), (730, '2 years'), (1825, '5 years')], help_text='Tag events older than this are considered stale: grids show fresh vs total counts and mark tags whose most recent event is older. Blank disables staleness.', null=True),
        ),
        migrations.RunPython(copy_stale_days_to_collections, migrations.RunPython.noop),
        migrations.RemoveField(
            model_name='settingsoverride',
            name='variant_tag_stale_days',
        ),
    ]
