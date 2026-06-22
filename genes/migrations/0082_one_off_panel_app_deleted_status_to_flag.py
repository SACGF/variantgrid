# One-off historical fix for issue #405.
#
# An earlier patched build wrote the sentinel string "deleted" into
# PanelAppPanel.status to mark panels that PanelApp had removed (404). The
# canonical representation is now the dedicated boolean PanelAppPanel.deleted
# (added in 0081), and status mirrors PanelApp's own values
# ("public" / "internal" / "promoted" / "retired") again.
#
# This migration carries those rows forward by flipping deleted=True for any
# row currently holding status="deleted". status is left untouched here; the
# next refresh from PanelApp will overwrite it with the real upstream value
# (or, for genuinely missing panels, the row stays with a stale status — which
# is fine, since deleted=True is what callers now check).
from django.db import migrations


def patch_deleted_status_to_flag(apps, schema_editor):
    PanelAppPanel = apps.get_model("genes", "PanelAppPanel")
    PanelAppPanel.objects.filter(status="deleted").update(deleted=True)


class Migration(migrations.Migration):

    dependencies = [
        ("genes", "0081_panelapppanel_deleted"),
    ]

    operations = [
        migrations.RunPython(patch_deleted_status_to_flag, migrations.RunPython.noop),
    ]
