from django.db import migrations

LEGACY_TABLE = "django_messages_message"
NEW_TABLE = "user_messages_message"

COLUMNS = (
    "id, subject, body, sent_at, read_at, replied_at, "
    "sender_deleted_at, recipient_deleted_at, parent_msg_id, recipient_id, sender_id"
)


def migrate_legacy_messages(apps, schema_editor):
    """
    On existing deployments, copy rows from the old ``django-messages`` table into the new
    ``user_messages`` table and drop the legacy table. On fresh installs the legacy table does not
    exist, so this is a no-op.
    """
    connection = schema_editor.connection
    with connection.cursor() as cursor:
        cursor.execute("SELECT to_regclass(%s)", [LEGACY_TABLE])
        if cursor.fetchone()[0] is None:
            return

        # Ordered by id so self-referential parent_msg rows (older = smaller id) insert first.
        cursor.execute(
            f"INSERT INTO {NEW_TABLE} ({COLUMNS}) "
            f"SELECT {COLUMNS} FROM {LEGACY_TABLE} ORDER BY id"
        )
        cursor.execute(
            f"SELECT setval(pg_get_serial_sequence('{NEW_TABLE}', 'id'), "
            f"COALESCE((SELECT MAX(id) FROM {NEW_TABLE}), 1))"
        )
        cursor.execute(f"DROP TABLE {LEGACY_TABLE}")


class Migration(migrations.Migration):

    dependencies = [
        ('user_messages', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(migrate_legacy_messages, migrations.RunPython.noop),
    ]
