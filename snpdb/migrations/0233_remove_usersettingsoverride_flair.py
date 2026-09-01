from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0232_one_off_update_user_awards"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="usersettingsoverride",
            name="flair",
        ),
    ]
