from django.db import migrations


class Migration(migrations.Migration):
    """ Throwaway: changes the CI test-DB cache key to exercise restore-keys. Never merged. """

    dependencies = [
        ("snpdb", "0018_one_off_split_shared_user_create_labs"),
    ]

    operations = []
