from django.core.management.base import BaseCommand

# Can be removed once environments migrated through snpdb/0089
from library.guardian_utils import (
    add_public_group_read_permission,
    assign_permission_to_user_and_groups,
)
from snpdb.models import TagColorsCollection


class Command(BaseCommand):
    def handle(self, *args, **options):
        for collection in TagColorsCollection.objects.all():
            if collection.user:
                assign_permission_to_user_and_groups(collection.user, collection)
            else:
                add_public_group_read_permission(collection)
