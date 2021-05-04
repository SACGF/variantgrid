from django.core.management.base import BaseCommand

from library.guardian_utils import assign_permission_to_user_and_groups, add_public_group_read_permission
from snpdb.models import CustomColumnsCollection


class Command(BaseCommand):
    """ This needs to be run on all VG systems to set initial permissions  """
    def handle(self, *args, **options):
        # User as NULL = public
        for ccc in CustomColumnsCollection.objects.filter(user__isnull=True):
            add_public_group_read_permission(ccc)

        for ccc in CustomColumnsCollection.objects.filter(user__isnull=False):
            assign_permission_to_user_and_groups(ccc.user, ccc)
