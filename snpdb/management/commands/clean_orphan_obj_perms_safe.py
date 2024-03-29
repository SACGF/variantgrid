from itertools import chain

from django.core.management.base import BaseCommand
from guardian.utils import get_user_obj_perms_model, get_group_obj_perms_model


class Command(BaseCommand):
    """
    clean_orphan_obj_perms command is a tiny wrapper around
    :func:`guardian.utils.clean_orphan_obj_perms`.

    Usage::

        $ python manage.py clean_orphan_obj_perms
        Removed 11 object permission entries with no targets

    """
    help = "Removes object permissions with not existing targets"

    def handle(self, **options):
        UserObjectPermission = get_user_obj_perms_model()
        GroupObjectPermission = get_group_obj_perms_model()

        deleted = 0
        # TODO: optimise
        for perm in chain(UserObjectPermission.objects.all().iterator(),
                          GroupObjectPermission.objects.all().iterator()):
            content_object = None
            try:
                content_object = perm.content_object
                if content_object is None:
                    print("Removing orphan %s (pk=%d)" % (perm, perm.pk))
            except:
                print(f"** Error retrieving content_object for {perm.pk}, Removing")

            if content_object is None:
                perm.delete()
                deleted += 1
        print("Total removed orphan object permissions instances: %d" %
                    deleted)
