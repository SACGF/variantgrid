from django.conf import settings
from django.contrib.auth.models import Group
from django.db import transaction

from analysis.tasks.karyomapping_tasks import create_genome_karyomapping_for_trio
from eventlog.notifications import AdminNotificationBuilder
from snpdb.models.models_lab_request import create_user_org_labs
from snpdb.models.models_user_settings import SettingsInitialGroupPermission, UserDataPrefix


def user_post_save_handler(sender, instance, **kwargs):
    """ Add new user to the public group """

    def add_user_to_group(group_name: str):
        group, _ = Group.objects.get_or_create(name=group_name)
        group.user_set.add(instance)

    if kwargs.get("created"):
        add_user_to_group(settings.PUBLIC_GROUP_NAME)
        add_user_to_group(settings.LOGGED_IN_USERS_GROUP_NAME)

        udp_kwargs = getattr(settings, "INITIAL_USER_DATA_PREFIX_KWARGS", None)
        if udp_kwargs:
            UserDataPrefix.objects.get_or_create(user=instance, **udp_kwargs)

        create_user_org_labs(instance)

        nb = AdminNotificationBuilder(message="User Created")
        nb.add_markdown("A new user has been created")
        nb.add_field("Username", instance.username)
        nb.add_field("Email", instance.email)
        nb.send()


def group_post_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        SettingsInitialGroupPermission.create_global_settings(instance)


def trio_post_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        celery_task = create_genome_karyomapping_for_trio.si(instance.pk)  # @UndefinedVariable
        transaction.on_commit(celery_task.apply_async)
