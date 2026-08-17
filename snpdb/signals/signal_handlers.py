from django.conf import settings
from django.contrib.auth.models import Group
from django.db import transaction

from analysis.tasks.karyomapping_tasks import create_genome_karyomapping_for_trio
from library.log_utils import AdminNotificationBuilder
from snpdb.models import Lab, Organization, SettingsInitialGroupPermission, UserDataPrefix
from snpdb.tasks.vcf_bed_file_task import create_backend_vcf_bed_intersections


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

        org_labs = getattr(settings, "USER_CREATE_ORG_LABS", None)
        if org_labs:
            for org_group_name, lab_pattern in org_labs.items():
                organization, _ = Organization.objects.get_or_create(group_name=org_group_name)
                lab_name = lab_pattern % instance.__dict__
                lab_group_name = f"{organization.group_name}/{lab_name.lower()}"
                lab, _ = Lab.objects.get_or_create(name=lab_name, organization=organization, group_name=lab_group_name)
                organization.add_member(instance)
                add_user_to_group(lab_group_name)

        nb = AdminNotificationBuilder(message="User Created")
        nb.add_markdown("A new user has been created")
        nb.add_field("Username", instance.username)
        nb.add_field("Email", instance.email)
        nb.send()


def group_post_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        SettingsInitialGroupPermission.create_global_settings(instance)


def backend_vcf_import_success_handler(*args, **kwargs):
    backend_vcf = kwargs["backend_vcf"]

    create_backend_vcf_bed_intersections(backend_vcf)


def trio_post_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        celery_task = create_genome_karyomapping_for_trio.si(instance.pk)  # @UndefinedVariable
        transaction.on_commit(celery_task.apply_async)
