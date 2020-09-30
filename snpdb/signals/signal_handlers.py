from django.conf import settings
from django.contrib.auth.models import Group

from analysis.tasks.karyomapping_tasks import create_genome_karyomapping_for_trio
from snpdb.models import UserDataPrefix, SettingsInitialGroupPermission
from snpdb.tasks.vcf_bed_file_task import create_backend_vcf_bed_intersections


def user_post_save_handler(sender, instance, **kwargs):
    """ Add new user to the public group """

    created = kwargs.get("created")
    if created:
        public_group, _ = Group.objects.get_or_create(name=settings.PUBLIC_GROUP_NAME)
        public_group.user_set.add(instance)

        public_group, _ = Group.objects.get_or_create(name=settings.LOGGED_IN_USERS_GROUP_NAME)
        public_group.user_set.add(instance)

        udp_kwargs = getattr(settings, "INITIAL_USER_DATA_PREFIX_KWARGS", None)
        if udp_kwargs:
            UserDataPrefix.objects.get_or_create(user=instance, **udp_kwargs)


def group_post_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        SettingsInitialGroupPermission.create_global_settings(instance)


def backend_vcf_import_success_handler(*args, **kwargs):
    backend_vcf = kwargs["backend_vcf"]

    create_backend_vcf_bed_intersections(backend_vcf)


def trio_save_handler(sender, instance, **kwargs):
    created = kwargs.get("created")
    if created:
        task = create_genome_karyomapping_for_trio.si(instance.pk)  # @UndefinedVariable
        task.apply_async()
