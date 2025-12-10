from django.db import transaction

from analysis.analysis_templates import auto_launch_analysis_templates_for_sample
from analysis.tasks.variant_tag_tasks import variant_tag_created_task, variant_tag_deleted_in_analysis_task, \
    analysis_tag_nodes_set_dirty
from library.guardian_utils import admin_bot


def variant_tag_create(sender, instance, created=False, **kwargs):
    if created:
        if instance.analysis:
            analysis_tag_nodes_set_dirty(instance.analysis, instance.tag, visible=True)
        # want to be as quick as we can so do analysis reload + liftover async
        # Need to launch this at end of transaction so we know VariantTag is in DB for celery job
        celery_task = variant_tag_created_task.si(instance.pk)
        transaction.on_commit(lambda: celery_task.apply_async())


def variant_tag_delete(sender, instance, **kwargs):
    if instance.analysis:
        analysis_tag_nodes_set_dirty(instance.analysis, instance.tag, visible=True)
        # want to be as quick as we can so do analysis reload + liftover async
        celery_task = variant_tag_deleted_in_analysis_task.si(instance.analysis_id, instance.tag_id)
        transaction.on_commit(lambda: celery_task.apply_async())


def vcf_import_success(sender, instance, **kwargs):
    user = admin_bot()
    for sample in instance.sample_set.all():
        auto_launch_analysis_templates_for_sample(user, sample)
