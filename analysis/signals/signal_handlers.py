from analysis.tasks.variant_tag_tasks import analysis_tag_created_task, analysis_tag_deleted_task


def variant_tag_create(sender, instance, created=False, **kwargs):
    if created:
        # Need to do liftover so do it all in an async task
        analysis_tag_created_task.si(instance.pk).apply_async()


def variant_tag_delete(sender, instance, **kwargs):
    analysis_tag_created_task.si(instance.pk)
    # This is mostly quick, launching other tasks so call method rather than
    analysis_tag_deleted_task.si(instance.analysis_id, instance.tag_id).apply_async()
