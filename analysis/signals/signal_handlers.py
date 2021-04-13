from django.db.models import Q

from analysis.models import TagNode, Analysis, Tag
from analysis.tasks.variant_tag_tasks import analysis_tag_created_task, analysis_tag_deleted_task


def _analysis_tag_nodes_set_dirty(analysis: Analysis, tag: Tag):
    """ Needs to be sync so version is bumped by the time client calls node_versions to see whats dirty """
    tag_filter = Q(tagnodetag__tag__isnull=True) | Q(tagnodetag__tag=tag)
    for node in TagNode.objects.filter(analysis=analysis).filter(tag_filter).distinct():
        node.queryset_dirty = True
        node.save()


def variant_tag_create(sender, instance, created=False, **kwargs):
    if created:
        if instance.analysis:
            _analysis_tag_nodes_set_dirty(instance.analysis, instance.tag)
            # want to be as quick as we can so do analysis reload + liftover async
            analysis_tag_created_task.si(instance.pk).apply_async()


def variant_tag_delete(sender, instance, **kwargs):
    if instance.analysis:
        _analysis_tag_nodes_set_dirty(instance.analysis, instance.tag)
        # want to be as quick as we can so do analysis reload + liftover async
        analysis_tag_deleted_task.si(instance.analysis_id, instance.tag_id).apply_async()
