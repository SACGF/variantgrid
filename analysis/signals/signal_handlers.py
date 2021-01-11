import logging
from time import time

from django.db.models import Q

from analysis.models import Analysis, Tag, TagNode
from analysis.models.nodes.node_utils import update_nodes


def _update_analysis_on_variant_tag_change(analysis: Analysis, tag: Tag):
    """ TagNodes relying on VariantTag need to be reloaded """

    start = time()
    tag_filter = Q(tag__isnull=True) | Q(tag=tag)
    # TagNode doesn't have any output, so no need to check children
    need_to_reload = False
    for node in TagNode.objects.filter(analysis=analysis).filter(tag_filter):
        node.queryset_dirty = True
        node.save()
        need_to_reload = True

    if need_to_reload:
        logging.info("Reloading %s (%d) as tag %s changed", analysis, analysis.pk, tag)
        update_nodes(analysis.pk)

    end = time()
    time_taken = end - start
    logging.info("Time taken %.3f", time_taken)


def variant_tag_create(sender, instance, created=False, **kwargs):
    if created:
        _update_analysis_on_variant_tag_change(instance.analysis, instance.tag)


def variant_tag_delete(sender, instance, **kwargs):
    _update_analysis_on_variant_tag_change(instance.analysis, instance.tag)
