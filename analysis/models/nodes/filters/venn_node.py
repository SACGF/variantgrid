import logging
import sys
import traceback
from typing import Optional, Tuple, Dict

import celery
from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q
from django.db.models.signals import post_delete
from django.dispatch import receiver
from lazy import lazy

from analysis.models.enums import SetOperations
from analysis.models.nodes.analysis_node import AnalysisNode, NodeStatus, NodeVersion
from library.database_utils import queryset_to_sql
from snpdb.models import VariantCollection, ProcessingStatus
from snpdb.variant_collection import write_variant_set_operation


class VennNode(AnalysisNode):
    LEFT_PARENT = 'left'
    RIGHT_PARENT = 'right'

    set_operation = models.CharField(max_length=1, choices=SetOperations.choices, default=SetOperations.INTERSECTION)
    # We need to keep track of which is left/right but use the original DAG
    # as we need to use existing code to go through graph
    left_parent = models.ForeignKey(AnalysisNode, null=True, related_name='left_parent_node', on_delete=SET_NULL)
    right_parent = models.ForeignKey(AnalysisNode, null=True, related_name='right_parent_node', on_delete=SET_NULL)
    min_inputs = 2
    max_inputs = 2
    uses_parent_queryset = False

    def add_parent(self, parent, *args, **kwargs):
        side = kwargs.pop("side")

        if side == VennNode.LEFT_PARENT:
            self.left_parent = parent
        elif side == VennNode.RIGHT_PARENT:
            self.right_parent = parent
        else:
            msg = f"{parent} add_parent not given 'side' information"
            raise ValueError(msg)

        return super().add_parent(parent, *args, **kwargs)

    def remove_parent(self, parent):
        if self.left_parent == parent:
            self.left_parent = None
        elif self.right_parent == parent:
            self.right_parent = None
        else:
            msg = f"Parent {parent} not on left or right side!"
            logging.error(msg)
            # Don't raise exception as we need to be able to delete broken nodes.

        logging.debug("END: left_parent = %s, right_parent = %s", self.left_parent, self.right_parent)
        return super().remove_parent(parent)

    def get_side_for_parent(self, parent):
        if self.left_parent == parent:
            return VennNode.LEFT_PARENT
        if self.right_parent == parent:
            return VennNode.RIGHT_PARENT
        msg = f"Parent {parent} not on left or right side!"
        logging.error(msg)
        raise ValueError(msg)

    def adjust_cloned_parents(self, old_new_map):
        """ Replaces left/right parents with values from dict. Keeps old if not present. """

        changed = False
        if self.left_parent:
            new_left_parent = old_new_map.get(self.left_parent.pk)
            if new_left_parent:
                self.left_parent = new_left_parent
                changed = True
        if self.right_parent:
            new_right_parent = old_new_map.get(self.right_parent.pk)
            if new_right_parent:
                self.right_parent = new_right_parent
                changed = True

        if changed:
            self.parents_changed = True
        return changed

    def get_connection_data(self, parent):
        values = super().get_connection_data(parent)
        side = self.get_side_for_parent(parent)
        values['side'] = side
        return values

    def get_rendering_args(self):
        return {"venn_flag": self.get_venn_flag()}

    @lazy
    def ordered_parents(self) -> Tuple[AnalysisNode, AnalysisNode]:
        """ Return left_parent, right_parent """
        parents = self.get_parent_subclasses()
        a, b = parents
        if self.left_parent.pk == a.pk and self.right_parent.pk == b.pk:
            pass
        elif self.left_parent.pk == b.pk and self.right_parent.pk == a.pk:
            a, b = b, a
        else:
            msg = "Left/Right all mixed up!"
            raise ValueError(msg)

        return a, b

    def get_venn_flag(self):
        return [t[0] for t in SetOperations.choices].index(self.set_operation)

    def get_cache_task_args_objs_set(self, force_cache=False):
        """ Override from AnalysisNode - returns Celery tasks which are called
            in node_utils.get_analysis_update_task before children are loaded  """

        task_args_objs_set = set()
        if self.is_valid():
            try:
                a, b = self.ordered_parents
                for intersection_type in self.get_vennodecache_intersection_types():
                    vennode_cache, created = VennNodeCache.objects.get_or_create(parent_a_node_version=a.node_version,
                                                                                 parent_b_node_version=b.node_version,
                                                                                 intersection_type=intersection_type)

                    if not created:
                        # Check variant collection is valid, otherwise set to None so it'll be regenerated
                        if vennode_cache.variant_collection.status == ProcessingStatus.SUCCESS:
                            logging.debug("Venn Cache %s still valid...", vennode_cache)
                        else:
                            status = vennode_cache.variant_collection.get_status_display()
                            logging.debug("Venn Cache %s had status %s, regenerating", vennode_cache, status)
                            vennode_cache.variant_collection.delete()
                            vennode_cache.variant_collection = None

                    if vennode_cache.variant_collection:
                        task = None
                        cache_args = None
                    else:
                        # Create all VariantCollections in a loop to avoid race condition
                        name = f"Venn Count for {vennode_cache}"
                        variant_collection = VariantCollection.objects.create(name=name)
                        vennode_cache.variant_collection = variant_collection
                        vennode_cache.save()
                        task = venn_cache_count
                        cache_args = (vennode_cache.pk, )

                    task_args_objs_set.add((task, cache_args, vennode_cache))
            except Exception as e:
                errors = f"get_cache_task_args_objs_set exception: {e}"
                self.errors = errors  # Setting error here is in effect permanent (never cleared)
                self.status = NodeStatus.ERROR
                self.save()
                logging.error(errors)

        return task_args_objs_set

    def get_vennodecache_intersection_types(self):
        INTERSECTIONS = {
            SetOperations.UNION: [VennNodeCache.A_ONLY, VennNodeCache.INTERSECTION, VennNodeCache.B_ONLY],
            SetOperations.A_NOT_B: [VennNodeCache.A_ONLY],
            SetOperations.INTERSECTION: [VennNodeCache.INTERSECTION],
            SetOperations.A_ONLY: [VennNodeCache.A_ONLY, VennNodeCache.INTERSECTION],
            SetOperations.B_NOT_A: [VennNodeCache.B_ONLY],
            SetOperations.SYMMETRIC_DIFFERENCE: [VennNodeCache.A_ONLY, VennNodeCache.B_ONLY],
            SetOperations.B_ONLY: [VennNodeCache.INTERSECTION, VennNodeCache.B_ONLY],
        }
        return INTERSECTIONS[SetOperations(self.set_operation)]

    def _get_node_q(self) -> Optional[Q]:
        raise ValueError("VennNode always uses cache - this should never be called!")

    def _get_node_cache_arg_q_dict(self) -> Dict[Optional[str], Q]:
        if self.set_operation == SetOperations.NONE:
            return {None: self.q_none()}

        a, b = self.ordered_parents
        variant_collections = []
        for intersection_type in self.get_vennodecache_intersection_types():
            vennode_cache = VennNodeCache.objects.get(parent_a_node_version=a.node_version,
                                                      parent_b_node_version=b.node_version,
                                                      intersection_type=intersection_type)

            if vennode_cache.variant_collection.status != ProcessingStatus.SUCCESS:
                raise ValueError(f"{vennode_cache} had status: {vennode_cache.variant_collection.get_status_display()}")
            variant_collections.append(vennode_cache.variant_collection)

        q = Q(variantcollectionrecord__variant_collection__in=variant_collections)
        return {None: q}

    def _get_method_summary(self):
        return self.get_set_operation_display()

    def get_node_name(self):
        return ""

    @staticmethod
    def get_help_text() -> str:
        return "Filter based on set intersections between 2 parent nodes"

    @staticmethod
    def get_node_class_label():
        return "Venn"


class VennNodeCache(models.Model):
    A_ONLY = 'A'
    INTERSECTION = 'I'
    B_ONLY = 'B'
    INTERSECTION_CHOICE = [
        (A_ONLY, "A_ONLY"),
        (INTERSECTION, "INTERSECTION"),
        (B_ONLY, "B_ONLY")
    ]
    parent_a_node_version = models.ForeignKey(NodeVersion, related_name='+', on_delete=CASCADE)
    parent_b_node_version = models.ForeignKey(NodeVersion, related_name='+', on_delete=CASCADE)
    intersection_type = models.CharField(max_length=1, choices=INTERSECTION_CHOICE, null=True)
    variant_collection = models.ForeignKey(VariantCollection, null=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ('parent_a_node_version', 'parent_b_node_version', 'intersection_type')

    def __str__(self):
        return f"A: ({self.parent_a_node_version}), B: ({self.parent_b_node_version}): {self.get_intersection_type_display()}"


@receiver(post_delete, sender=VennNodeCache)
def post_delete_intersection_cache(sender, instance, **kwargs):  # pylint: disable=unused-argument
    try:
        if instance.variant_collection:
            instance.variant_collection.delete_related_objects()
            instance.variant_collection.delete()
    except VariantCollection.DoesNotExist:
        pass  # OK as deleted elsewhere (eg version was bumped and old ones cleaned up)


@celery.shared_task
def venn_cache_count(vennode_cache_id):
    print(f"venn_cache_count: {vennode_cache_id}")
    try:
        vennode_cache = VennNodeCache.objects.get(pk=vennode_cache_id)
    except VennNodeCache.DoesNotExist:
        return  # obsolete
    vc = vennode_cache.variant_collection
    vc.status = ProcessingStatus.PROCESSING
    vc.save()

    try:
        a = AnalysisNode.objects.get_subclass(pk=vennode_cache.parent_a_node_version.node_id,
                                              version=vennode_cache.parent_a_node_version.version)
        b = AnalysisNode.objects.get_subclass(pk=vennode_cache.parent_b_node_version.node_id,
                                              version=vennode_cache.parent_b_node_version.version)
        a_qs = a.get_queryset()
        b_qs = b.get_queryset()

        a_sql = queryset_to_sql(a_qs.values_list('id'))
        b_sql = queryset_to_sql(b_qs.values_list('id'))

        if vennode_cache.intersection_type == VennNodeCache.A_ONLY:
            write_variant_set_operation(vc, a_sql, b_sql, 'EXCEPT')
        elif vennode_cache.intersection_type == VennNodeCache.INTERSECTION:
            write_variant_set_operation(vc, a_sql, b_sql, 'INTERSECT')
        elif vennode_cache.intersection_type == VennNodeCache.B_ONLY:
            write_variant_set_operation(vc, b_sql, a_sql, 'EXCEPT')

        # TODO: Do we need to set something here so it's ready?
        logging.debug("Done writing variant collection")
        vc.status = ProcessingStatus.SUCCESS
    except Exception as e:
        logging.error("Got exception: %s", e)
        exc_traceback = sys.exc_info()[2]
        tb = ''.join(traceback.format_tb(exc_traceback))
        logging.error(tb)
        vc.status = ProcessingStatus.ERROR

    vc.save()
