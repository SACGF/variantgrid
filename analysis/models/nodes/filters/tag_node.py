import operator
import time
from functools import reduce
from typing import List

from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q
from lazy import lazy

from analysis.models.enums import TagNodeMode
from analysis.models.models_variant_tag import VariantTag
from analysis.models.nodes.analysis_node import AnalysisNode
from snpdb.models import Tag


class TagNode(AnalysisNode):
    ANALYSIS_TAGS_NAME = "Tagged Variants"
    parent_input = models.BooleanField(default=True)
    exclude = models.BooleanField(default=False)
    mode = models.CharField(max_length=1, choices=TagNodeMode.choices, default=TagNodeMode.THIS_ANALYSIS)

    def modifies_parents(self):
        return True

    @lazy
    def tag_ids(self) -> List[str]:
        # This is called when Node is being initialised to set the name
        if self.pk is None:
            return []

        qs = Tag.objects.filter(tagnodetag__tag_node=self).order_by("id").values_list("id", flat=True)
        return list(qs)

    def _get_node_q(self) -> Q:
        # Pull in tags from this analysis - use variant query
        # VariantTags are same build as analysis, so use this not Allele as it avoids a race condition where
        # tagging a variant w/o an Allele takes a few seconds to create one via liftover pipelines
        variants_with_tags = VariantTag.objects.filter(analysis=self.analysis)
        if self.tag_ids:
            variants_with_tags = variants_with_tags.filter(tag__in=self.tag_ids)
        q_list = [Q(pk__in=variants_with_tags.values_list("variant_id"))]

        if self.mode == TagNodeMode.ALL_TAGS:
            tags_qs = VariantTag.filter_for_user(self.analysis.user)
            # We already have tags from this analysis, no need to retrieve again
            tags_qs = tags_qs.exclude(analysis=self.analysis)
            # Builds from different analyses (maybe diff builds) - so do query using Allele
            q_list.append(VariantTag.variants_for_build_q(self.analysis.genome_build, tags_qs, self.tag_ids))

        q = reduce(operator.or_, q_list)
        if self.exclude:
            q = ~q
        return q

    def get_node_name(self):
        if self.visible:
            description_list = []
            if self.exclude:
                description_list.append("Exclude")

            if self.mode == TagNodeMode.ALL_TAGS:
                description_list.append("Global")
            else:
                description_list.append("Analysis")

            if self.tag_ids:
                description_list.append(', '.join(self.tag_ids))
            else:
                description_list.append("Tags")

            description = " ".join(description_list)
        else:
            description = self.ANALYSIS_TAGS_NAME  # Has to be set to this

        return description

    @staticmethod
    def get_help_text() -> str:
        return "Filter to variants that have been tagged"

    def save_clone(self):
        tag_ids = self.tag_ids  # Save before clone
        copy = super().save_clone()
        for tag_id in tag_ids:
            copy.tagnodetag_set.create(tag_id=tag_id)
        return copy

    @staticmethod
    def get_node_class_label():
        return "Tags"

    def _get_method_summary(self):
        return f"Tagged {', '.join(self.tag_ids)} ({self.get_mode_display()})"

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.tag_ids:
            css_classes.extend([f"tagged-{tag_id}" for tag_id in self.tag_ids])
        return css_classes

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.parent_input:
            return 1
        return 0

    @staticmethod
    def get_analysis_tags_node(analysis):
        from analysis.tasks.node_update_tasks import update_node_task

        node, created = TagNode.objects.get_or_create(analysis=analysis,
                                                      name=TagNode.ANALYSIS_TAGS_NAME,
                                                      parent_input=False,
                                                      mode=TagNodeMode.THIS_ANALYSIS,
                                                      visible=False)
        if created:
            # Should be fast, so do sync (not as celery job)
            update_node_task(node.pk, node.version)
        return node


class TagNodeTag(models.Model):
    """ Stores multi-select values """
    tag_node = models.ForeignKey(TagNode, on_delete=CASCADE)
    tag = models.ForeignKey(Tag, null=True, blank=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ("tag_node", "tag")
