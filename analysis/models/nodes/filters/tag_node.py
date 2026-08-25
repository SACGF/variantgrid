import operator
from datetime import datetime, timedelta
from functools import cached_property, reduce
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.query_utils import Q
from django.utils import timezone
from django.utils.timezone import localtime

from analysis.models.enums import TagNodeInput, TagNodeMode
from analysis.models.models_variant_tag import VariantTag
from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin, NodeVersion
from analysis.models.nodes.node_display import NodeIcon
from snpdb.models import Tag


class TagNode(AnalysisNode):
    ANALYSIS_TAGS_NAME = "Tagged Variants"
    node_input = models.CharField(max_length=1, choices=TagNodeInput.choices, default=TagNodeInput.PARENT_TAGGED)
    mode = models.CharField(max_length=1, choices=TagNodeMode.choices, default=TagNodeMode.THIS_ANALYSIS)
    tagged_within_days = models.IntegerField(null=True, blank=True)

    def modifies_parents(self):
        return True

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if self.mode == TagNodeMode.ALL_TAGS:
            msg = "Global tags are a snapshot"
            if node_version := NodeVersion.objects.filter(node=self, version=self.version).first():
                msg += f" taken {localtime(node_version.created):%d %b %Y %H:%M}"
            msg += " - tags added or removed in other analyses since then won't be reflected here." \
                   " Press save to refresh."
            warnings.append(msg)
        return warnings

    @cached_property
    def tag_ids(self) -> list[str]:
        # This is called when Node is being initialised to set the name
        if self.pk is None:
            return []

        qs = Tag.objects.filter(tagnodetag__tag_node=self).order_by("id").values_list("id", flat=True)
        return list(qs)

    @property
    def tagged_within_cutoff(self) -> Optional[datetime]:
        """ Anchored to the node version's save date (not query time) so re-running an old analysis
            reproduces what it showed when saved - matching the global tags snapshot semantics """
        if self.tagged_within_days is None:
            return None
        node_version = NodeVersion.objects.filter(node=self, version=self.version).first()
        anchor = node_version.created if node_version else timezone.now()
        return anchor - timedelta(days=self.tagged_within_days)

    def _get_node_q(self) -> Q:
        cutoff = self.tagged_within_cutoff
        # Pull in tags from this analysis - use variant query
        # VariantTags are same build as analysis, so use this not Allele as it avoids a race condition where
        # tagging a variant w/o an Allele takes a few seconds to create one via liftover pipelines
        variants_with_tags = VariantTag.objects.filter(analysis=self.analysis)
        if self.tag_ids:
            variants_with_tags = variants_with_tags.filter(tag__in=self.tag_ids)
        if cutoff:
            variants_with_tags = variants_with_tags.filter(created__gte=cutoff)
        # Tagging is done manually so this will only ever be small - much faster to convert to list
        variant_ids = list(variants_with_tags.values_list("variant_id", flat=True))
        q_list = [Q(pk__in=variant_ids)]

        if self.mode == TagNodeMode.ALL_TAGS:
            tags_qs = VariantTag.filter_for_user(self.analysis.user)
            # We already have tags from this analysis, no need to retrieve again
            tags_qs = tags_qs.exclude(analysis=self.analysis)
            if cutoff:
                tags_qs = tags_qs.filter(created__gte=cutoff)
            # Builds from different analyses (maybe diff builds) - so do query using Allele
            q_list.append(VariantTag.variants_for_build_q(self.analysis.genome_build, tags_qs, self.tag_ids))

        q = reduce(operator.or_, q_list)
        if self.node_input == TagNodeInput.PARENT_NOT_TAGGED:
            q = ~q
        return q

    def get_node_name(self):
        if self.visible:
            description_list = []
            if self.node_input == TagNodeInput.PARENT_NOT_TAGGED:
                description_list.append("Exclude")

            if self.mode == TagNodeMode.ALL_TAGS:
                description_list.append("Global")
            else:
                description_list.append("Analysis")

            if self.tag_ids:
                description_list.append(', '.join(self.tag_ids))
            else:
                description_list.append("Tags")

            if self.tagged_within_days is not None:
                description_list.append(f"≤ {self.tagged_within_days}d")

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

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-tags")

    def _get_method_summary(self):
        summary = f"Tagged {', '.join(self.tag_ids)} ({self.get_mode_display()})"
        if self.tagged_within_days is not None:
            summary += f", tagged within {self.tagged_within_days} days" \
                       f" (since {localtime(self.tagged_within_cutoff):%d %b %Y})"
        return summary

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
        if self.node_input == TagNodeInput.TAGGED_VARIANTS:
            return 0
        return 1

    @staticmethod
    def get_analysis_tags_node(analysis):
        from analysis.tasks.node_update_tasks import update_node_task

        node, created = TagNode.objects.get_or_create(analysis=analysis,
                                                      name=TagNode.ANALYSIS_TAGS_NAME,
                                                      node_input=TagNodeInput.TAGGED_VARIANTS,
                                                      mode=TagNodeMode.THIS_ANALYSIS,
                                                      visible=False)
        if created:
            # Should be fast, so do sync (not as celery job)
            update_node_task(node.pk, node.version)
        return node


class TagNodeTag(NodeAuditLogMixin, models.Model):
    """ Stores multi-select values """
    tag_node = models.ForeignKey(TagNode, on_delete=CASCADE)
    tag = models.ForeignKey(Tag, null=True, blank=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ("tag_node", "tag")

    def _get_node(self):
        return self.tag_node

    def __str__(self):
        return f"TagNodeTag {self.tag_node_id}: {self.tag}"


auditlog.register(TagNode)
auditlog.register(TagNodeTag)
