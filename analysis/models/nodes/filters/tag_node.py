from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q
from django_extensions.db.models import TimeStampedModel

from analysis.models.nodes.analysis_node import Analysis, AnalysisNode
from snpdb.models import Tag, Variant


class TagNode(AnalysisNode):
    ANALYSIS_TAGS_NAME = "Tagged Variants"
    analysis_wide = models.BooleanField(default=False)
    tag = models.ForeignKey(Tag, null=True, blank=True, on_delete=SET_NULL)

    def modifies_parents(self):
        return True

    def _get_node_q(self) -> Q:
        # If we filter with varianttag__analysis=self.analysis we get a row per tag, we want to get a row per variant
        variants_with_tags = VariantTag.objects.filter(analysis=self.analysis)
        if self.tag:
            variants_with_tags = variants_with_tags.filter(tag=self.tag)
        variants_list = list(variants_with_tags.values_list("variant_id", flat=True))  # Much faster converting to list
        return Q(pk__in=variants_list)

    def get_node_name(self):
        if self.visible:
            if self.tag is not None:
                description = f"Tagged {self.tag.id}"
            else:
                description = "All Tags"
        else:
            description = self.ANALYSIS_TAGS_NAME  # Has to be set to this

        return description

    @staticmethod
    def get_node_class_label():
        return "Tags"

    def _get_method_summary(self):
        return f"Tagged {self.tag.id}"

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.tag:
            css_classes.append(f"tagged-{self.tag.id}")
        return css_classes

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.analysis_wide:
            return 0
        return 1

    @staticmethod
    def get_analysis_tags_node(analysis):
        node, created = TagNode.objects.get_or_create(analysis=analysis,
                                                      name=TagNode.ANALYSIS_TAGS_NAME,
                                                      analysis_wide=True,
                                                      visible=False)
        if created:
            node.load()  # Should be fast, so no need for celery job
        return node


class VariantTag(TimeStampedModel):
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    # Most recent node where it was added
    node = models.ForeignKey(AnalysisNode, null=True, on_delete=SET_NULL)  # Keep even if node deleted
    user = models.ForeignKey(User, on_delete=CASCADE)

    def can_write(self, user):
        return self.analysis.can_write(user)
