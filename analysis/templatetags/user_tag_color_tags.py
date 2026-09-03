import json
from collections import defaultdict

from django import template
from django.utils.safestring import mark_safe

from analysis.models import VariantTag
from analysis.models.nodes.node_counts import get_node_count_colors, get_tag_node_count_colors
from annotation.models import AnnotationVersion
from library import tag_utils
from library.django_utils import get_field_counts
from snpdb.models import GenomeBuild
from snpdb.models.models_enums import TagFilter
from snpdb.models.models_user_settings import UserSettings
from snpdb.utils import get_tag_sort_order_by_tag, get_tag_styles_and_colors
from snpdb.variant_queries import get_variant_queryset_for_gene_symbol

register = template.Library()


def render_user_tag_styles(prefix, user_tag_style):
    """ CSS rules for UserTagColor - .<prefix><tag> > .user-tag-colored """
    css_string = ''
    for tag, data in user_tag_style:
        if data:
            data_css_lines = []
            for k, v in data.items():
                data_css_lines.append(f"{k}: {v} !important;")

            data_string = '\n'.join(data_css_lines)
            string = """
    .%s%s>.user-tag-colored {
        %s
    }
            """
            css_string += string % (prefix, tag, data_string)
    return css_string


class VariableCSSRGBNode(template.Node):
    """ Renders CSS rules for UserTagColor """

    def __init__(self, prefix, user_tag_style):
        self.prefix = template.Variable(prefix)
        self.user_tag_style = template.Variable(user_tag_style)

    def render(self, context):
        prefix = self.prefix.resolve(context)
        user_tag_style = self.user_tag_style.resolve(context)
        return render_user_tag_styles(prefix, user_tag_style)


class VariantTagsJSNode(template.Node):

    def __init__(self, nodes):
        self.variable = template.Variable(nodes)

    def render(self, context):
        analysis = self.variable.resolve(context)

        variant_tags = defaultdict(list)
        variant_tags_qs = VariantTag.objects.filter(analysis=analysis).values_list('variant__id', 'tag__id')
        for variant_id, tag_id in variant_tags_qs:
            variant_tags[variant_id].append(tag_id)
        return json.dumps(variant_tags)


@register.tag
def render_rgb_css(_parser, token):
    return VariableCSSRGBNode(*tag_utils.get_passed_objects(token))


@register.simple_tag(takes_context=True)
def render_node_count_colors_css(context):
    """ Legend swatch colours - the built in filters, plus one per tag in the user's tag colours """
    prefix = 'node-count-legend-'
    user = context["user"]
    css = render_user_tag_styles(prefix, get_node_count_colors("background-color"))
    css += render_user_tag_styles(prefix, get_tag_node_count_colors(user, "background-color"))
    return mark_safe(css)


@register.tag
def render_variant_tags_dict(_parser, token):
    return VariantTagsJSNode(tag_utils.get_passed_object(token))


@register.simple_tag(takes_context=True)
def render_variant_tag_order(context):
    """ {tag_id: sort_order} for JS tag sorting - see sortVariantTags in grid.js """
    return mark_safe(json.dumps(get_tag_sort_order_by_tag(context["user"])))


@register.inclusion_tag("analysis/tags/render_tag_styles_and_formatter.html", takes_context=True)
def render_tag_styles_and_formatter(context, tag_colors_collection=None):
    """ Also relies on global.js being included """
    user = context["user"]
    user_tag_styles, _ = get_tag_styles_and_colors(user, tag_colors_collection=tag_colors_collection)

    return {"user_tag_styles": user_tag_styles,
            "url_name_visible": context["url_name_visible"]}


@register.inclusion_tag("analysis/tags/tag_colors_collection_link.html", takes_context=True)
def tag_colors_collection_link(context):
    """ Explains where tag colours/sort order come from, linking to where they can be changed """
    user_settings = UserSettings.get_for_user(context["user"])
    return {"tag_colors_collection": user_settings.tag_colors}


@register.inclusion_tag("analysis/tags/tag_counts_summary.html", takes_context=True)
def tag_counts_summary(context, genome_build: GenomeBuild = None, gene_symbol=None,
                       tag_counts=None, selected=None):
    """ Pill + count toggles that filter the grid below them - the page wires them up with
        setupTagCountsSummary(). tag_counts is a (tag, count) list - pass it in if the page has
        already counted them, it's an expensive count """
    if tag_counts is None:
        tag_kwargs = {}
        if gene_symbol:
            annotation_version = AnnotationVersion.latest(genome_build)
            gene_variant_qs = get_variant_queryset_for_gene_symbol(gene_symbol, annotation_version,
                                                                   traverse_aliases=True)
            tag_kwargs["variant_qs"] = gene_variant_qs
        variant_tags_qs = VariantTag.get_for_build(genome_build=genome_build, **tag_kwargs)
        tag_counts = get_field_counts(variant_tags_qs, "tag").items()

    sort_order_by_tag = get_tag_sort_order_by_tag(context["user"])
    tag_counts = sorted(tag_counts, key=lambda tc: (sort_order_by_tag.get(tc[0], 0), tc[0]))
    return {
        # The label is what the analysis grid takes as its extra_filters @see TagFilter
        "tag_counts": [(tag, TagFilter.label(tag), count) for tag, count in tag_counts],
        "selected": selected or [],
    }


@register.inclusion_tag("analysis/tags/tag_counts.html")
def tag_counts(tag_counts_dict):
    return {"tag_counts": tag_counts_dict}
