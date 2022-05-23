import json
from collections import defaultdict

from django import template

from analysis.models import VariantTag
from analysis.models.nodes.node_counts import get_node_count_colors
from annotation.models import AnnotationVersion
from library import tag_utils
from library.django_utils import get_field_counts
from snpdb.models import UserTagColors, GenomeBuild
from snpdb.variant_queries import get_variant_queryset_for_gene_symbol

register = template.Library()


class AbstractCSSRGBNode(template.Node):
    """ Renders CSS rule for UserTagColor """

    def render_user_tag_styles(self, prefix, user_tag_style):
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


class VariableCSSRGBNode(AbstractCSSRGBNode, template.Node):

    def __init__(self, prefix, user_tag_style):
        self.prefix = template.Variable(prefix)
        self.user_tag_style = template.Variable(user_tag_style)

    def render(self, context):
        prefix = self.prefix.resolve(context)
        user_tag_style = self.user_tag_style.resolve(context)
        return self.render_user_tag_styles(prefix, user_tag_style)


class ArgsCSSRGBNode(AbstractCSSRGBNode, template.Node):

    def __init__(self, prefix, user_tag_style):
        self.prefix = prefix
        self.user_tag_style = user_tag_style

    def render(self, context):
        return self.render_user_tag_styles(self.prefix, self.user_tag_style)


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


@register.tag
def render_node_count_colors_css(_parser, _token):
    prefix = 'node-count-legend-'
    tag_rgb = get_node_count_colors("background-color")
    return ArgsCSSRGBNode(prefix, tag_rgb)


@register.tag
def render_variant_tags_dict(_parser, token):
    return VariantTagsJSNode(tag_utils.get_passed_object(token))


@register.inclusion_tag("analysis/tags/render_tag_styles_and_formatter.html", takes_context=True)
def render_tag_styles_and_formatter(context):
    """ Also relies on global.js being included """
    user = context["user"]
    user_tag_styles, _ = UserTagColors.get_tag_styles_and_colors(user)

    return {"user_tag_styles": user_tag_styles,
            "url_name_visible": context["url_name_visible"]}


@register.inclusion_tag("analysis/tags/tag_counts_filter.html", takes_context=True)
def tag_counts_filter(context, genome_build: GenomeBuild,
                      click_func=None, show_all_func=None, gene_symbol=None, any_tag_button=True):
    tag_kwargs = {}
    if gene_symbol:
        annotation_version = AnnotationVersion.latest(genome_build)
        gene_variant_qs = get_variant_queryset_for_gene_symbol(gene_symbol, annotation_version,
                                                               traverse_aliases=True)
        tag_kwargs["variant_qs"] = gene_variant_qs
    variant_tags_qs = VariantTag.get_for_build(genome_build=genome_build, **tag_kwargs)
    tag_counts = sorted(get_field_counts(variant_tags_qs, "tag").items())
    return {
        "any_tag_button": any_tag_button,
        "tag_counts": tag_counts,
        "click_func": click_func,
        "show_all_func": show_all_func,
    }
