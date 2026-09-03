import logging
from typing import Optional

from django.db.models import Count
from django.db.models.query_utils import Q

from analysis.exceptions import NonFatalNodeError
from annotation.models import AnnotationVersion, GeneAnnotation
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity
from classification.enums import ClinicalSignificance, SomaticClinicalSignificance
from classification.models import Classification
from snpdb.models import Tag, Variant
from snpdb.models.models_enums import BuiltInFilters, TagFilter
from snpdb.models.models_user_settings import AbstractNodeCountSettings
from snpdb.utils import get_tag_styles_and_colors

# Tags the user hasn't given a colour still need a visible legend swatch
DEFAULT_TAG_NODE_COUNT_COLOR = "#888888"


def get_omim_q(annotation_version: AnnotationVersion) -> Q:
    """ Genes with OMIM terms, as a subquery against the small GeneAnnotation table - going through
        'variantannotation__gene__geneannotation' joins the huge VariantAnnotation table through Gene """
    gene_annotation_qs = GeneAnnotation.objects.filter(omim_terms__isnull=False)
    if gene_annotation_version_id := annotation_version.gene_annotation_version_id:
        gene_annotation_qs = gene_annotation_qs.filter(version_id=gene_annotation_version_id)
    return Q(variantannotation__gene__in=gene_annotation_qs.values_list("gene_id", flat=True))


def get_analysis_tagged_variant_ids(analysis, tag_ids: list[str]) -> list[int]:
    """ Variants carrying any of tag_ids in this analysis. Tagging is manual so this stays small -
        small enough to pass around as a list, @see TagNode.tagged_variants_q which does the same """
    # Go via Variant rather than VariantTag as analysis.models.models_variant_tag imports AnalysisNode,
    # which imports this module
    tagged_qs = Variant.objects.filter(varianttag__analysis=analysis, varianttag__tag__in=tag_ids)
    return list(tagged_qs.values_list("pk", flat=True).distinct())


def is_extra_filter(extra_filters) -> bool:
    """ Whether extra_filters narrows the node - "default" (and None) mean show the whole thing """
    return bool(TagFilter.get_tag_ids(extra_filters)) or extra_filters in dict(BuiltInFilters.FILTER_CHOICES)


def get_extra_filters_q(analysis, extra_filters) -> Q:
    if tag_ids := TagFilter.get_tag_ids(extra_filters):
        return Q(pk__in=get_analysis_tagged_variant_ids(analysis, tag_ids))

    user = analysis.user
    annotation_version = analysis.annotation_version
    if extra_filters == BuiltInFilters.CLINVAR:
        # ClinVar's 3 classification axes - a variant is significant if any of them says so
        q = Q(clinvar__highest_pathogenicity__gte=ClinVarPathogenicity.LIKELY_PATHOGENIC) \
            | Q(clinvar__somatic_tier__in=SomaticClinicalSignificance.TIER_1_AND_2_VALUES) \
            | Q(clinvar__highest_oncogenicity__gte=ClinVarOncogenicity.LIKELY_ONCOGENIC)
    elif extra_filters == BuiltInFilters.OMIM:
        q = get_omim_q(annotation_version)
    elif extra_filters in [BuiltInFilters.CLASSIFIED, BuiltInFilters.CLASSIFIED_PATHOGENIC,
                           BuiltInFilters.CLASSIFIED_TIER_1_2]:
        clinical_significance_list = None
        somatic_clinical_significance_list = None
        if extra_filters == BuiltInFilters.CLASSIFIED_PATHOGENIC:
            clinical_significance_list = [ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.PATHOGENIC]
        elif extra_filters == BuiltInFilters.CLASSIFIED_TIER_1_2:
            somatic_clinical_significance_list = SomaticClinicalSignificance.TIER_1_AND_2_VALUES
        q = Classification.get_variant_q(user, annotation_version.genome_build, clinical_significance_list,
                                         somatic_clinical_significance_list=somatic_clinical_significance_list)
    elif extra_filters == BuiltInFilters.IMPACT_HIGH_OR_MODERATE:
        q = Q(variantannotation__impact__in=(PathogenicityImpact.HIGH, PathogenicityImpact.MODERATE))
    elif extra_filters == BuiltInFilters.COSMIC:
        q = Q(variantannotation__cosmic_id__isnull=False)
    else:
        logging.warning("get_extra_filters_q, unknown filter '%s'", extra_filters)
        q = Q(pk__isnull=False)  # No op
    return q


def get_node_extra_filters_q(node, extra_filters) -> Optional[Q]:
    """ The node's own filter for an 'extra_filters' selection, or None if it selects everything.
        A global TagNode's tags reach outside the analysis, so it decides its own tag scope -
        everything else is analysis-scoped, which is what the DAG's node counts are """
    if not is_extra_filter(extra_filters):
        return None
    if tag_ids := TagFilter.get_tag_ids(extra_filters):
        if (node_tag_q := node.get_extra_filters_tag_q(tag_ids)) is not None:
            return node_tag_q
    return get_extra_filters_q(node.analysis, extra_filters)


def get_extra_filters_count(node, extra_filters) -> Optional[int]:
    """ How many rows the node shows under extra_filters, or None if we can't say cheaply. There
        being no filter is one of those - the node's own count is what covers that.
        NodeCount is reached through the node version rather than imported - analysis_node imports us """
    if not is_extra_filter(extra_filters):
        return None
    tag_ids = TagFilter.get_tag_ids(extra_filters)
    node_tag_q = node.get_extra_filters_tag_q(tag_ids) if tag_ids else None
    if node_tag_q is None:
        # Stored counts are analysis-scoped, so they only speak for an analysis-scoped filter
        if node_count := node.node_version.nodecount_set.filter(label=extra_filters).first():
            return node_count.count
    if tag_ids:
        # Tags narrow the node to a handful of variants, so counting them exactly here is cheap
        q = node_tag_q if node_tag_q is not None else get_extra_filters_q(node.analysis, extra_filters)
        try:
            return node.get_queryset(inner_query_distinct=True).filter(q).count()
        except NonFatalNodeError:
            pass  # An ancestor isn't ready - the count comes back when the node reloads
    return None


def get_node_count_colors(css_property):
    """ Returns a list of tuples with last element being a dict,
        css_property of "color" = [('ClinVar', {color: #ff0000}), etc] """

    node_count_colors = []
    for label, color in BuiltInFilters.COLORS:
        node_count_colors.append((label, {css_property: color}))

    return node_count_colors


def get_tag_node_count_colors(user, css_property):
    """ Tag node counts are drawn in the user's tag colours - same shape as get_node_count_colors() """
    _, user_tag_colors = get_tag_styles_and_colors(user)
    return [(TagFilter.label(tag_id), {css_property: color or DEFAULT_TAG_NODE_COUNT_COLOR})
            for tag_id, color in user_tag_colors.items()]


def get_node_counts_mine_and_available(analysis):
    """ (mine, available built in filters, available tags) - mine keeps the user's order, the rest stay
        in their natural order. Tags are split out as there are far more of them than built in filters """
    node_count_types = analysis.get_node_count_types()

    my_choices = [x[0] for x in node_count_types]
    filter_choices = [x[0] for x in BuiltInFilters.CHOICES]
    tag_choices = [TagFilter.label(tag_id) for tag_id in Tag.live_qs().order_by("pk").values_list("pk", flat=True)]

    def _node_count_list(choices):
        node_counts_list = []
        for node_count in choices:
            if description := AbstractNodeCountSettings.get_node_count_description(node_count):
                node_counts_list.append({"pk": node_count,
                                         "css_classes": 'node-count-legend-' + node_count,
                                         "tag_id": TagFilter.get_tag_id(node_count),
                                         "description": description})
        return node_counts_list

    mine = set(my_choices)
    return (_node_count_list(my_choices),
            _node_count_list([c for c in filter_choices if c not in mine]),
            _node_count_list([c for c in tag_choices if c not in mine]))


def get_node_counts_and_labels_dict(node, counts_to_get):

    # Need to do inner query as distinct needs to be applied
    # before aggregate functions
    qs = node.get_queryset(inner_query_distinct=True)
    return _aggregate_node_counts(qs, node.analysis, counts_to_get)


def get_tagged_variant_ids_by_label(analysis, tag_labels) -> dict[str, list[int]]:
    """ The tagged variants behind each tag node count label. The same for every node in the analysis,
        so look them up once and pass them to get_tag_node_counts_dict() per node """
    return {label: get_analysis_tagged_variant_ids(analysis, [TagFilter.get_tag_id(label)])
            for label in tag_labels}


def get_tag_node_counts_dict(node, tagged_variant_ids_by_label: dict[str, list[int]]) -> dict[str, int]:
    """ Tags change without bumping node versions, so these are recounted in place. Restricting the
        node queryset to the (small) tagged variants first makes this far cheaper than a node reload """
    all_tagged_variant_ids = set()
    for variant_ids in tagged_variant_ids_by_label.values():
        all_tagged_variant_ids.update(variant_ids)

    if not all_tagged_variant_ids:
        return dict.fromkeys(tagged_variant_ids_by_label, 0)

    qs = node.get_queryset(inner_query_distinct=True).filter(pk__in=all_tagged_variant_ids)
    aggregate_kwargs = {label: Count("pk", filter=Q(pk__in=variant_ids), empty_result_set_value=0)
                        for label, variant_ids in tagged_variant_ids_by_label.items()}
    return _aggregate(qs, aggregate_kwargs)


def _aggregate_node_counts(qs, analysis, counts_to_get) -> dict[str, int]:
    aggregate_kwargs = {}
    for count_type in counts_to_get:
        if count_type == BuiltInFilters.TOTAL:
            q = None
        else:
            q = get_extra_filters_q(analysis, count_type)
        aggregate_kwargs[count_type] = Count("pk", filter=q, empty_result_set_value=0)
    return _aggregate(qs, aggregate_kwargs)


def _aggregate(qs, aggregate_kwargs) -> dict[str, int]:
    # empty_result_set_value=0 only works for Django >= 4, so we handle None manually
    node_counts = qs.aggregate(**aggregate_kwargs)
    return {k: v if v is not None else 0 for k, v in node_counts.items()}
