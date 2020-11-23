from django.contrib.postgres.aggregates.general import StringAgg
from django.http.response import Http404
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
import re

from analysis import grids
from analysis.views.analysis_permissions import get_node_subclass_or_non_fatal_exception
from analysis.views.node_json_view import NodeJSONGetView
from annotation.vcf_files.vcf_export_utils import grid_export_vcf
from library.constants import WEEK_SECS
from library.jqgrid_export import grid_export_csv
from snpdb.models import Sample
from snpdb.models.models_variant import Variant


@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridHandler(NodeJSONGetView):

    def _get_data(self, request):
        grid = variant_grid_from_request(request)
        return grid.get_data(request)


@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridConfig(NodeJSONGetView):

    def _get_data(self, request, node_id, version_id, extra_filters):
        node = get_node_subclass_or_non_fatal_exception(request.user, node_id, version=version_id)
        errors = node.get_errors()

        if errors:
            ret = {"errors": errors}
        else:
            grid = grids.VariantGrid(request.user, node, extra_filters)
            ret = grid.get_config(as_json=False)
        return ret


def format_items_iterator(analysis, sample_ids, items):
    """ A few things are done in JS formatters, eg changing -1 to missing values (? in grid) and tags
        We can't just add tags via node queryset (in monkey patch func above) as we'll get an issue with
        tacked on zygosity columns etc not being in GROUP BY or aggregate func. So, just patch items via iterator """

    SAMPLE_FIELDS = ["allele_depth", "allele_frequency", "read_depth", "genotype_quality", "phred_likelihood"]

    variant_tags_qs = Variant.objects.filter(varianttag__analysis=analysis)
    variant_tags_qs = variant_tags_qs.annotate(tags=StringAgg("varianttag__tag", delimiter=', ', distinct=True))
    variant_tags = dict(variant_tags_qs.values_list("id", "tags"))

    for item in items:
        for sample_id in sample_ids:
            for f in SAMPLE_FIELDS:
                sample_field = f"sample_{sample_id}_ov_{f}"
                val = item.get(sample_field)
                if val and val == -1:
                    item[sample_field] = "."

        variant_id = item["id"]
        tags = variant_tags.get(variant_id)
        if tags:
            item["tags"] = tags
        yield item


def node_grid_export(request):
    grid = variant_grid_from_request(request, sort_by_contig_and_position=True)

    # TODO: Change filename to use set operation between samples
    basename = f"node_{grid.node.pk}"
    if grid.name:
        name = re.sub(r"\s", "_", grid.name)
        basename = f"{basename}_{name}"

    sample_ids = grid.node.get_sample_ids()
    _, _, items = grid.get_items(request)
    items = format_items_iterator(grid.node.analysis, sample_ids, items)

    export_type = request.GET["export_type"]
    colmodels = grid.get_colmodels()

    if export_type == 'csv':
        return grid_export_csv(basename, colmodels, items)
    elif export_type == 'vcf':
        genome_build = grid.node.analysis.genome_build
        values_qs = Sample.objects.filter(id__in=sample_ids).values_list("id", "name")
        sample_names_by_id = dict(values_qs)
        return grid_export_vcf(basename, genome_build, colmodels, items, sample_ids, sample_names_by_id)
    raise Http404(f"unknown export type: '{export_type}'")


def variant_grid_from_request(request, **kwargs):
    extra_filters = request.GET["extra_filters"]
    node_id = request.GET["node_id"]
    version_id = request.GET["version_id"]
    node = get_node_subclass_or_non_fatal_exception(request.user, node_id, version=version_id)
    return grids.VariantGrid(request.user, node, extra_filters, **kwargs)
