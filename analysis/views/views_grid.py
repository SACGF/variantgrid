import logging
import operator
import os
import re
import time
from collections import Counter
from io import StringIO
from typing import Optional
from urllib.parse import urlencode

from django.contrib.postgres.aggregates.general import StringAgg
from django.core.cache import cache
from django.http.response import Http404, StreamingHttpResponse, HttpResponseRedirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from vcf import Writer, Reader
from vcf.model import _Substitution, _Record, make_calldata_tuple, _Call

from analysis import grids
from analysis.analysis_templates import get_cohort_analysis, get_sample_analysis
from analysis.models import AnalysisNode, AnalysisTemplate, SampleNode
from analysis.views.analysis_permissions import get_node_subclass_or_non_fatal_exception
from analysis.views.node_json_view import NodeJSONGetView, NodeJSONViewMixin
from annotation.models import VariantTranscriptAnnotation
from genes.models import CanonicalTranscriptCollection
from library.constants import WEEK_SECS
from library.django_utils import get_model_fields
from library.jqgrid.jqgrid_export import grid_export_csv, StashFile
from library.utils import name_from_filename
from snpdb.models import Sample, ColumnVCFInfo, VCFInfoTypes, Zygosity, VCF, Cohort
from snpdb.models.models_variant import Variant
from snpdb.vcf_export_utils import get_vcf_header_from_contigs


@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridHandler(NodeJSONViewMixin):
    def get(self, request, *args, **kwargs):
        """ This can be a really expensive operation (ie a few mins)
            And users can sometimes click multiple times, causing the DB to get slow running the same query
            multiple times, interfering with itself - so make a per-user lock, and redirect any further calls
            which should hopefully hit the cache next time
        """
        LOCK_EXPIRE = 60 * 10  # 10 mins
        url = f"{request.path}?" + urlencode(request.GET.dict())
        lock_id = f"{url}_{request.user}"
        if cache.add(lock_id, "true", LOCK_EXPIRE):  # Acquire lock
            try:
                logging.info("Got the lock...")
                response = self.get_response(request, *args, **kwargs)
            finally:
                cache.delete(lock_id)  # release lock
        else:
            logging.info("Don't have lock - going to sleep then retry...")
            time.sleep(2)
            response = HttpResponseRedirect(url)
        return response

    def _get_node(self, request, **kwargs) -> AnalysisNode:
        return _node_from_request(request)

    def _get_redirect(self, request, node):
        """ If node can be represented by another, use that to re-use cache """
        ret = None
        grid_node_id, grid_node_version = node.get_grid_node_id_and_version()
        if grid_node_id != node.pk:
            url = reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id})
            params = request.GET.dict()
            params.update({"node_id": grid_node_id, "version_id": grid_node_version})
            url += "?" + urlencode(params)
            ret = HttpResponseRedirect(url)
        return ret

    def _get_data(self, request, node, **kwargs):
        grid = _variant_grid_from_request(request, node)
        return grid.get_data(request)


@method_decorator([cache_page(WEEK_SECS), vary_on_cookie], name='get')
class NodeGridConfig(NodeJSONGetView):
    def _get_node(self, request, **kwargs) -> AnalysisNode:
        return get_node_subclass_or_non_fatal_exception(request.user, kwargs["node_id"], version=kwargs["node_version"])

    def _get_data(self, request, node, **kwargs):
        errors = node.get_errors(flat=True)

        if errors:
            ret = {"errors": errors}
        else:
            grid = grids.VariantGrid(request.user, node, kwargs["extra_filters"])
            ret = grid.get_config(as_json=False)
        return ret


def format_items_iterator(sample_ids, items, variant_tags_dict: Optional[dict] = None):
    """ A few things are done in JS formatters, e.g. changing -1 to missing values (? in grid) and tags
        We can't just add tags via node queryset (in monkey patch func above) as we'll get an issue with
        tacked on zygosity columns etc not being in GROUP BY or aggregate func. So, just patch items via iterator

        variant_tags_dict: key = variant_id, value = tags (for this analysis) """

    if variant_tags_dict is None:
        variant_tags_dict = {}

    SAMPLE_FIELDS = ["allele_depth", "allele_frequency", "read_depth", "genotype_quality", "phred_likelihood"]

    for item in items:
        for sample_id in sample_ids:
            for f in SAMPLE_FIELDS:
                sample_field = f"sample_{sample_id}_ov_{f}"
                val = item.get(sample_field)
                if val and val == -1:
                    item[sample_field] = "."

        if tags_global := item["tags_global"]:
            tag_counts = Counter(tags_global.split("|"))
            summarised_tags = []
            for tag, count in sorted(tag_counts.items(), key=operator.itemgetter(1), reverse=True):
                if count > 1:
                    summarised_tags.append(f"{tag} x {count}")
                else:
                    summarised_tags.append(tag)

            item["tags_global"] = ", ".join(summarised_tags)

        variant_id = item["id"]
        if tags := variant_tags_dict.get(variant_id):
            item["tags"] = tags
        yield item


def _replace_transcripts_iterator(request, grid, ctc: CanonicalTranscriptCollection, items):
    """ This uses a large amount of RAM - reading a whole  """

    variant_transcript_annotation_variant_id_field = "variant_id"

    # Work out what fields
    transcript_replace_fields = {variant_transcript_annotation_variant_id_field: "id"}

    transcript_fields = set(get_model_fields(VariantTranscriptAnnotation, ignore_fields=["id", "version", "variant"]))
    annotation_prefix = "variantannotation__"
    annotation_prefix_len = len(annotation_prefix)
    for f in grid.get_field_names():
        if f.startswith(annotation_prefix):
            suffix = f[annotation_prefix_len:]
            tf = suffix.split("__", 1)[0]
            if tf in transcript_fields:
                transcript_replace_fields[suffix] = f

    # We only need things from VariantTranscriptAnnotation - so join there directly
    variants_qs = grid.get_queryset(request).values_list("id")
    version = grid.node.analysis.annotation_version.variant_annotation_version
    ct_qs = ctc.canonicaltranscript_set
    transcript_versions = ct_qs.values_list("transcript_version", flat=True)
    vta_qs = VariantTranscriptAnnotation.objects.filter(version=version, variant__in=variants_qs,
                                                        transcript_version__in=transcript_versions)
    transcript_values = vta_qs.values(*transcript_replace_fields.keys())

    # Read into a massive dictionary
    transcript_items_by_id = {}

    def transcript_items():
        for transcript_data in transcript_values:
            transcript_item = {}
            for before, after in transcript_replace_fields.items():
                transcript_item[after] = transcript_data[before]
            yield transcript_item

    for item in grid.iter_format_items(transcript_items()):
        transcript_items_by_id[item["id"]] = item

    # Loop through items and changeroo
    for item in items:
        variant_id = item["id"]
        if transcript_data := transcript_items_by_id.get(variant_id):
            item.update(transcript_data)
        yield item


def get_node_export_basename(node: AnalysisNode) -> str:
    """ For CSV/VCF etc """
    name_parts = []
    if samples := node.get_samples():
        if len(samples) == 1:
            name_parts.append(samples[0].name)

    name_parts.append(f"analysis_{node.analysis.pk}")

    node_label = node.get_node_class_label()
    if not node_label.endswith("Node"):
        node_label += "Node"
    name_parts.append(node_label)
    name_parts.append(str(node.pk))

    if node.name:
        name_underscores = re.sub(r"\s", "_", node.name)
        name_parts.append(name_underscores)
    return "_".join(name_parts)


def cohort_grid_export(request, cohort_id, export_type):
    EXPORT_TYPES = {"csv", "vcf"}

    cohort = Cohort.get_for_user(request.user, cohort_id)
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")

    analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT")
    analysis = get_cohort_analysis(cohort, analysis_template)
    node = analysis.analysisnode_set.get_subclass(output_node=True)  # Should only be 1
    basename = "_".join([name_from_filename(cohort.name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(cohort.genome_build)])
    return _node_grid_export(request, node, export_type, basename=basename, grid_kwargs={"paging": False})


def sample_grid_export(request, sample_id, export_type):
    EXPORT_TYPES = {"csv", "vcf"}

    sample = Sample.get_for_user(request.user, sample_id)
    if export_type not in EXPORT_TYPES:
        raise ValueError(f"{export_type} must be one of: {EXPORT_TYPES}")

    analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_SAMPLE")
    analysis = get_sample_analysis(sample, analysis_template)
    node = SampleNode.objects.get(analysis=analysis, output_node=True)  # Should only be 1
    basename = "_".join([name_from_filename(sample.name), "annotated", f"v{analysis.annotation_version.pk}",
                         str(sample.genome_build)])
    return _node_grid_export(request, node, export_type, basename=basename, grid_kwargs={"paging": False})



def node_grid_export(request, analysis_id):
    export_type = request.GET["export_type"]
    use_canonical_transcripts = request.GET.get("use_canonical_transcripts")

    node = _node_from_request(request)
    canonical_transcript_collection = None
    if use_canonical_transcripts:
        # Whether to use it or not is set server-side. Just use client to see what they wanted
        if ctc := node.analysis.canonical_transcript_collection:
            canonical_transcript_collection = ctc
        else:
            logging.warning("Grid request had 'use_canonical_transcripts' but analysis did not.")

    variant_tags_qs = Variant.objects.filter(varianttag__analysis=node.analysis)
    variant_tags_qs = variant_tags_qs.annotate(tags=StringAgg("varianttag__tag", delimiter=', ', distinct=True))
    variant_tags_dict = dict(variant_tags_qs.values_list("id", "tags"))

    return _node_grid_export(request, node, export_type, canonical_transcript_collection, variant_tags_dict)


def _node_grid_export(request, node, export_type, canonical_transcript_collection=None, variant_tags_dict=None,
                      basename: str = None, grid_kwargs: dict = None):

    if grid_kwargs is None:
        grid_kwargs = {}
    else:
        grid_kwargs = grid_kwargs.copy()

    if export_type == 'vcf':
        grid_kwargs["sort_by_contig_and_position"] = True
        grid_kwargs["af_show_in_percent"] = False

    grid = _variant_grid_from_request(request, node, **grid_kwargs)

    if basename is None:
        basename = get_node_export_basename(node)
    sample_ids = node.get_sample_ids()
    _, _, items = grid.get_items(request)

    if canonical_transcript_collection:
        basename += f"_{canonical_transcript_collection}"
        items = _replace_transcripts_iterator(request, grid, canonical_transcript_collection, items)

    items = format_items_iterator(sample_ids, items, variant_tags_dict)

    colmodels = grid.get_colmodels()

    if export_type == 'csv':
        return grid_export_csv(basename, colmodels, items)
    elif export_type == 'vcf':
        genome_build = node.analysis.genome_build
        values_qs = Sample.objects.filter(id__in=sample_ids).values_list("id", "name")
        sample_names_by_id = dict(values_qs)
        return _grid_export_vcf(basename, genome_build, colmodels, items, sample_ids, sample_names_by_id)
    raise Http404(f"unknown export type: '{export_type}'")


def _node_from_request(request) -> AnalysisNode:
    node_id = request.GET["node_id"]
    version_id = int(request.GET["version_id"])
    return get_node_subclass_or_non_fatal_exception(request.user, node_id, version=version_id)


def _variant_grid_from_request(request, node: AnalysisNode, **kwargs):
    extra_filters = request.GET.get("extra_filters")
    return grids.VariantGrid(request.user, node, extra_filters, **kwargs)


def _grid_export_vcf(filename, genome_build, colmodels, items, sample_ids, sample_names_by_id):
    samples = [sample_names_by_id[s_id] for s_id in sample_ids]

    info_dict = _get_colmodel_info_dict(colmodels)
    vcf_template_file = _colmodels_to_vcf_header(genome_build, info_dict, samples)
    vcf_reader = Reader(vcf_template_file, strict_whitespace=True)

    pseudo_buffer = StashFile()

    vcf_writer = Writer(pseudo_buffer, vcf_reader)

    def iter_row_writer():

        for obj in items:
            record = _grid_item_to_vcf_record(info_dict, obj, sample_ids, samples)
            vcf_writer.write_record(record)
            yield pseudo_buffer.value

    response = StreamingHttpResponse(iter_row_writer(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{filename}.vcf"'
    return response


def _grid_item_to_vcf_record(info_dict, obj, sample_ids, sample_names):  # , get_genotype_from_expanded_zygosity):
    CHROM = obj.get("locus__contig__name", ".")
    POS = obj.get("locus__position", ".")
    ID = obj.get("variantannotation__dbsnp_rs_id")
    REF = obj.get("locus__ref__seq", ".")
    ALT = obj.get("alt__seq", ".")
    QUAL = '.'  # QUAL = obj.get("annotation__quality", ".")
    FILTER = None
    INFO = {}

    for info_id, data in info_dict.items():
        col = data['column__variant_column']
        val = obj.get(col)
        if val:
            INFO[info_id] = val

    FORMAT = None
    MY_FORMAT = ['GT', 'AD', 'AF', 'PL', 'DP', 'GQ']
    CallData = make_calldata_tuple(MY_FORMAT)
    sample_indexes = {}
    samples = []

    if sample_ids:
        FORMAT = ':'.join(MY_FORMAT)

    alts = [_Substitution(ALT)]
    ALT = alts
    record = _Record(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes)

    if sample_ids:
        for i, (sample_id, sample) in enumerate(zip(sample_ids, sample_names)):
            ad = obj[f"{sample_id}_samples_allele_depth"]
            zygosity = obj[f"{sample_id}_samples_zygosity"]
            gt = Zygosity.get_genotype_from_expanded_zygosity(zygosity)
            dp = obj[f"{sample_id}_samples_read_depth"]
            af = obj[f"{sample_id}_samples_allele_frequency"]
            # GQ/PL/FT are optional now
            # TODO: Ideally, we'd not write them out
            pl = obj.get(f"{sample_id}_samples_phred_likelihood", ".")
            gq = obj.get(f"{sample_id}_samples_genotype_quality", ".")
            # TODO: Need to grab information for reference base to be able to properly fill in this data.
            data_args = {'AD': ['.', ad],
                         'GT': gt,
                         'PL': ['.', pl],
                         'DP': ['.', dp],
                         'GQ': ['.', gq],
                         'AF': ['.', af]}

            data = CallData(**data_args)
            call = _Call(record, sample, data)
            samples.append(call)
            sample_indexes[sample] = i

        record.samples = samples

    return record


def _get_column_vcf_info():
    column_vcf_info = {}
    for cvi in ColumnVCFInfo.objects.all().values('column__variant_column', 'info_id', 'number', 'type', 'description'):
        cvi['type'] = VCFInfoTypes(cvi['type']).label
        index = cvi['column__variant_column']
        column_vcf_info[index] = cvi
    return column_vcf_info


def _get_colmodel_info_dict(colmodels):
    column_vcf_info = _get_column_vcf_info()

    info_dict = {}
    for c in colmodels:
        name = c['name']
        col_info = column_vcf_info.get(name)
        if col_info:
            col_info['number'] = col_info['number'] or '.'

            info_id = col_info['info_id']
            info_dict[info_id] = col_info
    return info_dict


def _colmodels_to_vcf_header(genome_build, info_dict, samples):
    """ returns file which contains header """

    header_lines = get_vcf_header_from_contigs(genome_build, info_dict, samples)
    return StringIO('\n'.join(header_lines))
