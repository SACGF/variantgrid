import csv
import operator
import re
from collections import Counter
from io import StringIO
from typing import Iterator, Optional
from vcf import Writer, Reader
from vcf.model import _Substitution, _Record, make_calldata_tuple, _Call

from analysis.grids import ExportVariantGrid
from analysis.models import AnalysisNode
from annotation.models import VariantTranscriptAnnotation
from genes.models import CanonicalTranscriptCollection
from library.django_utils import get_model_fields
from library.django_utils.jqgrid_view import grid_export_csv
from library.utils import StashFile
from patients.models_enums import Zygosity
from snpdb.models import Sample, ColumnVCFInfo, VCFInfoTypes
from snpdb.vcf_export_utils import get_vcf_header_from_contigs


def node_grid_get_export_iterator(request, node, export_type, canonical_transcript_collection=None,
                                  variant_tags_dict=None, basename: str = None, grid_kwargs: dict = None) -> tuple[str, Iterator[str]]:

    if grid_kwargs is None:
        grid_kwargs = {}
    else:
        grid_kwargs = grid_kwargs.copy()

    if export_type == 'vcf':
        # TODO: If we use the contig at a time method in ExportVariantGrid we can remove the sort by contig
        grid_kwargs["order_by"] = ["locus__contig__name", "locus__position"]
        grid_kwargs["af_show_in_percent"] = False

    extra_filters = request.GET.get("extra_filters")
    grid = ExportVariantGrid(request.user, node, extra_filters, **grid_kwargs)

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
        file_iterator = grid_export_csv(colmodels, items)
    elif export_type == 'vcf':
        genome_build = node.analysis.genome_build
        values_qs = Sample.objects.filter(id__in=sample_ids).values_list("id", "name")
        sample_names_by_id = dict(values_qs)
        file_iterator = _grid_export_vcf(genome_build, colmodels, items, sample_ids, sample_names_by_id)
    else:
        raise ValueError(f"unknown export type: '{export_type}'")

    filename = f"{basename}.{export_type}"
    return filename, file_iterator

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


def _grid_export_vcf(genome_build, colmodels, items, sample_ids, sample_names_by_id) -> Iterator[str]:
    samples = [sample_names_by_id[s_id] for s_id in sample_ids]

    use_accession = False
    info_dict = _get_colmodel_info_dict(colmodels)
    vcf_template_file = _colmodels_to_vcf_header(genome_build, info_dict, samples, use_accession=use_accession)
    vcf_reader = Reader(vcf_template_file, strict_whitespace=True)

    pseudo_buffer = StashFile()

    vcf_writer = Writer(pseudo_buffer, vcf_reader)
    # Need to pass escapechar
    vcf_writer.writer = csv.writer(
        pseudo_buffer,
        delimiter="\t",
        lineterminator="\n",
        quoting=csv.QUOTE_NONE,
        escapechar="\\",
    )

    def iter_row_writer():
        for obj in items:
            record = _grid_item_to_vcf_record(info_dict, obj, sample_ids, samples, use_accession=use_accession)
            vcf_writer.write_record(record)
            yield pseudo_buffer.value

    return iter_row_writer


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


def _colmodels_to_vcf_header(genome_build, info_dict, samples, use_accession=True):
    """ returns file which contains header """

    header_lines = get_vcf_header_from_contigs(genome_build, info_dict, samples, use_accession=use_accession)
    return StringIO('\n'.join(header_lines))



def _grid_item_to_vcf_record(info_dict, obj, sample_ids, sample_names, use_accession=True):
    if use_accession:
        CHROM = obj.get("locus__contig__refseq_accession", ".")
    else:
        CHROM = obj.get("locus__contig__name", ".")

    POS = obj.get("locus__position", ".")
    ID = obj.get("variantannotation__dbsnp_rs_id")
    REF = obj.get("locus__ref__seq", ".")
    ALT = obj.get("alt__seq", ".")
    QUAL = '.'  # QUAL = obj.get("annotation__quality", ".")
    FILTER = None
    INFO = {}

    INFO_REPLACE = {
        ";": ",:",  # semi-colon used as INFO delimiter
        ",": "|",  # commas forbidden except as list
    }

    for info_id, data in info_dict.items():
        col = data['column__variant_column']
        if val := obj.get(col):
            if isinstance(val, str):
                for old, new in INFO_REPLACE.items():
                    val = val.replace(old, new)
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
            sample_prefix = f"sample_{sample_id}_samples"
            ad = obj[f"{sample_prefix}_allele_depth"]
            zygosity = obj[f"{sample_prefix}_zygosity"]
            gt = Zygosity.get_genotype_from_expanded_zygosity(zygosity)
            dp = obj[f"{sample_prefix}_read_depth"]
            af = obj[f"{sample_prefix}_allele_frequency"]
            # GQ/PL/FT are optional now
            # TODO: Ideally, we'd not write them out
            pl = obj.get(f"{sample_prefix}_phred_likelihood", ".")
            gq = obj.get(f"{sample_prefix}_genotype_quality", ".")
            # These are all number=1
            data_args = {'AD': [ad],
                         'GT': gt,
                         'PL': [pl],
                         'DP': [dp],
                         'GQ': [gq],
                         'AF': [af]}

            data = CallData(**data_args)
            call = _Call(record, sample, data)
            samples.append(call)
            sample_indexes[sample] = i

        record.samples = samples

    return record


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
