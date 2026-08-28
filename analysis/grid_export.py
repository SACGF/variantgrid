import operator
import re
from collections import Counter
from collections.abc import Callable, Iterator
from typing import Optional

from analysis.grids import ExportVariantGrid
from analysis.models import AnalysisNode
from annotation.models import VariantTranscriptAnnotation
from genes.models import CanonicalTranscriptCollection
from library.django_utils import get_model_fields
from library.django_utils.grid_export import EXPORT_ROWS_PER_CHUNK, grid_export_csv
from library.genomics.vcf_writer import VCFWriter
from library.utils import StashFile, iter_fixed_chunks
from patients.models_enums import Zygosity
from snpdb.grid_columns.variant_columns import format_rows
from snpdb.models import Sample, VariantGridColumn
from snpdb.vcf_export_columns import COLUMN_VCF_INFO
from snpdb.vcf_export_utils import get_vcf_header_from_contigs

TRANSCRIPT_REPLACE_BATCH_SIZE = 1000


def node_grid_get_export_iterator(request, node, export_type, canonical_transcript_collection=None,
                                  variant_tags_dict=None, basename: str = None, grid_kwargs: dict = None,
                                  row_wrapper: Callable = None) -> tuple[str, Iterator[str]]:
    """ row_wrapper: wraps the rows iterator, e.g. so a Celery task can report progress per row -
        the file iterator yields a chunk of rows at a time so is no use for counting """

    if grid_kwargs is None:
        grid_kwargs = {}
    else:
        grid_kwargs = grid_kwargs.copy()

    if export_type == 'vcf':
        # ExportVariantGrid emits genome build contig order, which is what the header declares
        grid_kwargs["af_show_in_percent"] = False

    extra_filters = request.GET.get("extra_filters")
    grid = ExportVariantGrid(request, node, extra_filters, **grid_kwargs)

    if basename is None:
        basename = get_node_export_basename(node)
    sample_ids = node.get_sample_ids()
    items = grid.export_rows()

    if canonical_transcript_collection:
        basename += f"_{canonical_transcript_collection}"
        items = _replace_transcripts_iterator(grid, canonical_transcript_collection, items)

    tag_stale_date = node.analysis.variant_tag_stale_date
    items = format_items_iterator(items, variant_tags_dict, tag_stale_date=tag_stale_date)
    if row_wrapper:
        items = row_wrapper(items)

    columns = grid.csv_columns()

    if export_type == 'csv':
        file_iterator = grid_export_csv(columns, items)
    elif export_type == 'vcf':
        genome_build = node.analysis.genome_build
        values_qs = Sample.objects.filter(id__in=sample_ids).values_list("id", "name")
        sample_names_by_id = dict(values_qs)
        file_iterator = _grid_export_vcf(genome_build, columns, items, sample_ids, sample_names_by_id)
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
    name_parts.append(f"v{node.version}")
    # A live-source node's rows depend on data that changes under it - say which version this file reflects
    if live_data_sources := node.live_data_sources:
        sources = "_".join(f"{k}.{v}" for k, v in sorted(live_data_sources.items()))
        name_parts.append(re.sub(r"\W", "_", sources))
    return "_".join(name_parts)


def _grid_export_vcf(genome_build, columns, items, sample_ids, sample_names_by_id) -> Iterator[str]:
    samples = [sample_names_by_id[s_id] for s_id in sample_ids]

    use_accession = False
    info_dict = _get_column_info_dict(columns)
    header_lines = get_vcf_header_from_contigs(genome_build, info_dict, samples, use_accession=use_accession)

    pseudo_buffer = StashFile()
    writer = VCFWriter(pseudo_buffer, header_lines)
    yield pseudo_buffer.value  # header

    for i, obj in enumerate(items, start=1):
        chrom, pos, vcf_id, ref, alt, info, fmt, sample_calls = \
            _grid_item_to_vcf_row(info_dict, obj, sample_ids, samples, use_accession=use_accession)
        writer.write_record(chrom, pos, ref, alt, vcf_id=vcf_id, info=info, fmt=fmt, sample_calls=sample_calls)
        if i % EXPORT_ROWS_PER_CHUNK == 0:
            yield pseudo_buffer.value
    if remaining := pseudo_buffer.value:
        yield remaining


def _get_column_vcf_info():
    columns = [c.column for c in COLUMN_VCF_INFO]
    variant_column_by_name = dict(
        VariantGridColumn.objects.filter(pk__in=columns).values_list("grid_column_name", "variant_column")
    )
    column_vcf_info = {}
    for c in COLUMN_VCF_INFO:
        if variant_column := variant_column_by_name.get(c.column):
            column_vcf_info[variant_column] = {
                "column__variant_column": variant_column,
                "info_id": c.info_id,
                "number": c.number,
                "type": c.type.label,
                "description": c.description,
            }
    return column_vcf_info


def _get_column_info_dict(columns):
    column_vcf_info = _get_column_vcf_info()

    info_dict = {}
    for c in columns:
        name = c['name']
        col_info = column_vcf_info.get(name)
        if col_info:
            col_info['number'] = col_info['number'] or '.'

            info_id = col_info['info_id']
            info_dict[info_id] = col_info
    return info_dict


VCF_INFO_REPLACE = {
    ";": ",:",  # semi-colon used as INFO delimiter
    ",": "|",  # commas forbidden except as list
}

VCF_SAMPLE_FORMAT = ['GT', 'AD', 'AF', 'PL', 'DP', 'GQ']


def _vcf_info_encode(val):
    if isinstance(val, str):
        for old, new in VCF_INFO_REPLACE.items():
            val = val.replace(old, new)
    return val


def _format_sample_value(value):
    """ Mirror how the value is rendered in a sample column ('.' for missing) """
    if value is None:
        return "."
    return str(value)


def _format_sample_call(gt, ad, af, pl, dp, gq) -> str:
    # GT leads whenever present; the remaining fields follow VCF_SAMPLE_FORMAT order
    parts = [gt] if gt else []
    parts.extend(_format_sample_value(v) for v in (ad, af, pl, dp, gq))
    return ":".join(parts)


def _grid_item_to_vcf_row(info_dict, obj, sample_ids, sample_names, use_accession=True):
    if use_accession:
        chrom = obj.get("locus__contig__refseq_accession", ".")
    else:
        chrom = obj.get("locus__contig__name", ".")

    pos = obj.get("locus__position", ".")
    vcf_id = obj.get("variantannotation__dbsnp_rs_id")
    if vcf_id is None:
        vcf_id = "."
    ref = obj.get("locus__ref__seq", ".")
    alt = obj.get("alt__seq", ".")

    info = {}
    for info_id, data in info_dict.items():
        col = data['column__variant_column']
        if val := obj.get(col):
            info[info_id] = _vcf_info_encode(val)

    fmt = None
    sample_calls = None
    if sample_ids:
        fmt = ':'.join(VCF_SAMPLE_FORMAT)
        sample_calls = []
        for sample_id in sample_ids:
            sample_prefix = f"sample_{sample_id}_samples"
            ad = obj[f"{sample_prefix}_allele_depth"]
            zygosity = obj[f"{sample_prefix}_zygosity"]
            gt = Zygosity.get_genotype_from_expanded_zygosity(zygosity)
            dp = obj[f"{sample_prefix}_read_depth"]
            af = obj[f"{sample_prefix}_allele_frequency"]
            # GQ/PL/FT are optional now
            pl = obj.get(f"{sample_prefix}_phred_likelihood", ".")
            gq = obj.get(f"{sample_prefix}_genotype_quality", ".")
            sample_calls.append(_format_sample_call(gt, ad, af, pl, dp, gq))

    return chrom, pos, vcf_id, ref, alt, info or None, fmt, sample_calls


def _summarise_tags_global(tags_global: str, stale_cutoff: Optional[str]) -> str:
    """ tags_global entries are 'tag:date' - see get_variantgrid_extra_annotate.
        stale_cutoff: ISO date - events on/after it count as fresh (None = no fresh counts) """
    totals = Counter()
    fresh = Counter()
    for entry in tags_global.split("|"):
        tag, sep, entry_date = entry.rpartition(":")
        if not sep:  # rpartition puts a separator-less entry in the tail
            tag, entry_date = entry_date, ""
        totals[tag] += 1
        if entry_date and (stale_cutoff is None or entry_date >= stale_cutoff):
            fresh[tag] += 1

    summarised_tags = []
    for tag, count in sorted(totals.items(), key=operator.itemgetter(1), reverse=True):
        summary = f"{tag} x {count}" if count > 1 else tag
        if stale_cutoff is not None and (count > 1 or fresh[tag] < count):
            summary += f" ({fresh[tag]} fresh)"
        summarised_tags.append(summary)
    return ", ".join(summarised_tags)


def format_items_iterator(items, variant_tags_dict: Optional[dict] = None, tag_stale_date=None):
    """ A few things are done in JS formatters, e.g. tags
        We can't just add tags via node queryset (in monkey patch func above) as we'll get an issue with
        tacked on zygosity columns etc not being in GROUP BY or aggregate func. So, just patch items via iterator

        variant_tags_dict: key = variant_id, value = tags (for this analysis)
        tag_stale_date: the analysis' variant_tag_stale_date - when set, tags_global gains fresh counts,
                        e.g. 'Artefact x 5 (2 fresh)' """

    if variant_tags_dict is None:
        variant_tags_dict = {}

    stale_cutoff = tag_stale_date.date().isoformat() if tag_stale_date else None
    for item in items:
        if tags_global := item["tags_global"]:
            item["tags_global"] = _summarise_tags_global(tags_global, stale_cutoff)

        variant_id = item["id"]
        if tags := variant_tags_dict.get(variant_id):
            item["tags"] = tags
        yield item


def _replace_transcripts_iterator(grid, ctc: CanonicalTranscriptCollection, items):
    """ Overwrite the representative ('pick') variantannotation__ values with the canonical transcript ones.

        Transcript rows are looked up in batches of variant IDs taken off the streamed items, so we hit the
        (version, variant, transcript_version) index directly rather than re-running the node queryset as an
        'IN' subquery, and only hold one batch in RAM """

    variant_transcript_annotation_variant_id_field = "variant_id"

    # Work out what fields
    transcript_replace_fields = {variant_transcript_annotation_variant_id_field: "id"}

    transcript_fields = set(get_model_fields(VariantTranscriptAnnotation, ignore_fields=["id", "version", "variant"]))
    annotation_prefix = "variantannotation__"
    annotation_prefix_len = len(annotation_prefix)
    for f in [rc.name for rc in grid.enabled_columns]:
        if f.startswith(annotation_prefix):
            suffix = f[annotation_prefix_len:]
            tf = suffix.split("__", 1)[0]
            if tf in transcript_fields:
                transcript_replace_fields[suffix] = f

    # We only need things from VariantTranscriptAnnotation - so join there directly
    version = grid.node.analysis.annotation_version.variant_annotation_version
    ct_qs = ctc.canonicaltranscript_set
    transcript_versions = ct_qs.values_list("transcript_version", flat=True)

    def _transcript_items_by_id(variant_ids: list[int]) -> dict[int, dict]:
        vta_qs = VariantTranscriptAnnotation.objects.filter(version=version, variant__in=variant_ids,
                                                           transcript_version__in=transcript_versions)

        def transcript_items():
            for transcript_data in vta_qs.values(*transcript_replace_fields.keys()):
                transcript_item = {}
                for before, after in transcript_replace_fields.items():
                    transcript_item[after] = transcript_data[before]
                yield transcript_item

        return {item["id"]: item for item in format_rows(grid.enabled_columns, transcript_items())}

    # Loop through items and changeroo
    for batch in iter_fixed_chunks(items, TRANSCRIPT_REPLACE_BATCH_SIZE):
        transcript_items_by_id = _transcript_items_by_id([item["id"] for item in batch])
        for item in batch:
            if transcript_data := transcript_items_by_id.get(item["id"]):
                item.update(transcript_data)
            yield item
