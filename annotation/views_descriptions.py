import json
from collections import defaultdict

from django.conf import settings
from django.shortcuts import render
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from annotation import annotsv_columns, vep_columns
from annotation.vep_annotation import VEPConfig
from library.constants import WEEK_SECS
from library.unit_percent import get_allele_frequency_formatter
from snpdb.grid_columns.composite_examples import (
    COMPOSITE_EXAMPLE_ROWS,
    SAMPLE_EXAMPLE_PREFIX,
    SAMPLE_EXAMPLE_ROW,
    SAMPLE_EXAMPLE_SAMPLE_NAME,
)
from snpdb.grid_columns.custom_columns import composite_rich_column
from snpdb.grid_columns.grid_sample_columns import SAMPLE_SORT_KEY_LABELS
from snpdb.grids import (
    AF_UNIT_COLUMNS,
    AbstractVariantGrid,
    get_standard_overrides,
    variant_grid_client_extra,
)
from snpdb.models import ColumnAnnotationLevel, UserSettings, VariantGridColumn
from snpdb.views.datatable_view import CellData, RichColumn, rich_column_json

# Composite sections lead the page in the order the level cards below them run
_COMPOSITE_LEVEL_ORDER = [
    None,  # the Variant cell - the identity columns have no annotation level of their own
    ColumnAnnotationLevel.SAMPLE_LEVEL,
    ColumnAnnotationLevel.DATABASE_LEVEL,
    ColumnAnnotationLevel.TRANSCRIPT_LEVEL,
    ColumnAnnotationLevel.VARIANT_LEVEL,
    ColumnAnnotationLevel.HGNC_LEVEL,
    ColumnAnnotationLevel.UNIPROT_LEVEL,
    ColumnAnnotationLevel.CLINVAR_LEVEL,
    ColumnAnnotationLevel.GENE_LEVEL,
]

# The sample call cell is built per sample by the analysis grid rather than catalogued, so its member
# table is written out here - (format column, label, description), zygosity first as the headline
# @see VariantGrid._get_grid_genotype_columns
_SAMPLE_CELL_DESCRIPTION = (
    "Everything the VCF says about this sample&#146;s call, in the one cell - the zygosity glyph, the "
    "allele frequency and depths, quality marks for GQ/PL and a warning where the call failed a filter"
)
_SAMPLE_MEMBER_COLUMNS = [
    ("samples_zygosity", "Sample Zygosity", "GT",
     "Zygosity call (HET/HOM_ALT/REF or Unknown)"),
    ("samples_allele_depth", "AD (Allele Depth)", "AD",
     "Unfiltered allele depth, the number of reads that support each of the reported alleles. All reads "
     "at the position (including reads that did not pass the variant caller&#146;s filters) are included "
     "in this number, except reads that were considered uninformative. Reads are considered uninformative "
     "when they do not provide enough statistical evidence to support one allele over another."),
    ("samples_allele_frequency", "AF (Allele Frequency)", "AF",
     'Unit value (0-1) of Allele Frequency. From "AF" FORMAT field if present in VCF otherwise calculated '
     'for each allele as <code>AD / (sum of all AD for that locus)</code>. See "VCF Info" tab on VCF page '
     'to see which was used for a particular file.'),
    ("samples_read_depth", "DP (Depth)", "DP",
     "Filtered depth, at the sample level. This gives you the number of filtered reads that support each "
     "of the reported alleles. You can check the variant caller&#146;s documentation to see which filters "
     "are applied by default. Only reads that passed the variant caller&#146;s filters are included in this "
     "number. However, unlike the AD calculation, uninformative reads are included in DP."),
    ("samples_genotype_quality", "GQ (Genotype Quality)", "GQ",
     "Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype."),
    ("samples_phred_likelihood", "PL (Phred Likelihood)", "PL",
     "Normalised, Phred-scaled likelihoods for genotypes as defined in the VCF specification"),
    ("samples_filters", "FT (Sample Filters)", "FT",
     "Which of the VCF&#146;s filters this sample&#146;s call failed, where the caller reports them per "
     "sample rather than per record. The cell only marks a call that failed one"),
]


def _sample_composite_section(format_afs) -> dict:
    """ The sample call cell. The analysis grid builds one of these per sample rather than reading it
        from the catalogue, so its column, member table and example row are all written out here """
    column = VariantGridColumn(grid_column_name="sample_call", variant_column="sample_call",
                               annotation_level=ColumnAnnotationLevel.SAMPLE_LEVEL,
                               label="Sample call", description=_SAMPLE_CELL_DESCRIPTION,
                               model_field=False, queryset_field=False)
    rich_column = RichColumn(
        key=None, name=f"{SAMPLE_EXAMPLE_PREFIX}samples_zygosity", width=140,
        label=f"{SAMPLE_EXAMPLE_SAMPLE_NAME} Zygosity", header_title=_SAMPLE_CELL_DESCRIPTION,
        client_renderer="VariantGridFormat.sampleZygosity",
        client_renderer_kwargs={"samplePrefix": SAMPLE_EXAMPLE_PREFIX},
        sort_menu=[{"label": label, "column": f"{SAMPLE_EXAMPLE_PREFIX}{c}"}
                   for c, label in SAMPLE_SORT_KEY_LABELS.items()],
    )
    rows = []
    for i, (format_column, label, source_field, description) in enumerate(_SAMPLE_MEMBER_COLUMNS):
        member_column = VariantGridColumn(grid_column_name=format_column, variant_column=format_column,
                                          annotation_level=ColumnAnnotationLevel.SAMPLE_LEVEL,
                                          label=label, description=description)
        rows.append({
            "column": member_column,
            "source": "VCF",
            "source_detail": "",
            "category": None,
            "source_field": source_field,
            "source_field_processing_description": None,
            "vep": None,
            "level": ColumnAnnotationLevel.SAMPLE_LEVEL.label,
            "role": "headline" if i == 0 else "sort",
        })
    row = format_afs(dict(SAMPLE_EXAMPLE_ROW), [f"{SAMPLE_EXAMPLE_PREFIX}samples_allele_frequency"])
    return {
        "column": column,
        "level": ColumnAnnotationLevel.SAMPLE_LEVEL.label,
        "level_css": ColumnAnnotationLevel.SAMPLE_LEVEL.label.lower(),
        "rows": rows,
        "sort_by": rows[0]["column"].label,
        "column_json": json.dumps(rich_column_json(rich_column, AbstractVariantGrid.default_column_width)),
        "row_json": json.dumps(row),
    }


def _member_role(member) -> str:
    """ What a member contributes to its cell - what it reads as, an alternative sort key, or hover detail """
    if member.sort_order == 0:
        return "headline"
    return "sort" if member.in_sort_menu else "hover"


@cache_page(WEEK_SECS)
@vary_on_cookie  # the information isn't actually different per user, but hack to avoid showing other user's email/notifications etc in the top right
def view_annotation_descriptions(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    variantgrid_columns_by_annotation_level = defaultdict(list)
    annotation_run_levels = [ColumnAnnotationLevel.TRANSCRIPT_LEVEL, ColumnAnnotationLevel.VARIANT_LEVEL]
    columns_by_annotation_level = {al.label: [] for al in annotation_run_levels}

    # Use a VEPConfig so we hide columns whose data file isn't configured for this build
    # (e.g. dbNSFP under T2T) - same source of truth as get_vep_command/BulkVEPVCFAnnotationInserter.
    vep_config = VEPConfig(genome_build)
    vep_referenced = vep_columns.all_variant_grid_column_ids()

    def _first_for_build(vgc_id, build_name):
        for c in vep_columns.for_variant_grid_column(vgc_id, vep_config=vep_config):
            # cosmic_version so cosmic_count describes the INFO field the installed COSMIC release
            # actually carries the count in (#1673)
            if c.applies_to(genome_build_name=build_name, cosmic_version=vep_config.cosmic_version):
                return c
        return None

    def _vep_source_detail(vep) -> str:
        details = []
        if vep.vep_plugin:
            details.append(f"Plugin: {vep.vep_plugin.label}")
        if vep.vep_custom:
            details.append("custom")
        return ", ".join(details)

    def _column_row(vgc):
        """ Where a column's value comes from - None for a VEP column this build doesn't populate """
        if vep := _first_for_build(vgc.pk, genome_build.name):
            row = {
                "source": "VEP",
                "source_detail": _vep_source_detail(vep),
                "category": vep.category,
                "source_field": vep.source_field,
                "source_field_processing_description": vep.source_field_processing_description,
                "vep": vep,
            }
        elif vgc.pk in vep_referenced:
            return None  # VEP column that isn't populated in this build
        elif annotsv := annotsv_columns.for_variant_grid_column(vgc.pk):
            row = {
                "source": "AnnotSV",
                "source_detail": "",
                "category": annotsv.category,
                "source_field": annotsv.source_field,
                "source_field_processing_description": annotsv.source_field_processing_description,
                "vep": None,
            }
        else:
            # Calculated during annotation import rather than read from an annotation tool's output
            row = {
                "source": "VariantGrid",
                "source_detail": "",
                "category": None,
                "source_field": None,
                "source_field_processing_description": None,
                "vep": None,
            }
        row["column"] = vgc
        row["level"] = vgc.get_annotation_level_display()
        return row

    # The grid formats unit AFs server side so its CSV matches what it draws - the example rows hold
    # unit values and go through the same formatter, so they follow the deployment's percent setting
    render_unit_af = get_allele_frequency_formatter(
        source_in_percent=False, dest_in_percent=settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT)

    def _format_afs(row: dict, af_columns) -> dict:
        for path in af_columns:
            if path in row:
                row[path] = render_unit_af(CellData(all_data=row, key=path))
        return row

    column_overrides = get_standard_overrides(settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT)
    variant_grid_columns = list(VariantGridColumn.objects.prefetch_related("composite_members__column"))
    VariantGridColumn.annotate_composite_membership(variant_grid_columns)

    def _example_row(vgc, members) -> dict:
        """ The client row the example cell draws. A member this build doesn't annotate is left out, so
            the example goes as sparse as the cell does in the grid """
        paths = {m.column.variant_column for m in members} | {vgc.variant_column}
        row = {path: value for path, value in COMPOSITE_EXAMPLE_ROWS[vgc.pk].items() if path in paths}
        return _format_afs(row, AF_UNIT_COLUMNS)

    def _composite_section(vgc, member_rows, row) -> dict:
        members = [member for member, _ in member_rows]
        rich_column = composite_rich_column(vgc, members, column_overrides)
        return {
            "column": vgc,
            "level": vgc.get_annotation_level_display() or "",
            "level_css": (vgc.get_annotation_level_display() or "").lower(),
            "rows": [row for _, row in member_rows],
            "sort_by": members[0].column.label,
            "column_json": json.dumps(rich_column_json(rich_column, AbstractVariantGrid.default_column_width)),
            "row_json": json.dumps(row),
        }

    # (composite, [(member, its row)]) for the cells this build has anything to draw in - a composite
    # whose members are all unannotated here isn't in the grid either, so the page doesn't list it
    drawn_composites = []
    for vgc in variant_grid_columns:
        if not vgc.is_composite:
            continue
        member_rows = [(m, dict(row, role=_member_role(m))) for m in vgc.composite_members.all()
                       if (row := _column_row(m.column)) is not None]
        if member_rows:
            drawn_composites.append((vgc, member_rows))

    # gnomADVariant() and the representative variant cell read the locus columns straight off the row,
    # so every example is drawn on the one fictional variant
    variant_row = next(_example_row(vgc, [m for m, _ in member_rows])
                       for vgc, member_rows in drawn_composites if vgc.pk == "variant")
    composite_sections = [
        _composite_section(vgc, member_rows,
                           {**variant_row, **_example_row(vgc, [m for m, _ in member_rows])})
        for vgc, member_rows in drawn_composites
    ]
    composite_sections.append(_sample_composite_section(_format_afs))
    composite_sections.sort(key=lambda s: (_COMPOSITE_LEVEL_ORDER.index(s["column"].annotation_level),
                                           s["column"].label.lower()))

    for vgc in sorted(variant_grid_columns, key=lambda c: c.pk):
        if vgc.is_composite or vgc.composite_of:
            continue  # members are documented in their composite's section
        if vgc.annotation_level in annotation_run_levels:
            if row := _column_row(vgc):
                columns_by_annotation_level[vgc.get_annotation_level_display()].append(row)
        else:
            variantgrid_columns_by_annotation_level[vgc.annotation_level].append(vgc)

    context = {
        "genome_build": genome_build,
        "composite_sections": composite_sections,
        "example_extra": variant_grid_client_extra(genome_build),
        "variantgrid_columns_by_annotation_level": variantgrid_columns_by_annotation_level,
        "columns_by_annotation_level": columns_by_annotation_level,
    }
    return render(request, "annotation/view_annotation_descriptions.html", context)
