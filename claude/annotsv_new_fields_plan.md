# Plan — Capture more AnnotSV fields (columns_version 4)

Follow-up to the AnnotSV ingestion landed in #1533. Issue #1040 (umbrella SV
annotation wishlist) calls out a number of AnnotSV columns we read from the
TSV but throw away today. This plan extends the AnnotSV ingestion to capture
them. Lands directly into **`columns_version` 4**, which has not been
deployed yet, so no backward-compat / per-version gating is required for
these fields — they are simply present on every v4 SV-pipeline row.

## Currently captured

Source of truth: `FULL_COLUMN_MAP` in
`annotation/vcf_files/bulk_annotsv_tsv_inserter.py` (and the corresponding
fields on `VariantAnnotation` in `annotation/models/models.py:1189-1206`).

| TSV column | Model field |
|---|---|
| `ACMG_class` | `annotsv_acmg_class` |
| `AnnotSV_ranking_score` | `annotsv_acmg_score` |
| `RE_gene` | `annotsv_re_gene` |
| `Repeat_type_left` / `Repeat_type_right` | `annotsv_repeat_type_left/right` |
| `SegDup_left` / `SegDup_right` | `annotsv_segdup_left/right` |
| `ENCODE_blacklist_left` / `ENCODE_blacklist_right` | `annotsv_encode_blacklist_left/right` |
| `ENCODE_blacklist_characteristics_left/right` | `annotsv_encode_blacklist_characteristics_left/right` |
| `B_{gain,loss,ins,inv}_AFmax` | `annotsv_b_{gain,loss,ins,inv}_af_max` |

## New fields to add (this plan)

Targeted at the user wishlist plus a few obvious neighbours that come for
free from the same TSV row.

### Wishlist columns

| AnnotSV TSV column | New model field | Type | Notes |
|---|---|---|---|
| `Frameshift` | `annotsv_frameshift` | `BooleanField(null=True)` | AnnotSV emits `yes`/`no` (or empty); coerce |
| `Exons_spanned` | `annotsv_exons_spanned` | `IntegerField(null=True)` | Count of exons fully spanned |
| `Dist_nearest_SS` | `annotsv_dist_nearest_ss` | `IntegerField(null=True)` | bp to nearest splice site |
| `Nearest_SS_type` | `annotsv_nearest_ss_type` | `TextField(null=True, blank=True)` | `5'` / `3'` |
| `OMIM_inheritance` | `annotsv_omim_inheritance` | `TextField(null=True, blank=True)` | comma/`&`-joined codes (AD/AR/XL/...) |
| `OMIM_morbid` | `annotsv_omim_morbid` | `BooleanField(null=True)` | yes/no |
| `OMIM_phenotype` | `annotsv_omim_phenotype` | `TextField(null=True, blank=True)` | free text, `;`-separated |
| `OMIM_ID` | `annotsv_omim_id` | `TextField(null=True, blank=True)` | one or more MIM numbers |

### ClinVar / pathogenic SV overlaps

AnnotSV does not have a single "ClinVar overlap" column. Pathogenic SV
overlap data lives in the `P_*` family — `P_{event}_{phen,hpo,source,coord}`
for `event` in `{gain, loss, ins, inv}` — 16 columns total. Sources include
ClinVar, dbVar, ClinGen, OMIM-morbid (the `_source` field carries the
source string).

These fields are **reference-only**: rendered on the variant detail page,
used to build URLs to the upstream record, never queried/sorted/grouped
on. The `_phen`/`_hpo` are long free-text and `_source` is `&`-joined,
so none of them are usefully filterable as grid columns.

Capture them as a **single nullable JSONB column** rather than 16 typed
TextFields:

| New model field | Type | Notes |
|---|---|---|
| `annotsv_pathogenic_overlaps` | `JSONField(null=True)` | One nested dict per SV event-type that has any non-empty value |

Wire format (only events with at least one non-empty value are included;
fields within an event with `NA`/empty are dropped):

```json
{
  "gain": {"source": "ClinVar&dbVar", "phen": "...", "hpo": "HP:0001234&HP:0005678", "coord": "chr3:127500000-128800000"},
  "loss": {"source": "ClinGen", "phen": "...", "coord": "chr5:1000-2000"}
}
```

Why JSONB rather than 16 TextFields:

- `VariantAnnotation` is hundreds of millions of rows, ~99% SNVs (NULL
  for all AnnotSV fields). On Postgres, a NULL JSONB costs 1 bit in the
  null bitmap and 0 payload bytes — same as a NULL TextField. 16 nullable
  TextFields cost 16 bits × N rows of bitmap; one JSONB costs 1 bit × N.
  At 200M rows that's ~375MB saved on the bitmap alone.
- Populated rows pay slightly more than equivalent typed columns
  (embedded keys), but only on the ~1% of rows that are SVs. Net win.
- Single `VariantGridColumn` entry with a custom renderer instead of 16.
- If we ever want to filter on `_source` (e.g. "has ClinVar pathogenic
  overlap") we add a derived boolean later — cheaper than carrying 16
  columns now on speculation.

The other AnnotSV fields stay as typed columns: `B_*_AFmax` are
filterable numerics, `OMIM_*` are short and at least `_morbid` /
`_inheritance` are realistic filter targets.

### Rule trail for the existing ACMG score

| AnnotSV TSV column | New model field | Type |
|---|---|---|
| `AnnotSV_ranking_criteria` | `annotsv_acmg_criteria` | `TextField(null=True, blank=True)` |

Useful to expose alongside the score (`annotsv_acmg_score`) as it's the
list of triggered ACMG-SV rules.

## Code changes

### 1. Model — `annotation/models/models.py`

Append to the AnnotSV block (after line 1206). Example for the wishlist
group:

```python
    # AnnotSV (cont.)
    annotsv_acmg_criteria = models.TextField(null=True, blank=True)        # AnnotSV_ranking_criteria
    annotsv_frameshift = models.BooleanField(null=True)                    # Frameshift (yes/no)
    annotsv_exons_spanned = models.IntegerField(null=True, blank=True)     # Exons_spanned
    annotsv_dist_nearest_ss = models.IntegerField(null=True, blank=True)   # Dist_nearest_SS
    annotsv_nearest_ss_type = models.TextField(null=True, blank=True)      # Nearest_SS_type
    annotsv_omim_inheritance = models.TextField(null=True, blank=True)     # OMIM_inheritance
    annotsv_omim_morbid = models.BooleanField(null=True)                   # OMIM_morbid
    annotsv_omim_phenotype = models.TextField(null=True, blank=True)       # OMIM_phenotype
    annotsv_omim_id = models.TextField(null=True, blank=True)              # OMIM_ID

    # AnnotSV pathogenic-SV overlap summary (ClinVar / dbVar / ClinGen / OMIM-morbid).
    # Reference-only: 16 P_{event}_{phen,hpo,source,coord} columns folded into one
    # JSONB blob keyed by event-type. See plan §"ClinVar / pathogenic SV overlaps".
    annotsv_pathogenic_overlaps = models.JSONField(null=True)
```

### 2. Migration — `annotation/migrations/0142_annotsv_more_fields.py`

Use `makemigrations` to generate; the operations should be a flat list of
`AddField` operations matching the model additions above. Dependency:
`("annotation", "0141_variantannotationversion_transcript_resolver")`.

### 3. TSV ingestion — `annotation/vcf_files/bulk_annotsv_tsv_inserter.py`

Extend `FULL_COLUMN_MAP` for the typed scalars, add `BOOL_FIELDS` + boolean
parsing, and assemble the pathogenic-overlap JSON separately in
`_row_to_update`:

```python
FULL_COLUMN_MAP: dict[str, str] = {
    # ... existing entries unchanged ...
    "AnnotSV_ranking_criteria": "annotsv_acmg_criteria",
    "Frameshift": "annotsv_frameshift",
    "Exons_spanned": "annotsv_exons_spanned",
    "Dist_nearest_SS": "annotsv_dist_nearest_ss",
    "Nearest_SS_type": "annotsv_nearest_ss_type",
    "OMIM_inheritance": "annotsv_omim_inheritance",
    "OMIM_morbid": "annotsv_omim_morbid",
    "OMIM_phenotype": "annotsv_omim_phenotype",
    "OMIM_ID": "annotsv_omim_id",
}

INT_FIELDS = {
    "annotsv_acmg_class",
    "annotsv_exons_spanned",
    "annotsv_dist_nearest_ss",
}
BOOL_FIELDS = {
    "annotsv_frameshift",
    "annotsv_omim_morbid",
}
# FLOAT_FIELDS unchanged

_BOOL_TRUE = {"yes", "true", "1"}
_BOOL_FALSE = {"no", "false", "0"}

_P_EVENTS = ("gain", "loss", "ins", "inv")
_P_SUBFIELDS = ("source", "phen", "hpo", "coord")


def _parse_value(field: str, raw: str) -> Optional[Any]:
    if raw is None:
        return None
    raw = raw.strip()
    if raw in EMPTY_VALUES:
        return None
    try:
        if field in INT_FIELDS:
            return int(raw)
        if field in FLOAT_FIELDS:
            return float(raw)
        if field in BOOL_FIELDS:
            low = raw.lower()
            if low in _BOOL_TRUE:
                return True
            if low in _BOOL_FALSE:
                return False
            return None
    except (TypeError, ValueError):
        return None
    return raw


def _extract_pathogenic_overlaps(row: dict[str, str]) -> Optional[dict]:
    """ Collapse AnnotSV's 16 P_{event}_{sub} columns into a nested dict.
        Empty/NA values are dropped; events with no surviving values are
        omitted; returns None if nothing remained. """
    result: dict[str, dict[str, str]] = {}
    for event in _P_EVENTS:
        event_data: dict[str, str] = {}
        for sub in _P_SUBFIELDS:
            raw = row.get(f"P_{event}_{sub}")
            if raw is None:
                continue
            value = raw.strip()
            if value in EMPTY_VALUES:
                continue
            event_data[sub] = value
        if event_data:
            result[event] = event_data
    return result or None
```

Then in `_row_to_update`, after the `FULL_COLUMN_MAP` loop:

```python
def _row_to_update(row: dict[str, str]) -> dict[str, Any]:
    update: dict[str, Any] = {}
    for tsv_col, model_field in FULL_COLUMN_MAP.items():
        if tsv_col not in row:
            continue
        value = _parse_value(model_field, row[tsv_col])
        if value is not None:
            update[model_field] = value
    overlaps = _extract_pathogenic_overlaps(row)
    if overlaps is not None:
        update["annotsv_pathogenic_overlaps"] = overlaps
    return update
```

### 4. Surface as VariantGridColumn entries — new `snpdb` migration

Mirror `snpdb/migrations/0122_new_gnomad_sv_overlap_column_vep_fields2.py`.
File: `snpdb/migrations/0XYZ_new_annotsv_variantgrid_columns.py`.

**Scope is all 24 AnnotSV fields, not just the new ones.** #1533 added the
14 original AnnotSV columns to `VariantAnnotation` but never registered
them as `VariantGridColumn`, so the `description` strings — which are
what `VariantGridColumn.get_column_descriptions()` returns into the
template's `annotation_description` dict — don't exist for them. This
migration backfills both the original 14 AND adds the new 10, with
populated `description` so all 24 rows on the variant detail page can
use `help=annotation_description.<field>` consistently.

```python
from django.db import migrations
from django.db.models import Max, F

from library.django_utils import bulk_insert_class_data


_ANNOTSV_LINK = "<a href='https://lbgi.fr/AnnotSV/' target='_blank'>AnnotSV</a>"

# All AnnotSV fields. The first block is the 14 columns shipped in #1533
# that never got VariantGridColumn entries; the second block is new in
# this change.
_ANNOTSV_COLUMNS = [
    # (grid_column_name, label, description)

    # ---- pre-existing #1533 fields (backfill) ----
    ("annotsv_acmg_class",
     "AnnotSV ACMG class",
     f"{_ANNOTSV_LINK} ACMG-style ranking class for SVs (1=benign, 2=likely benign, 3=VUS, 4=likely pathogenic, 5=pathogenic)."),
    ("annotsv_acmg_score",
     "AnnotSV ACMG score",
     f"{_ANNOTSV_LINK} ACMG-style ranking score: sum of points across triggered rules."),
    ("annotsv_re_gene",
     "Regulatory element gene",
     "Regulatory element gene(s) overlapping the SV, per AnnotSV's regulatory-element bundle (RE_gene)."),
    ("annotsv_repeat_type_left",
     "Repeat type (left)",
     "Repeat type at the SV left breakpoint, per UCSC RepeatMasker (Repeat_type_left)."),
    ("annotsv_repeat_type_right",
     "Repeat type (right)",
     "Repeat type at the SV right breakpoint, per UCSC RepeatMasker (Repeat_type_right)."),
    ("annotsv_segdup_left",
     "SegDup (left)",
     "UCSC segmental duplication feature overlapping the SV left breakpoint (SegDup_left)."),
    ("annotsv_segdup_right",
     "SegDup (right)",
     "UCSC segmental duplication feature overlapping the SV right breakpoint (SegDup_right)."),
    ("annotsv_encode_blacklist_left",
     "ENCODE blacklist (left)",
     "ENCODE blacklist region overlapping the SV left breakpoint."),
    ("annotsv_encode_blacklist_right",
     "ENCODE blacklist (right)",
     "ENCODE blacklist region overlapping the SV right breakpoint."),
    ("annotsv_encode_blacklist_characteristics_left",
     "ENCODE blacklist characteristics (left)",
     "ENCODE blacklist characteristic at the SV left breakpoint (e.g. Low Mappability, High Signal Region)."),
    ("annotsv_encode_blacklist_characteristics_right",
     "ENCODE blacklist characteristics (right)",
     "ENCODE blacklist characteristic at the SV right breakpoint."),
    ("annotsv_b_gain_af_max",
     "AnnotSV benign gain AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV gain (duplication/insertion-gain) sources (B_gain_AFmax)."),
    ("annotsv_b_loss_af_max",
     "AnnotSV benign loss AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV loss (deletion) sources (B_loss_AFmax)."),
    ("annotsv_b_ins_af_max",
     "AnnotSV benign ins AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV insertion sources (B_ins_AFmax)."),
    ("annotsv_b_inv_af_max",
     "AnnotSV benign inv AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV inversion sources (B_inv_AFmax)."),

    # ---- new in this change ----
    ("annotsv_acmg_criteria",
     "AnnotSV ACMG criteria",
     f"Triggered {_ANNOTSV_LINK} ACMG-SV rules (e.g. 1A, 2H) that summed into the ranking score."),
    ("annotsv_frameshift",
     "AnnotSV frameshift",
     "Whether the SV introduces a frameshift in any overlapped transcript."),
    ("annotsv_exons_spanned",
     "AnnotSV exons spanned",
     "Number of exons fully spanned by the SV."),
    ("annotsv_dist_nearest_ss",
     "Distance to nearest splice site",
     "Distance in bp from the SV breakpoint to the nearest splice site (Dist_nearest_SS)."),
    ("annotsv_nearest_ss_type",
     "Nearest splice site type",
     "Type of the nearest splice site (5' donor or 3' acceptor)."),
    ("annotsv_omim_inheritance",
     "OMIM inheritance",
     "OMIM inheritance pattern(s) for genes overlapped by the SV (AD, AR, XL, etc.)."),
    ("annotsv_omim_morbid",
     "OMIM morbid",
     "Whether any gene overlapped by the SV is listed in the OMIM morbid map."),
    ("annotsv_omim_phenotype",
     "OMIM phenotype",
     "OMIM phenotype text for genes overlapped by the SV."),
    ("annotsv_omim_id",
     "OMIM ID",
     "OMIM MIM number(s) for genes overlapped by the SV."),
    ("annotsv_pathogenic_overlaps",
     "Pathogenic SV overlaps",
     f"Per-event-type (gain/loss/ins/inv) summary of pathogenic SV overlaps from {_ANNOTSV_LINK}'s pathogenic-SV bundle: source (ClinVar / dbVar / ClinGen / OMIM-morbid), reported phenotype, HPO terms, and reference SV coordinates."),
]


def _add_annotsv_columns(apps, _schema_editor):
    rows = [
        {
            "grid_column_name": name,
            "variant_column": f"variantannotation__{name}",
            "annotation_level": "V",
            "width": None,
            "label": label,
            "description": description,
            "model_field": True,
            "queryset_field": True,
        }
        for name, label, description in _ANNOTSV_COLUMNS
    ]
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", rows)])

    # Append to "All columns". sort_order_max is computed across whatever
    # SV-adjacent columns already exist in the collection; this is the
    # first time any annotsv_ column gets registered, so we anchor off
    # the SV overlap block instead.
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")
    all_columns = CustomColumnsCollection.objects.get(name="All columns")
    sv_qs = all_columns.customcolumn_set.filter(
        column__grid_column_name__startswith="gnomad_sv_overlap",
        column__annotation_level="V",
    )
    sort_order_max = sv_qs.aggregate(Max("sort_order"))["sort_order__max"] or 0
    all_columns.customcolumn_set.filter(sort_order__gt=sort_order_max).update(
        sort_order=F("sort_order") + len(_ANNOTSV_COLUMNS)
    )
    for i, (name, _label, _desc) in enumerate(_ANNOTSV_COLUMNS, start=1):
        CustomColumn.objects.create(
            custom_columns_collection=all_columns,
            sort_order=sort_order_max + i,
            column_id=name,
        )


def _remove_annotsv_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(
        grid_column_name__in=[name for name, _label, _desc in _ANNOTSV_COLUMNS]
    ).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "<previous head>"),  # fill in with `showmigrations snpdb | tail -1`
    ]

    operations = [
        migrations.RunPython(_add_annotsv_columns, reverse_code=_remove_annotsv_columns),
    ]
```

### 5. Variant detail page — `variantopedia/templates/variantopedia/variant_details.html`

The existing AnnotSV block lives at lines ~805-841. There are two
sub-areas:

- **Always-visible header** (lines 805-808): ACMG class + score.
- **Collapsed `#annotsv-context` panel** (lines 813-841): regulatory /
  repeat / SegDup / ENCODE / benign-AF rows, each gated by `{% if %}` so
  it only renders when populated.

**While we're in here, also add `help=annotation_description.<field>` to
every existing AnnotSV `{% labelled %}` row** (the 14 rows from #1533).
The snpdb migration in §4 backfills `description` for those columns,
so the help-icon hover-cards will populate automatically. This is part
of the same "fill in `annotation_description` for all AnnotSV fields"
scope.

Extend each:

**Header — append the ACMG criteria string under the score** (visible
alongside the class/score it explains; back-fill `help=` on the two
pre-existing rows too):

```django
{% if variant_annotation.annotsv_acmg_class is not None or variant_annotation.annotsv_acmg_score is not None %}
    {% labelled id='annotsv_acmg_class' label='AnnotSV ACMG class' help=annotation_description.annotsv_acmg_class hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_acmg_class|default_if_none:"-" }}{% endlabelled %}
    {% labelled id='annotsv_acmg_score' label='AnnotSV ACMG score' help=annotation_description.annotsv_acmg_score hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_acmg_score|default_if_none:"-" }}{% endlabelled %}
    {% if variant_annotation.annotsv_acmg_criteria %}
        {% labelled id='annotsv_acmg_criteria' label='AnnotSV ACMG criteria' help=annotation_description.annotsv_acmg_criteria hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_acmg_criteria }}{% endlabelled %}
    {% endif %}
{% endif %}
```

**Inside `#annotsv-context`, after the existing benign-AF rows
(line 840)** — add a "Gene impact" group, an "OMIM" group, and a
"Pathogenic SV overlaps" group. Each row is `{% if %}`-gated so SNV /
SV-without-AnnotSV rows render nothing extra. Every row passes
`help=annotation_description.<field>` so the existing description
hover-card mechanism picks up the strings registered in §4.

```django
{# Gene impact — Frameshift / Exons spanned / Nearest splice site #}
{% if variant_annotation.annotsv_frameshift is not None %}
    {% labelled id='annotsv_frameshift' label='AnnotSV frameshift' help=annotation_description.annotsv_frameshift hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_frameshift|yesno:"yes,no" }}{% endlabelled %}
{% endif %}
{% if variant_annotation.annotsv_exons_spanned is not None %}
    {% labelled id='annotsv_exons_spanned' label='AnnotSV exons spanned' help=annotation_description.annotsv_exons_spanned hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_exons_spanned }}{% endlabelled %}
{% endif %}
{% if variant_annotation.annotsv_dist_nearest_ss is not None or variant_annotation.annotsv_nearest_ss_type %}
    {% labelled id='annotsv_nearest_ss' label='AnnotSV nearest splice site' help=annotation_description.annotsv_dist_nearest_ss hint="tiny" value_css="text-monospace" %}
        {{ variant_annotation.annotsv_nearest_ss_type|default_if_none:"-" }} ({{ variant_annotation.annotsv_dist_nearest_ss|default_if_none:"-" }} bp)
    {% endlabelled %}
{% endif %}

{# OMIM #}
{% if variant_annotation.annotsv_omim_id or variant_annotation.annotsv_omim_inheritance or variant_annotation.annotsv_omim_phenotype or variant_annotation.annotsv_omim_morbid is not None %}
    {% if variant_annotation.annotsv_omim_id %}
        {% labelled id='annotsv_omim_id' label='OMIM ID' help=annotation_description.annotsv_omim_id hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_omim_id }}{% endlabelled %}
    {% endif %}
    {% if variant_annotation.annotsv_omim_phenotype %}
        {% labelled id='annotsv_omim_phenotype' label='OMIM phenotype' help=annotation_description.annotsv_omim_phenotype hint="tiny" %}{{ variant_annotation.annotsv_omim_phenotype }}{% endlabelled %}
    {% endif %}
    {% if variant_annotation.annotsv_omim_inheritance %}
        {% labelled id='annotsv_omim_inheritance' label='OMIM inheritance' help=annotation_description.annotsv_omim_inheritance hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_omim_inheritance }}{% endlabelled %}
    {% endif %}
    {% if variant_annotation.annotsv_omim_morbid is not None %}
        {% labelled id='annotsv_omim_morbid' label='OMIM morbid' help=annotation_description.annotsv_omim_morbid hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_omim_morbid|yesno:"yes,no" }}{% endlabelled %}
    {% endif %}
{% endif %}

{# Pathogenic SV overlaps — JSONB blob.
   Inline summary (event + source); click "Details" to expand the full
   phenotype / HPO / coords behind a per-event Bootstrap collapse. #}
{% if variant_annotation.annotsv_pathogenic_overlaps %}
    {% labelled label="Pathogenic SV overlaps" help=annotation_description.annotsv_pathogenic_overlaps hint="tiny" %}
        <a class='toggle-link' data-toggle="collapse" href='#annotsv-p-overlaps'>Show pathogenic SV overlaps</a>
    {% endlabelled %}
    <div id="annotsv-p-overlaps" class="collapse">
        {% for event, fields in variant_annotation.annotsv_pathogenic_overlaps.items %}
            {% with collapse_id="annotsv-p-overlap-"|add:event %}
                {% labelled label="Pathogenic"|add:" "|add:event hint="tiny" value_css="text-monospace" %}
                    {% if fields.source %}<strong>{{ fields.source }}</strong>{% endif %}
                    {% if fields.coord %} — {{ fields.coord }}{% endif %}
                    <a class='toggle-link ml-2' data-toggle="collapse" href='#{{ collapse_id }}'>Details</a>
                    <div id="{{ collapse_id }}" class="collapse mt-1">
                        <dl class="row mb-0">
                            {% if fields.phen %}
                                <dt class="col-sm-3">Phenotype</dt>
                                <dd class="col-sm-9">{{ fields.phen }}</dd>
                            {% endif %}
                            {% if fields.hpo %}
                                <dt class="col-sm-3">HPO</dt>
                                <dd class="col-sm-9">{{ fields.hpo }}</dd>
                            {% endif %}
                            {% if fields.coord %}
                                <dt class="col-sm-3">Coords</dt>
                                <dd class="col-sm-9">{{ fields.coord }}</dd>
                            {% endif %}
                            {% if fields.source %}
                                <dt class="col-sm-3">Source</dt>
                                <dd class="col-sm-9">{{ fields.source }}</dd>
                            {% endif %}
                        </dl>
                    </div>
                {% endlabelled %}
            {% endwith %}
        {% endfor %}
    </div>
{% endif %}
```

**Back-fill `help=annotation_description.<field>` on the 11 pre-existing
`#annotsv-context` rows from #1533** (lines 814-840). Pattern, applied
to each existing row:

```django
{# before #}
{% labelled id='annotsv_re_gene' label='Regulatory element gene' hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_re_gene }}{% endlabelled %}
{# after #}
{% labelled id='annotsv_re_gene' label='Regulatory element gene' help=annotation_description.annotsv_re_gene hint="tiny" value_css="text-monospace" %}{{ variant_annotation.annotsv_re_gene }}{% endlabelled %}
```

Apply the same `help=annotation_description.<field>` insertion to:
`annotsv_repeat_left`, `annotsv_repeat_right`, `annotsv_segdup_left`,
`annotsv_segdup_right`, `annotsv_encode_blacklist_left`,
`annotsv_encode_blacklist_right`, `annotsv_b_gain_af_max`,
`annotsv_b_loss_af_max`, `annotsv_b_ins_af_max`, `annotsv_b_inv_af_max`.

Notes for the implementer:

- Bootstrap 4 — per `CLAUDE.md`, use `data-toggle` / `href` (not the BS5
  `data-bs-toggle` / `data-bs-target`).
- Two-level disclosure: outer collapse (`#annotsv-p-overlaps`) wraps the
  whole list; per-event "Details" toggles expand the long-form
  phen/hpo/coords for one event-type at a time. Coords appear inline
  *and* in the details panel deliberately — the inline copy is for quick
  scanning, the details panel is for copy-paste / URL-building.
- Collapse IDs are namespaced per event (`annotsv-p-overlap-gain`, etc.)
  so multiple panels can be open at once and don't clash with the
  enclosing `#annotsv-context` collapse.
- All new rows live inside the existing `#annotsv-context` collapse so
  they don't add to the always-visible variant detail height for SVs
  with sparse AnnotSV output.
- The same nested-collapse pattern can be applied to any other long-form
  AnnotSV field if it gets too noisy — e.g. wrap a long
  `annotsv_acmg_criteria` string the same way if the rule trail is
  unwieldy. Default in this plan keeps `annotsv_acmg_criteria` inline
  since it's typically short.
- Verify in a browser on an SV variant once data is populated — start
  the dev server (`python3 manage.py runserver`) and load a variant
  detail page for a known SV. Test both the outer toggle and at least
  one per-event "Details" toggle. Confirm SNV variant detail pages are
  unchanged.

### 6. Bump `columns_version` definitions where relevant

Most of these fields don't need entries in `vep_columns.py` — they're
populated by the AnnotSV stage, not by VEP. The columns_version 4 bump
already gates the SV pipeline as a whole; nothing more is required there.

### 7. Tests — `annotation/tests/test_annotsv.py`

Add coverage for:

- `test_full_column_map_parses_new_fields` — feed a synthetic row through
  `_row_to_update` covering each new typed column (incl. boolean coercion
  of `Frameshift`/`OMIM_morbid` and integer coercion of `Exons_spanned`,
  `Dist_nearest_SS`).
- `test_bool_field_na` — `NA` → `None`, `yes` → `True`, `no` → `False`,
  bogus value → `None`.
- `test_pathogenic_overlaps_assembly` — synthetic row with a mix of
  populated and `NA` `P_*_*` fields; assert the resulting JSON drops
  empties, omits events with nothing populated, and returns `None` when
  all 16 cells are `NA`.
- Extend the existing fixture TSV
  `annotation/tests/test_data/annotsv/test_grch37_sv.annotated.tsv`
  (single file — there is no grch38 AnnotSV TSV in the repo today; the
  end-to-end test in `test_annotsv.py` only exercises grch37) with the
  new column names in the header (incl. all 16 `P_*_*` headers) and at
  least one row carrying real values for both the typed scalars and a
  couple of populated P_* events, so the end-to-end import test asserts
  the typed fields land and `annotsv_pathogenic_overlaps` is the
  expected nested dict. The end-to-end test pairs this with a
  VEP-annotated VCF (`TEST_SV_VCF_GRCH37 = ...test_columns_version3_grch37_sv.vep_annotated.vcf`);
  swap to the new `test_columns_version4_grch37_sv.vep_annotated.vcf`
  fixture when v4 lands (currently untracked in `git status`).

## Order of work

1. Add fields to `VariantAnnotation` model.
2. Generate + commit migration `0142_annotsv_more_fields`.
3. Extend `FULL_COLUMN_MAP` + parser in `bulk_annotsv_tsv_inserter.py`.
4. Extend test TSV (`test_grch37_sv.annotated.tsv`) and `test_annotsv.py`.
5. New `snpdb` `VariantGridColumn` migration registering **all 24
   AnnotSV fields** (the 14 from #1533 plus the 10 new ones) with
   populated `description` strings. This is what powers
   `annotation_description.<field>` in the template.
6. Add the new fields to `variantopedia/templates/variantopedia/variant_details.html`
   inside the existing `#annotsv-context` collapse, with a nested
   Bootstrap collapse per pathogenic-overlap event-type. **Also
   back-fill `help=annotation_description.<field>` on the 11
   pre-existing `#annotsv-context` rows from #1533.**
7. Run `python3 manage.py test --keepdb annotation.tests.test_annotsv`.
8. Smoke-test the variant detail page on an SV (dev server) and on an SNV
   (confirm no regression). Test outer + per-event Bootstrap toggles.
9. Backfill recipe in PR description: `python3 manage.py annotsv_run
   --annotation-run <pk> --skip-run` re-imports existing on-disk TSVs
   without re-running AnnotSV.

## Open questions

- **Boolean vs text for `Frameshift` / `OMIM_morbid`?** Plan stores them
  as nullable booleans. AnnotSV occasionally emits unexpected tokens; the
  parser falls back to `None` in that case. If we'd rather preserve the
  raw string, switch them to TextField.
- **Split-line ingestion (#1533 follow-up)** is still deferred. Several of
  these fields (`Exons_spanned`, `Frameshift`, `Nearest_SS_type`,
  `Dist_nearest_SS`) are *more* meaningful on split-lines (per-gene) than
  full-lines. For this plan we accept the full-line value — it's the
  worst-case across overlapped genes — and revisit when split-line
  persistence lands.

## Cross-checks against the prior #1533 commit (`c34f1b555`)

Re-read the original AnnotSV change diff to make sure nothing it set up
also needs touching here:

- **`VariantAnnotationVersion.annotsv_code` / `annotsv_bundle`** — version
  pins added in #1533. **Untouched** here. New *fields* don't change the
  AnnotSV binary or bundle version, so no version bump and no automatic
  re-annotation trigger is needed. See "Backfill" below for how existing
  rows get the new columns populated.
- **`AnnotationRun.annotsv_tsv_filename` / `annotsv_error` /
  `annotsv_imported`** — pipeline plumbing added in #1533. **Untouched**.
- **`annotation/tasks/annotate_variants.py`** (runs AnnotSV after VEP) —
  no changes needed; it writes the TSV path onto the run, the inserter
  reads it.
- **`annotation/vcf_files/import_vcf_annotations.py`** (calls
  `import_annotsv_tsv` best-effort, logs errors) — no changes needed.
- **`annotation/management/commands/annotsv_run.py`** — backfill command
  (re-runs AnnotSV and/or re-imports the TSV for a single AnnotationRun).
  *Used* by this change for backfill, not modified.
- **`variantgrid/settings/components/annotation_settings.py`** —
  `ANNOTATION_ANNOTSV_*` settings. **Untouched**.
- **`annotation/annotsv_annotation.py`** — command-line builder + version
  check. **Untouched**.
- **`ColumnVCFInfo`** — #1533 did not register any AnnotSV fields here.
  This plan also leaves it out — none of the 24 AnnotSV fields are
  currently expected to round-trip through VCF export. Add later if a
  consumer needs them.
- **`VariantGridColumn` registration covers all 24 AnnotSV fields** —
  #1533 added the 14 original AnnotSV columns to `VariantAnnotation`
  but never registered them as `VariantGridColumn`. This change fixes
  that gap (see §4) so the whole AnnotSV column family becomes
  grid-pickable and `annotation_description.<field>` resolves for every
  row.
- **`annotation_description.<field>` template help** — back-filled on
  the 11 pre-existing `#annotsv-context` rows from #1533 (see §5)
  alongside the new rows. No AnnotSV row on the variant detail page is
  left without a help hover-card after this change.

### Backfill

`VariantAnnotation` rows already annotated under #1533 will have NULL
for all the new columns until re-imported. The `annotsv_run` management
command from #1533 supports this:

```bash
# Re-import existing TSVs (AnnotSV not re-run) for one annotation run.
python3 manage.py annotsv_run --annotation-run <pk> --skip-run
```

`--skip-run` skips invoking AnnotSV and just re-runs `import_annotsv_tsv`
against the on-disk TSV, which now reads the new columns. For bulk
backfill, iterate over `AnnotationRun.objects.filter(pipeline_type='C',
annotsv_imported=True)` PKs.

Worth adding a one-liner to the PR description so deployers know the
backfill recipe.

## Resolved

- **`P_*_coord` storage** — folded into `annotsv_pathogenic_overlaps`
  JSONB along with `_phen` / `_hpo` / `_source`. All 16 fields are
  reference-only (variant detail page rendering, URL building); none are
  usefully filterable as grid columns. Single nullable JSONB on a 200M-row
  table is cheaper than 16 nullable TextFields (1 null-bitmap bit vs 16
  per row, ~375MB saved). If filtering on `_source` is wanted later, add
  a derived boolean.
- **No `SVAnnotation` split** — keeping AnnotSV columns on
  `VariantAnnotation`. NULL columns are effectively free in Postgres at
  any row count; the join cost + parallel versioning lifecycle of a split
  table outweigh the storage saving. Revisit only if SV-specific
  annotation grows substantially (CADD-SV / SVAFotate / DeepSVP from
  #1040).
