# Refactor `ColumnVEPField` → Python registry

Issue: [#1213](https://github.com/SACGF/variantgrid/issues/1213) — *"Refactor ColumnVEPField to be simpler data structures"*

This plan is for review by Dave + codex before implementation.

> **TL;DR** — `ColumnVEPField` is static config, not runtime data. Replace the model + table with a Python registry (`annotation/vep_columns.py`) of frozen dataclasses, then delete the model. Drops 237 DB rows + 10+ historical data migrations of churn into one editable file with enums, IDE support, and import-time validation.

---

## 1. Why this isn't database data

Symptoms from issue #1213:

- "Combinatorial explosion" of build duplicates (237 rows; ~64 are pure 37/38/T2T copies of the same logical mapping).
- `Mastermind_counts` has 6 rows (3 destinations × 2 builds) for what is conceptually one mapping.
- `gnomad_sv_overlap_copied_to_gnomad_af` is a fake row whose `column` is meaningless — it exists only because `column TextField unique` blocks two rows from sharing a destination.

Root cause: this is a **registry of static config**, not data the app reads or writes at runtime.

Evidence:

- Every row was inserted by a migration. There's no admin form, no view that creates one, no API.
- Non-migration call sites read the model in 6 places (`vep_annotation.py:50, 230, 357`; `bulk_vep_vcf_annotation_inserter.py:132, 268, 850`; `views.py:466`; `models_columns.py:50`; `custom_columns.py:18`). All read; none write.
- 10+ historical migrations exist *purely* to add/edit rows: `0024_new_vep_annotation_gnomad3`, `0027_new_column_spliaceai_gene_vep_config`, `0053_new_vep_columns`, `0082_new_vep_110_columns_v3`, `0094_gnomad_sv_vep_fields`, `0097_modify_gnomad_sv_overlap_column_vep_fields`, `0115_one_off_vep_column_fields_vep112`, `0116_one_off_update_vep_112_conservation_sv`, `0121_vep_columns_37_38_only`, `0122_vep_gnomad_t2t`. Adding T2T or VEP 113 fields will keep producing more.
- Only one fact wants to be in the DB: "for this `VariantGridColumn`, is there an applicable VEP mapping under the current `columns_version`?" That's a tiny set lookup against the registry — no FK needed.

Putting this in Python gets us:

| | DB model (today) | Python registry |
|---|---|---|
| Add a new VEP field | new migration `_one_off_*` | edit one file |
| Editor support / typos | string FK ids; caught at runtime | enums + dataclass; caught by linter/IDE |
| Versioning | synthesised `min_/max_*_version` joins | `git log annotation/vep_columns.py` |
| Per-build duplication | one row per build | `genome_builds=frozenset({...})` |
| Multi-destination mapping | extra rows w/ fake `column` | `variant_grid_columns=("a","b","c")` |
| Validation | none beyond DB constraints | import-time assertions + unit test |

---

## 2. Proposed shape

### 2.1 `annotation/vep_columns.py`

```python
from dataclasses import dataclass
from typing import Optional

from annotation.models.models_enums import (
    ColumnAnnotationCategory,
    VariantAnnotationPipelineType,
    VEPCustom,
    VEPPlugin,
)


@dataclass(frozen=True)
class VEPColumnDef:
    """ Maps a VEP CSQ / plugin / custom output field → one or more VariantGridColumn destinations.
        Empty `genome_builds` / `pipeline_types` mean "applies to all". """

    source_field: Optional[str]
    variant_grid_columns: tuple[str, ...]
    category: ColumnAnnotationCategory
    vep_plugin: Optional[VEPPlugin] = None
    vep_custom: Optional[VEPCustom] = None
    source_field_has_custom_prefix: bool = False
    genome_builds: frozenset[str] = frozenset()
    pipeline_types: frozenset[VariantAnnotationPipelineType] = frozenset()
    min_columns_version: Optional[int] = None
    max_columns_version: Optional[int] = None
    min_vep_version: Optional[int] = None
    max_vep_version: Optional[int] = None
    summary_stats: Optional[str] = None
    source_field_processing_description: Optional[str] = None

    # ---- behavioural helpers (replace ColumnVEPField properties / static methods) ----

    @property
    def vep_info_field(self) -> Optional[str]:
        """ Replicates ColumnVEPField.vep_info_field. """
        vif = self.source_field
        if self.vep_custom and self.source_field_has_custom_prefix:
            prefix = self.vep_custom.label
            vif = f"{prefix}_{vif}" if vif else prefix
        return vif

    @property
    def columns_version_description(self) -> str:
        limits = []
        if self.min_columns_version:
            limits.append(f"column version >= {self.min_columns_version}")
        if self.max_columns_version:
            limits.append(f"column version <= {self.max_columns_version}")
        return " and ".join(limits)

    def applies_to(
        self,
        *,
        genome_build_name: Optional[str] = None,
        pipeline_type: Optional[VariantAnnotationPipelineType] = None,
        columns_version: Optional[int] = None,
        vep_version: Optional[int] = None,
    ) -> bool:
        if genome_build_name is not None and self.genome_builds and genome_build_name not in self.genome_builds:
            return False
        if pipeline_type is not None and self.pipeline_types and pipeline_type not in self.pipeline_types:
            return False
        if columns_version is not None:
            if self.min_columns_version is not None and columns_version < self.min_columns_version:
                return False
            if self.max_columns_version is not None and columns_version > self.max_columns_version:
                return False
        if vep_version is not None:
            if self.min_vep_version is not None and vep_version < self.min_vep_version:
                return False
            if self.max_vep_version is not None and vep_version > self.max_vep_version:
                return False
        return True


# ---- the registry ----

VEP_COLUMNS: tuple[VEPColumnDef, ...] = (
    # --- dbNSFP plugin (37/38 only) ---
    VEPColumnDef(
        source_field="Aloft_Confidence",
        variant_grid_columns=("aloft_high_confidence",),
        category=ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
        vep_plugin=VEPPlugin.DBNSFP,
        genome_builds=frozenset({"GRCh37", "GRCh38"}),
        min_columns_version=2,
        source_field_processing_description="Most damaging transcript prediction chosen, and Ensembl transcript stored.",
    ),
    # ... (one entry per logical mapping; ~140 entries total, down from 237 rows)

    # --- Mastermind plugin (one source field, three destinations) ---
    VEPColumnDef(
        source_field="Mastermind_counts",
        variant_grid_columns=(
            "mastermind_count_1_cdna",
            "mastermind_count_2_cdna_prot",
            "mastermind_count_3_aa_change",
        ),
        category=ColumnAnnotationCategory.LITERATURE,
        vep_plugin=VEPPlugin.MASTERMIND,
        genome_builds=frozenset({"GRCh37", "GRCh38"}),
        pipeline_types=frozenset({VariantAnnotationPipelineType.STANDARD}),
    ),

    # --- gnomAD SV overlap (the "copied_to_gnomad_af" anomaly, fixed natively) ---
    VEPColumnDef(
        source_field="SV_overlap_AF",
        variant_grid_columns=("gnomad_sv_overlap_af", "gnomad_af"),
        category=ColumnAnnotationCategory.FREQUENCY_DATA,
        vep_custom=VEPCustom.GNOMAD_SV,
        source_field_has_custom_prefix=True,
        pipeline_types=frozenset({VariantAnnotationPipelineType.STRUCTURAL_VARIANT}),
    ),

    # --- Existing variation produces both COSMIC and dbSNP ids ---
    VEPColumnDef(
        source_field="Existing_variation",
        variant_grid_columns=("cosmic_id", "dbsnp_rs_id"),
        category=ColumnAnnotationCategory.EXTERNAL_ID,
    ),
    # ...
)


# ---- query helpers (replace ColumnVEPField.filter / get / get_source_fields) ----

def filter_for(
    *,
    genome_build_name: Optional[str] = None,
    pipeline_type: Optional[VariantAnnotationPipelineType] = None,
    columns_version: Optional[int] = None,
    vep_version: Optional[int] = None,
    vep_plugin: Optional[VEPPlugin] = None,
    vep_custom: Optional[VEPCustom] = None,
) -> tuple[VEPColumnDef, ...]:
    return tuple(
        c for c in VEP_COLUMNS
        if c.applies_to(
            genome_build_name=genome_build_name,
            pipeline_type=pipeline_type,
            columns_version=columns_version,
            vep_version=vep_version,
        )
        and (vep_plugin is None or c.vep_plugin == vep_plugin)
        and (vep_custom is None or c.vep_custom == vep_custom)
    )


def source_fields_for(**kwargs) -> list[str]:
    """ Distinct, sorted source_fields for the given filter. """
    return sorted({c.source_field for c in filter_for(**kwargs) if c.source_field})


def for_variant_grid_column(vgc_id: str) -> tuple[VEPColumnDef, ...]:
    return tuple(c for c in VEP_COLUMNS if vgc_id in c.variant_grid_columns)


def all_variant_grid_column_ids() -> frozenset[str]:
    return frozenset(vgc for c in VEP_COLUMNS for vgc in c.variant_grid_columns)
```

### 2.2 Import-time validation

In `annotation/apps.py:AnnotationConfig.ready()`:

```python
def ready(self):
    if not settings.UNIT_TEST:
        return  # validation runs as part of system check below
    self._validate_vep_columns()

@staticmethod
def _validate_vep_columns():
    from annotation.vep_columns import VEP_COLUMNS, all_variant_grid_column_ids
    from snpdb.models import VariantGridColumn
    known = set(VariantGridColumn.objects.values_list("pk", flat=True))
    referenced = all_variant_grid_column_ids()
    missing = referenced - known
    if missing:
        raise ImproperlyConfigured(f"vep_columns.py references unknown VariantGridColumn: {sorted(missing)}")
    # Uniqueness: (source_field, plugin, custom, version-bounds) must not collide
    seen = set()
    for c in VEP_COLUMNS:
        k = (c.source_field, c.vep_plugin, c.vep_custom,
             c.min_columns_version, c.max_columns_version,
             c.min_vep_version, c.max_vep_version)
        if k in seen:
            raise ImproperlyConfigured(f"Duplicate VEPColumnDef key: {k}")
        seen.add(k)
```

Plus a unit test (`annotation/tests/test_vep_columns.py`) that asserts the same on every test run.

---

## 3. Call-site rewrites

### 3.1 `annotation/vep_annotation.py:46-53`

Before:
```python
def _get_dbnsfp_plugin_command(genome_build: GenomeBuild, vc: VEPConfig):
    dbnsfp_data_path = vc["dbnsfp"]
    q = ColumnVEPField.get_columns_version_q(vc.columns_version)
    fields = ColumnVEPField.get_source_fields(genome_build, q, vep_plugin=VEPPlugin.DBNSFP)
    joined_columns = ",".join(fields)
    return f"dbNSFP,{dbnsfp_data_path},{joined_columns}"
```

After:
```python
def _get_dbnsfp_plugin_command(genome_build: GenomeBuild, vc: VEPConfig):
    dbnsfp_data_path = vc["dbnsfp"]
    fields = vep_columns.source_fields_for(
        genome_build_name=genome_build.name,
        columns_version=vc.columns_version,
        vep_plugin=VEPPlugin.DBNSFP,
    )
    return f"dbNSFP,{dbnsfp_data_path},{','.join(fields)}"
```

### 3.2 `annotation/vep_annotation.py:228-241`

Before:
```python
for vep_custom, prefix in dict(VEPCustom.choices).items():
    try:
        q = ColumnVEPField.get_q(genome_build, vc.vep_version, vc.columns_version, pipeline_type)
        if cvf_list := list(ColumnVEPField.get(genome_build, q, vep_custom=vep_custom)):
            ...
```

After:
```python
for vep_custom in VEPCustom:
    try:
        cvf_list = vep_columns.filter_for(
            genome_build_name=genome_build.name,
            pipeline_type=pipeline_type,
            columns_version=vc.columns_version,
            vep_version=vc.vep_version,
            vep_custom=vep_custom,
        )
        if cvf_list:
            prefix = vep_custom.label
            prefix_lc = prefix.lower()
            if cfg := vc[prefix_lc]:
                cmd.extend(_get_custom_params_list(cvf_list, prefix, cfg))
            ...
```

### 3.3 `annotation/vep_annotation.py:355-361`

Before:
```python
q_cvf = ColumnVEPField.get_columns_version_q(vep_config.columns_version)
if cvf := ColumnVEPField.objects.filter(q_cvf, variant_grid_column='gnomad_af', genome_build=genome_build).first():
    gnomad_filename = vep_config[cvf.get_vep_custom_display().lower()]
```

After:
```python
candidates = [c for c in vep_columns.for_variant_grid_column("gnomad_af")
              if c.applies_to(genome_build_name=genome_build.name,
                              columns_version=vep_config.columns_version)]
if candidates:
    cvf = candidates[0]
    gnomad_filename = vep_config[cvf.vep_custom.label.lower()]
```

### 3.4 `annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:132-145, 240-270`

Before:
```python
cvf_qs = ColumnVEPField.filter(self.genome_build,
                               self.vep_config.vep_version,
                               self.vep_config.columns_version,
                               self.annotation_run.pipeline_type)
self._setup_vep_fields_and_db_columns(validate_columns, cvf_qs)
...
for cvf in cvf_qs.order_by("source_field"):
    ...
    self.source_field_to_columns[cvf.vep_info_field].add(cvf.variant_grid_column_id)
...
other_cvf_qs = ColumnVEPField.objects.all().difference(cvf_qs)
vep_fields_not_this_version = set(other_cvf_qs.values_list("variant_grid_column_id", flat=True))
ignore_columns.update(vep_fields_not_this_version)
```

After:
```python
cvf_list = vep_columns.filter_for(
    genome_build_name=self.genome_build.name,
    pipeline_type=self.annotation_run.pipeline_type,
    columns_version=self.vep_config.columns_version,
    vep_version=self.vep_config.vep_version,
)
self._setup_vep_fields_and_db_columns(validate_columns, cvf_list)
...
for cvf in sorted(cvf_list, key=lambda c: c.source_field or ""):
    ...
    for vgc_id in cvf.variant_grid_columns:
        self.source_field_to_columns[cvf.vep_info_field].add(vgc_id)
...
in_scope = {vgc for c in cvf_list for vgc in c.variant_grid_columns}
all_known = vep_columns.all_variant_grid_column_ids()
ignore_columns.update(all_known - in_scope)
```

`SVOverlapProcessor.__init__` (line 849) similarly rewrites `cvf_qs.values_list("variant_grid_column_id", ...)` → flatten `c.variant_grid_columns` from a filtered list.

### 3.5 `snpdb/grid_columns/custom_columns.py:17-19`

Before:
```python
q_cvf = ColumnVEPField.get_columns_version_q(annotation_version.variant_annotation_version.columns_version)
cvf_qs = ColumnVEPField.objects.filter(q_cvf)
q_columns_this_version = Q(column__columnvepfield__isnull=True) | Q(column__columnvepfield__in=cvf_qs)
```

After:
```python
columns_version = annotation_version.variant_annotation_version.columns_version
in_version_vgcs = {
    vgc for c in vep_columns.filter_for(columns_version=columns_version)
    for vgc in c.variant_grid_columns
}
ever_referenced = vep_columns.all_variant_grid_column_ids()
# columns that are either not VEP-driven at all, OR are VEP-driven and applicable to this version
q_columns_this_version = ~Q(column__in=ever_referenced) | Q(column__in=in_version_vgcs)
```

### 3.6 `snpdb/models/models_columns.py:48-52`

Before:
```python
@cached_property
def columns_version_description(self) -> str:
    q = Q(min_columns_version__isnull=False) | Q(max_columns_version__isnull=False)
    if cvf := self.columnvepfield_set.filter(q).first():
        return cvf.columns_version_description
    return ""
```

After:
```python
@cached_property
def columns_version_description(self) -> str:
    from annotation.vep_columns import for_variant_grid_column
    for c in for_variant_grid_column(self.pk):
        if c.min_columns_version is not None or c.max_columns_version is not None:
            return c.columns_version_description
    return ""
```

### 3.7 `annotation/views.py:460-480` (descriptions page)

Before:
```python
vep_qs = ColumnVEPField.filter_for_build(genome_build)
for vgc in VariantGridColumn.objects.all().order_by("grid_column_name"):
    if vgc.annotation_level in vep_annotation_levels:
        if vep := vep_qs.filter(variant_grid_column=vgc).first():
            columns_and_vep_by_annotation_level[vgc.get_annotation_level_display()][vgc] = vep
```

After:
```python
def _first_for_build(vgc_id, build_name):
    for c in vep_columns.for_variant_grid_column(vgc_id):
        if c.applies_to(genome_build_name=build_name):
            return c
    return None

for vgc in VariantGridColumn.objects.all().order_by("grid_column_name"):
    if vgc.annotation_level in vep_annotation_levels:
        if vep := _first_for_build(vgc.pk, genome_build.name):
            columns_and_vep_by_annotation_level[vgc.get_annotation_level_display()][vgc] = vep
```

### 3.8 Template — `annotation/templates/annotation/view_annotation_descriptions.html`

The template reads `columnvepfield.category` / `get_category_display` / `source_field` / `source_field_processing_description` / `vep_plugin` / `get_vep_plugin_display` / `vep_custom` / `columns_version_description`. The dataclass has all of these as fields/properties — but `get_category_display` and `get_vep_plugin_display` are model-method names that don't exist on a dataclass.

Options:
- **(a)** Rename in template to `category.label` / `vep_plugin.label` (Django enums expose `.label`). Cleanest.
- **(b)** Add `get_category_display = lambda self: self.category.label` properties to the dataclass for zero template changes.

Vote: **(a)** — small touch-up, more honest.

---

## 4. Migrations

### 4.1 Model + table — stay in place for now

**No `DeleteModel` migration in this refactor.** The `ColumnVEPField` model definition and its data stay untouched. Application code stops reading the table after PR1 lands, but the data is preserved as a passive backup while production validates the registry. Removal is a separate follow-up once Dave is satisfied (see §8 decision 5).

### 4.2 Existing data migrations (`0094`, `0115`, `0121`, `0122`, …)

**Stay untouched.** They continue to populate `ColumnVEPField` on fresh DBs as before. Harmless — nothing reads from the table at runtime once PR1 lands.

### 4.3 Bootstrap of `vep_columns.py` content

> **Authoring rule for the assistant — strict:**
>
> **First pass = dumb transcription. No collapsing. No cleverness.**
>
> Walk `column_vep_field.json` top to bottom and emit **one `VEPColumnDef` per JSON row, in JSON order**. Field-for-field copy. Single-element `genome_builds=frozenset({"GRCh37"})` etc. — no merging across rows yet. The synthetic `column` field is dropped (verified unread in §1.3); everything else transcribes 1:1.
>
> The result will be 237 entries, the same shape as the dump. **That is the goal of pass 1.** Run `manage.py check_vep_columns_equivalence` (§4.4). Only when it returns OK is pass 1 complete.
>
> **Pass 2 (separate commit) = cleanup.** Only then do we collapse multi-build duplicates, fold multi-destination groups (Mastermind, SV_overlap → gnomad_af, Existing_variation → cosmic+dbsnp), and reorganise into logical sections. After every collapse, re-run the equivalence check; it must stay green. The commit history then reads as "transcribe → verify → collapse → verify", which is the audit trail Dave is asking for.

So the two-stage shape:

| Pass | Goal | Entry count | Equivalence check |
|---|---|---|---|
| 1 | Byte-for-byte transcription, JSON order preserved | 237 | Must pass |
| 2 | Collapse duplicates, group logically | ~140 | Must still pass after each collapse |

### 4.4 Equivalence check — one-off management command

Run after **every** change to `vep_columns.py` — pass 1 (transcription) and every collapse step in pass 2. The command loads `column_vep_field.json` *and* the new `VEP_COLUMNS` registry and asserts they describe the same set of mappings. Every JSON row must round-trip back out of the registry, byte-identical in meaning, before the JSON is allowed to be deleted.

```python
# annotation/management/commands/check_vep_columns_equivalence.py
import json
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from annotation.vep_columns import VEP_COLUMNS


class Command(BaseCommand):
    help = ("One-off: load column_vep_field.json and verify the new VEP_COLUMNS registry "
            "describes an identical set of (source_field, vep_plugin, vep_custom, "
            "source_field_has_custom_prefix, category, processing_description, summary_stats, "
            "version-bounds, genome_build, pipeline_type, variant_grid_column) tuples. "
            "Delete me + column_vep_field.json once green.")

    def add_arguments(self, parser):
        parser.add_argument("--json", default="column_vep_field.json",
                            help="Path to the JSON dump of the old ColumnVEPField table.")

    def handle(self, *args, **opts):
        json_path = Path(opts["json"])
        if not json_path.exists():
            raise CommandError(f"{json_path} not found")

        rows = json.loads(json_path.read_text())

        # ---- expand the JSON: each row is one (mapping, build, pipeline, vgc) tuple ----
        json_tuples = set()
        for r in rows:
            json_tuples.add((
                r["source_field"],
                r["vep_plugin"],
                r["vep_custom"],
                r["source_field_has_custom_prefix"],
                r["category"],
                r["source_field_processing_description"],
                r["summary_stats"],
                r["min_columns_version"],
                r["max_columns_version"],
                r["min_vep_version"],
                r["max_vep_version"],
                r["genome_build_id"],          # may be None
                r["pipeline_type"],            # may be None
                r["variant_grid_column_id"],
            ))

        # ---- expand the registry the same way ----
        registry_tuples = set()
        for c in VEP_COLUMNS:
            builds = sorted(c.genome_builds) or [None]
            pipelines = [p.value for p in sorted(c.pipeline_types, key=lambda p: p.value)] or [None]
            for build in builds:
                for pipeline in pipelines:
                    for vgc in c.variant_grid_columns:
                        registry_tuples.add((
                            c.source_field,
                            c.vep_plugin.value if c.vep_plugin else None,
                            c.vep_custom.value if c.vep_custom else None,
                            c.source_field_has_custom_prefix,
                            c.category.value,
                            c.source_field_processing_description,
                            c.summary_stats,
                            c.min_columns_version,
                            c.max_columns_version,
                            c.min_vep_version,
                            c.max_vep_version,
                            build,
                            pipeline,
                            vgc,
                        ))

        only_in_json = json_tuples - registry_tuples
        only_in_registry = registry_tuples - json_tuples

        if not only_in_json and not only_in_registry:
            self.stdout.write(self.style.SUCCESS(
                f"OK: {len(json_tuples)} mappings match across JSON and registry."
            ))
            return

        if only_in_json:
            self.stdout.write(self.style.ERROR(
                f"\n{len(only_in_json)} mapping(s) in JSON but missing from registry:"))
            for t in sorted(only_in_json, key=lambda x: tuple(str(v) for v in x)):
                self.stdout.write(f"  {t}")
        if only_in_registry:
            self.stdout.write(self.style.ERROR(
                f"\n{len(only_in_registry)} mapping(s) in registry but missing from JSON:"))
            for t in sorted(only_in_registry, key=lambda x: tuple(str(v) for v in x)):
                self.stdout.write(f"  {t}")

        raise CommandError("Registry does not match JSON. See diff above.")
```

The semantics: explode both sides into the same flat tuple representation — one tuple per `(mapping × build × pipeline × destination)` — and compare as sets. This is exactly the "are these two configurations identical?" question and surfaces every kind of drift (missing build, missing destination, wrong category, etc.).

Workflow:
1. **Pass 1 — dumb transcribe.** Build `vep_columns.py` as one entry per JSON row, JSON order. Run `python3 manage.py check_vep_columns_equivalence`. Must return OK before doing anything else. This is its own commit.
2. **Pass 2 — collapse.** Each collapse (per-build merge, multi-destination merge, logical reorganise) is a separate commit. After every commit, re-run the equivalence check. It must stay green.
3. Land PR1 once pass 2 is complete and the call sites have been swapped over.
4. Once PR1 is in production, Dave deletes `column_vep_field.json`. The management command becomes dead code; PR2 (DeleteModel) also deletes the command file.

I'll commit the bootstrap script as `scripts/one_offs/dump_to_vep_columns.py` (read-only, never re-run after the migration lands).

### 4.5 Roll-out

**Single PR.** Add `annotation/vep_columns.py` + validation hook + unit tests + `check_vep_columns_equivalence` command. Keep `ColumnVEPField` model + table + data untouched. Switch all 6 call sites to read from the registry. Deploy. Verify in production (descriptions page renders, annotation runs work, custom columns grid loads).

A second cleanup PR (delete model + dump + commands + bootstrap script) happens **later, on Dave's signal**, once production confirms the registry is correct. Tracked as a follow-up note (§8 decision 5), not part of this refactor.

---

## 5. Anomalies — outcomes

| Anomaly | Today | After |
|---|---|---|
| Per-build duplicates (Aloft, dbNSFP, …) | 64 duplicate rows across 37/38/T2T | One `VEPColumnDef` with `genome_builds=frozenset({...})` (or empty = all) |
| `Mastermind_counts` (3 × 2) | 6 rows | 1 entry, `variant_grid_columns=(mm1, mm2, mm3)`, `genome_builds={GRCh37, GRCh38}` |
| `gnomad_sv_overlap_copied_to_gnomad_af` placeholder | 1 fake row | Folded: SV_overlap_AF entry's `variant_grid_columns=("gnomad_sv_overlap_af", "gnomad_af")` |
| `Existing_variation` → cosmic + dbsnp | 2 rows | 1 entry, two destinations |
| Empty source SV_overlap (`coords` / `percent`) | 2 rows | 1 entry, two destinations |
| Adding T2T or VEP 113 | New `_one_off` migration | Edit one Python file |

---

## 6. Risks

1. **Loss of FK integrity from `VariantGridColumn`** — replaced by import-time + test-time validation that every `variant_grid_columns` entry exists in the DB. Same effective guarantee, slightly later in the boot cycle. (`apps.ready()` doesn't have DB access guaranteed; we run the check in the management `system_check` framework instead, which `manage.py check` and CI both invoke.)
2. **Editor workflow** — these aren't edited by anyone outside dev; deferring to a PR is correct for this kind of registry.
3. **Migration `0094` reverse path** uses `vep_plugin=...` to delete by query against the old model. Still works — historical migrations operate on historical state.
4. **`SVOverlapProcessor` semantics** — see §3.4. After collapsing the `_copied_to_gnomad_af` anomaly, `gnomad_af` now appears in `self.sv_fields` for the SV pipeline, which matches what the code already does (it writes the chosen overlap's `gnomad_af` value). Worth a smoke test against a CNV annotation run, but no logic change.
5. **`columnvepfield_set` reverse accessor** in `models_columns.py:50` disappears with the model. Migration §3.6 swaps it for `for_variant_grid_column(self.pk)`.
6. **Two-PR window** — between PR1 and PR2 the DB still has a `ColumnVEPField` table that's no longer read. Harmless. Don't merge `0129` until PR1 is verified in production.

---

## 7. Test plan

1. **Static** — `python3 manage.py check` exercises the validation. Add `annotation/tests/test_vep_columns.py` that asserts:
   - Every referenced `variant_grid_column` exists.
   - Every `VEPColumnDef` is unique under the dedupe key.
   - For each of (GRCh37, GRCh38, T2T-CHM13v2.0) × (STANDARD, STRUCTURAL_VARIANT) × (columns_version=2, 3) × (vep_version=110, 111, 112), the resolved set of `(source_field → variant_grid_columns)` matches the snapshot from `column_vep_field.json` filtered through the equivalent SQL `Q` filters. This is the **regression net** for the bootstrap.
2. **Inserter** — run existing `annotation/tests/` suite with `--keepdb`. The `bulk_vep_vcf_annotation_inserter` tests are the most sensitive.
3. **Descriptions view** — render `view_annotation_descriptions` for each build in dev and visually diff against the current page.
4. **CNV pipeline** — re-run a known SV annotation against staging and diff `gnomad_sv_overlap_*` + `gnomad_af` columns.

---

## 8. Decisions

1. **`genome_builds=frozenset()` means "all builds".** If a future build needs per-build tweaks, we make those edits then. Matches existing `null` semantics (~52 rows already use it).
2. **Template rename.** Swap `get_category_display` / `get_vep_plugin_display` → `.label` in `view_annotation_descriptions.html` (§3.7 option a). No shim on the dataclass.
3. **Bootstrap script lifecycle.** `scripts/one_offs/dump_to_vep_columns.py` is committed alongside its output for auditability. Once everything is done (registry green, PR1 + PR2 landed, JSON dump deleted), the bootstrap script gets deleted too.
4. **`cvf.column` was a human-readable handle for management scripts** (lets you do `ColumnVEPField.objects.get(column="aloft_pred_38").update(...)` instead of pasting PKs). Acceptable to lose — management edits now happen by editing `vep_columns.py` directly.
5. **Don't `DeleteModel` yet.** Skip §4.1 and PR2 for now. The model and table stay in place, ignored by application code, while testing happens in production. The `column_vep_field.json` dump also stays put.

   **Follow-up note (file in `claude/plans/` once testing is green):**
   > After production has run on the registry-only path for long enough to trust it: (a) delete `column_vep_field.json`; (b) delete `annotation/management/commands/check_vep_columns_equivalence.py`; (c) delete `scripts/one_offs/dump_to_vep_columns.py`; (d) write `0129_delete_columnvepfield.py` and remove `ColumnVEPField` from `annotation/models/models.py`.

---

## 9. Files touched

**New**
- `annotation/vep_columns.py` (registry + helpers; transcribed record-by-record from `column_vep_field.json`)
- `annotation/tests/test_vep_columns.py` (validation + snapshot tests)
- `annotation/management/commands/check_vep_columns_equivalence.py` (one-off equivalence check vs. JSON dump; deleted in cleanup follow-up)
- `scripts/one_offs/dump_to_vep_columns.py` (bootstrap; deleted in cleanup follow-up)

**Modified**
- `annotation/vep_annotation.py` — 3 call sites swapped to registry
- `annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py` — inserter setup + `SVOverlapProcessor`
- `annotation/views.py` — descriptions view
- `annotation/templates/annotation/view_annotation_descriptions.html` — `get_*_display` → `.label`
- `annotation/apps.py` — register system check
- `snpdb/models/models_columns.py` — `columns_version_description`
- `snpdb/grid_columns/custom_columns.py` — version-applicability `Q`

**Untouched (this refactor)**
- `annotation/models/models.py` — `ColumnVEPField` stays defined
- All historical `annotation/migrations/00**_*.py` that populate `ColumnVEPField`
- The `ColumnVEPField` table data — preserved as a passive backup

**Cleanup follow-up (after Dave signs off in production)**
- Delete `column_vep_field.json`
- Delete `annotation/management/commands/check_vep_columns_equivalence.py`
- Delete `scripts/one_offs/dump_to_vep_columns.py`
- Add `annotation/migrations/0129_delete_columnvepfield.py`
- Remove `ColumnVEPField` from `annotation/models/models.py`

