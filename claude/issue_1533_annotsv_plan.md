# Issue #1533 — Integrate AnnotSV into SV annotation pipeline

GitHub: https://github.com/SACGF/variantgrid/issues/1533
Parent (umbrella SV-annotation comparison): https://github.com/SACGF/variantgrid/issues/1040

## Goal

Enrich the `STRUCTURAL_VARIANT` annotation pipeline with the SV-specific
fields requested in the umbrella issue #1040 (regulatory elements,
repeats, SegDups, ENCODE blacklist, OMIM inheritance, gene morbidity,
ACMG/AMP-style SV ranking, etc.). The wishlist in #1040 corresponds
almost 1:1 to AnnotSV's TSV output schema, so this issue adds **AnnotSV
as a post-VEP stage** of the existing SV pipeline rather than try to
reimplement those columns via VEP `--custom` BED files.

ClassifyCNV is **already implemented inside AnnotSV** (it's the engine
behind `ACMG_class` / `AnnotSV_ranking_score`), so this iteration covers
ClassifyCNV-equivalent output as part of the AnnotSV stage rather than
as a separate tool.

Out of scope for this iteration (parked for follow-ups, none of which
ship with AnnotSV): CADD-SV, SVAFotate, DeepSVP. The integration shape
we choose for AnnotSV should make adding any of these as an additional
stage straightforward.

## Current SV pipeline (what's there today)

`AnnotationRun` (`annotation/models/models.py:852`) has a
`pipeline_type` enum
(`VariantAnnotationPipelineType` — `annotation/models/models_enums.py:73`)
with values `STANDARD` ("S") and `STRUCTURAL_VARIANT` ("C"). The
scheduler creates one `AnnotationRun` per `(annotation_range_lock,
pipeline_type)` combination
(`annotation/tasks/annotation_scheduler_task.py:46-57`).

`annotate_variants` (`annotation/tasks/annotate_variants.py:51`) for an
SV run does:

1. `dump_and_annotate_variants` — writes a VCF of unannotated symbolic
   variants in the range, runs VEP via `get_vep_command(...,
   pipeline_type=STRUCTURAL_VARIANT)`
   (`annotation/vep_annotation.py:143`). For the SV pipeline, all
   plugins are stripped (`vep_annotation.py:219-248`), and
   `--max_sv_size` is passed when configured.
2. `import_vcf_annotations` — parses VEP CSQ output via
   `BulkVEPVCFAnnotationInserter`
   (`annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py`). For
   SV runs it instantiates an `SVOverlapProcessor`
   (line 146-149) that computes overlap percent against the gnomAD-SV
   custom annotation BED and writes
   `gnomad_sv_overlap_af / _name / _percent / _coords` on
   `VariantAnnotation` (`annotation/models/models.py:1138-1143`).

So the SV pipeline today is essentially: **VEP (no plugins) + a
gnomAD-SV BED overlap pass**. Everything in the issue #1040 wishlist is
new territory.

## Proposed architecture

Add an **AnnotSV stage** that runs on the same VCF dump as VEP, in the
same `AnnotationRun`. The simplest framing:

```
dump VCF  ──►  VEP (existing)            ──►  vcf_annotated_filename
          └─►  AnnotSV (new)             ──►  annotsv_tsv_filename
                                              │
                          import_vcf_annotations (existing CSQ parse)
                                          +
                          import_annotsv_tsv (new, joined by VCF ID)
                                          ▼
                      VariantAnnotation rows (existing model, new fields)
```

Both stages consume the dump VCF and write into the same
`VariantAnnotation` row keyed off `Variant`. AnnotSV's "full" lines map
to whole-SV columns; "split" lines (per-overlapped-gene) map naturally
into `VariantTranscriptAnnotation` (or a new `VariantSVGeneOverlap`
model — see "Open question 1" below).

We do **not** convert AnnotSV TSV → VCF INFO. AnnotSV doesn't natively
emit VCF, and round-tripping TSV→VCF→ingest is more brittle than
reading the TSV directly. AnnotSV column names are stable across
releases (see `https://lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf`).

### Why a stage on the existing pipeline_type, not a new pipeline_type

Pros of keeping it in the same `AnnotationRun`:
- Single dump VCF, single set of variants, single completion signal.
- Reuses range-lock book-keeping, retry/error handling, status display.
- The wishlist columns logically belong on the same row as the gnomAD-SV
  overlap fields.

Cons:
- AnnotSV failures block the run unless we make the stage best-effort
  (recommended: log + continue on AnnotSV failure, set a flag on the
  `AnnotationRun` so it's visible in admin).
- A version bump of AnnotSV's annotation sources requires reannotation;
  this is the same problem we already have with VEP versions and is
  handled by `VariantAnnotationVersion` increments.

A separate `pipeline_type=ANNOTSV` would let us reannotate AnnotSV
independently of VEP, but at the cost of doubling SV runs and
complicating the admin UI. Recommendation: **stay with one
`STRUCTURAL_VARIANT` run, add stage state**. If we ever need
independent reannotation, splitting later is a straightforward
migration.

## Detailed changes

### 1. AnnotSV install + settings

AnnotSV ships as TCL scripts plus a large annotations bundle
(`Annotations_Human/`). Treat it like VEP — install once per host,
point at it via settings.

`variantgrid/settings/components/annotation_settings.py`:

```python
ANNOTATION_ANNOTSV_ENABLED = False  # opt-in per env
ANNOTATION_ANNOTSV_BIN = "/data/annotation/AnnotSV/bin/AnnotSV"          # the TCL entry point
ANNOTATION_ANNOTSV_ANNOTATIONS_DIR = "/data/annotation/AnnotSV/share/AnnotSV/Annotations_Human"
ANNOTATION_ANNOTSV_GENOME_BUILD = {                          # AnnotSV's build flag
    "GRCh37": "GRCh37",
    "GRCh38": "GRCh38",
}
ANNOTATION_ANNOTSV_EXTRA_ARGS: list[str] = []                # e.g. ["-SVminSize", "50"]
ANNOTATION_ANNOTSV_TIMEOUT_SECONDS = 60 * 60                 # SV volume is low; cap anyway
```

`-tx` is set per-run from the project's `AnnotationConsortium`
(`REFSEQ` → `RefSeq`, `ENSEMBL` → `ENSEMBL`) so AnnotSV's transcript
provider matches VEP's. Note that AnnotSV's transcript snapshot is
pinned by the annotations bundle and is **not synchronised with VEP's
cache version** — only the provider aligns. Per-gene rows therefore
key off `Gene` (not transcript), and the bundle version
(`annotsv_bundle`) acts as the version pin for that snapshot.

A version sentinel goes on `VariantAnnotationVersion` (see #4).

### 2. New stage: `run_annotsv`

New module `annotation/annotsv_annotation.py`:

```python
import logging
import os
import subprocess
from django.conf import settings
from snpdb.models import GenomeBuild


def get_annotsv_command(vcf_filename: str, output_dir: str,
                        genome_build: GenomeBuild,
                        annotation_consortium: AnnotationConsortium) -> list[str]:
    build_arg = settings.ANNOTATION_ANNOTSV_GENOME_BUILD[genome_build.name]
    tx_arg = "RefSeq" if annotation_consortium == AnnotationConsortium.REFSEQ else "ENSEMBL"
    cmd = [
        settings.ANNOTATION_ANNOTSV_BIN,
        "-SVinputFile", vcf_filename,
        "-outputDir", output_dir,
        "-genomeBuild", build_arg,
        "-annotationsDir", settings.ANNOTATION_ANNOTSV_ANNOTATIONS_DIR,
        "-tx", tx_arg,                # match VEP's transcript provider
        "-SVinputInfo", "1",          # keep INFO from input VCF
        "-includeCI", "0",            # don't expand by CIPOS/CIEND
        "-overwrite", "1",
    ]
    cmd.extend(settings.ANNOTATION_ANNOTSV_EXTRA_ARGS)
    return cmd


def run_annotsv(vcf_filename: str, output_dir: str,
                genome_build: GenomeBuild) -> tuple[str, int, str, str]:
    """ Returns (tsv_filename, return_code, stdout, stderr). """
    os.makedirs(output_dir, exist_ok=True)
    cmd = get_annotsv_command(vcf_filename, output_dir, genome_build)
    proc = subprocess.run(
        cmd,
        capture_output=True, text=True,
        timeout=settings.ANNOTATION_ANNOTSV_TIMEOUT_SECONDS,
        check=False,
    )
    # AnnotSV writes <basename>.annotated.tsv inside output_dir
    base = os.path.splitext(os.path.basename(vcf_filename))[0]
    tsv_filename = os.path.join(output_dir, f"{base}.annotated.tsv")
    if proc.returncode != 0 or not os.path.exists(tsv_filename):
        logging.warning("AnnotSV failed (rc=%s): %s", proc.returncode, proc.stderr[:2000])
    return tsv_filename, proc.returncode, proc.stdout, proc.stderr
```

Hook it into `dump_and_annotate_variants`
(`annotation/tasks/annotate_variants.py:104`) **after** the successful
VEP block, gated on `ANNOTATION_ANNOTSV_ENABLED` and
`pipeline_type == STRUCTURAL_VARIANT`. AnnotSV failure is logged onto
the run but does not raise — the VEP-only result is still useful.

```python
if (settings.ANNOTATION_ANNOTSV_ENABLED
        and annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT):
    annotsv_dir = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR,
                               f"annotsv_{annotation_run.pk}")
    tsv, rc, _, stderr = run_annotsv(vcf_dump_filename, annotsv_dir, genome_build)
    if rc == 0 and os.path.exists(tsv):
        annotation_run.annotsv_tsv_filename = tsv
    else:
        annotation_run.annotsv_error = (stderr or "")[:100_000]
    annotation_run.save()
```

### 3. Schema changes on `AnnotationRun`

Add on `annotation/models/models.py:852`:

```python
# stage tracking for AnnotSV
annotsv_tsv_filename = models.TextField(null=True)
annotsv_error = models.TextField(null=True)
annotsv_imported = models.BooleanField(default=False)
```

Update `get_status` (line 917) to surface AnnotSV-only failure in the
admin grid (e.g. a new `AnnotationStatus.PARTIAL` value), or — simpler
— leave overall status driven by VEP and expose AnnotSV state as its
own column in `AnnotationRunColumns` (`annotation/grids.py:13`).

### 4. New annotation version fields

`VariantAnnotationVersion` already gates VEP/columns versions. Add two
AnnotSV version pins — kept **separate and explicit** because the
binary and the annotations bundle have independent release cadences and
upgrades will be rare; an explicit mismatch is preferable to a single
collapsed version string:

```python
# annotation/models/models.py — VariantAnnotationVersion
annotsv_code = models.TextField(null=True, blank=True)     # e.g. "3.5.8" — `AnnotSV -version`
annotsv_bundle = models.TextField(null=True, blank=True)   # annotations bundle stamp
```

**Backwards-compatible defaults.** Both fields default to `null` /
unset, and `ANNOTATION_ANNOTSV_ENABLED=False` in the shipped defaults.
Upgrading VariantGrid alone — without installing AnnotSV or flipping
the enable flag — does **not** change `VariantAnnotationVersion` and
therefore does **not** trigger reannotation of existing variants. The
version pins only become populated (and only become reannotation
triggers) on deployments that explicitly opt in by:

1. Installing AnnotSV + the annotations bundle.
2. Setting `ANNOTATION_ANNOTSV_ENABLED=True` and the path settings in
   the deployment-specific settings file.
3. Running the management command that creates a new
   `VariantAnnotationVersion` populating `annotsv_code` and
   `annotsv_bundle` from the installed binary/bundle.

`vep_check_command_line_version_match` has a counterpart
`annotsv_check_command_line_version_match` that runs `AnnotSV -version`
and compares with `VariantAnnotationVersion.annotsv_code`. The bundle
version is read from the annotations directory (AnnotSV ships a
`Annotations_Human/Users/configfile` / version stamp; the importer
captures the string and compares against `annotsv_bundle`). The check
is **skipped entirely** when `ANNOTATION_ANNOTSV_ENABLED=False` so a
default install never tries to call a binary that isn't there. Once
opted in, bumping either field triggers reannotation via the standard
`VariantAnnotationVersion` path.

### 5. Schema changes on `VariantAnnotation` (whole-SV columns)

These are the AnnotSV "full" columns from issue #1040. Add to
`AbstractVariantAnnotation` or to `VariantAnnotation` directly — they
are SV-only so adding to `VariantAnnotation` (and not the per-transcript
table) keeps `STANDARD` rows from carrying NULLs they'll never use.

```python
# annotation/models/models.py — VariantAnnotation (SV-only block,
# adjacent to existing gnomad_sv_overlap_* fields)
annotsv_acmg_class = models.IntegerField(null=True, blank=True)        # 1..5
annotsv_acmg_score = models.FloatField(null=True, blank=True)

# Regulatory / dark / repeat / segdup region annotations
re_gene = models.TextField(null=True, blank=True)
repeat_type_left = models.TextField(null=True, blank=True)
repeat_type_right = models.TextField(null=True, blank=True)
segdup_left = models.TextField(null=True, blank=True)
segdup_right = models.TextField(null=True, blank=True)
encode_blacklist_left = models.TextField(null=True, blank=True)
encode_blacklist_right = models.TextField(null=True, blank=True)

# Population freq (AnnotSV bundles gnomAD-SV / 1kGP / IMH / DGV)
annotsv_gnomad_pmax_af = models.FloatField(null=True, blank=True)
annotsv_1kgp_max_af = models.FloatField(null=True, blank=True)
annotsv_imh_max_af = models.FloatField(null=True, blank=True)
```

If/when CoLoRSdb is added, it slots in here as
`colorsdb_af / _ac / _ac_hemi / _nhomalt / _maxaf`. Same for
DeepSVP's `del_z / dup_z / cnv_z`.

### 6. Per-gene SV overlap rows (AnnotSV "split" lines) — **deferred**

AnnotSV emits one **full** line per SV, then N **split** lines, one per
overlapped gene, each carrying gene-scoped morbidity / OMIM /
inheritance / exons-spanned / frameshift columns.

**Not in this iteration.** AnnotSV's gene/transcript snapshot is
bundle-pinned, so its gene symbols and IDs come from a different release
than the project's `Gene` table (driven by VEP / cdot). Persisting split
rows now would require a gene-release mapping layer we don't have, and
would risk dangling FKs / silent symbol drift.

For this iteration we ingest **only the "full" lines** onto
`VariantAnnotation`. Split lines are parsed for diagnostic logging but
not stored. A follow-up issue will design the gene-release mapping and
introduce a `VariantSVGeneOverlap`-equivalent model once the mapping
strategy is settled.

### 7. TSV ingestion

New module `annotation/vcf_files/bulk_annotsv_tsv_inserter.py`. Mirrors
the design of `BulkVEPVCFAnnotationInserter` but reads a TSV with
pandas (or `csv.DictReader`) and bulk-inserts via existing
`upsert`/`bulk_create` paths:

```python
import pandas as pd
from annotation.models import VariantAnnotation


# AnnotSV TSV → our model field name. Names below are the AnnotSV
# README v3.x labels — verify exact casing on the deployed bundle.
FULL_COLUMN_MAP = {
    "ACMG_class": "annotsv_acmg_class",
    "AnnotSV_ranking_score": "annotsv_acmg_score",
    "RE_gene": "re_gene",
    "Repeat_type_left": "repeat_type_left",
    "Repeat_type_right": "repeat_type_right",
    "SegDup_left": "segdup_left",
    "SegDup_right": "segdup_right",
    "ENCODE_blacklist_characteristics_left": "encode_blacklist_left",
    "ENCODE_blacklist_characteristics_right": "encode_blacklist_right",
    "GnomAD_pLI": None,        # already on Gene
    "B_gain_AFmax": "annotsv_gnomad_pmax_af",  # depends on bundle
    "1000g_max_AF": "annotsv_1kgp_max_af",
    "IMH_AF": "annotsv_imh_max_af",
}

# Split lines are NOT persisted in this iteration — AnnotSV's gene
# snapshot is on a different release than our Gene table, and we don't
# have a mapping yet. Defer until follow-up.


def import_annotsv_tsv(annotation_run):
    if not annotation_run.annotsv_tsv_filename:
        return
    df = pd.read_csv(annotation_run.annotsv_tsv_filename, sep="\t",
                     dtype=str, na_values=["", "NA", "."])

    # AnnotSV ID format encodes the input VCF ID; we wrote IDs as
    # variant_id during dump (see _unannotated_variants_to_vcf), so
    # we can join straight back to Variant.pk.
    full = df[df["Annotation_mode"] == "full"]
    _bulk_update_full(annotation_run, full)
    annotation_run.annotsv_imported = True
    annotation_run.save()
```

Called from `import_vcf_annotations` immediately after the VEP CSQ
import path completes, and only when
`annotation_run.annotsv_tsv_filename` is set.

### 8. UI surfacing

Existing variant detail templates already render the
`gnomad_sv_overlap_*` fields when present. Extend the SV detail panel
to display:
- AnnotSV ACMG class + score (with the standard 1–5 colour scale).
- Region columns (RE_gene, repeats, SegDup, ENCODE blacklist) in a
  collapsible "Genomic context" section.

Per-gene split-line rows are deferred — see §6.

`annotation/templates/annotation/variantannotation_detail.html` (or
the SV-specific partial — confirm location during implementation) is
the touchpoint. No new datatable views required for the first cut.

## Migration & rollout

1. **Schema** — one Django migration adds AnnotSV fields to
   `AnnotationRun`, `VariantAnnotationVersion`, and `VariantAnnotation`.
   All new fields are nullable / default off, so an upgrade with no
   AnnotSV install is a no-op for existing rows. No new per-gene table
   this iteration (deferred — see §6).
2. **Settings** — `ANNOTATION_ANNOTSV_ENABLED=False` in default,
   opt-in per env. Shariant / SAP can enable when the bundle is
   installed; vgtest stays off until verified. **A VG version bump
   alone must not flip any AnnotSV setting on** — the enable flag and
   path overrides only ever come from the deployment-specific settings
   file, never from the shipped defaults.
3. **Backfill** — `annotsv_code` / `annotsv_bundle` start `null` and
   only become populated when a deployment opts in and creates a new
   `VariantAnnotationVersion`. Until then, no reannotation is
   triggered. Once opted in, bumping either field triggers reannotation
   the standard way. SV volume is low (orders of magnitude fewer than
   SNVs), so a full backfill is cheap.
4. **Failure isolation** — AnnotSV stage is best-effort. A failed
   AnnotSV stage leaves `annotsv_imported=False` and an error string;
   ops can rerun that stage alone via a management command:

   ```
   python3 manage.py annotsv_run --annotation-run <id>
   ```

   (Thin wrapper that re-invokes `run_annotsv` + `import_annotsv_tsv`
   for one run.)
5. **Tests** — fixture TSV in `annotation/tests/test_data/annotsv/`,
   unit test for `bulk_annotsv_tsv_inserter` parsing both full and
   split rows, and a smoke test in `annotation/tests/test_annotsv.py`
   that mocks `subprocess.run` and asserts the importer wires data to
   the right `Variant` rows.

## Decisions (resolved 2026-04-29)

1. **Per-gene overlaps storage** — **deferred**. AnnotSV's split lines
   reference a bundle-pinned gene release that doesn't match our `Gene`
   table; persisting them needs a gene-release mapping we don't have
   yet. This iteration ingests only the AnnotSV "full" lines onto
   `VariantAnnotation`. A follow-up issue will design the mapping and
   introduce a dedicated `VariantSVGeneOverlap`-style table — kept
   separate from `VariantGeneOverlap` so AnnotSV vs VEP overlap calls
   can be compared.
2. **ACMG class** — annotation-only. Stored on `VariantAnnotation`,
   surfaced in the variant detail UI; **not** auto-populating any
   classification evidence keys. Classification app can pull it
   explicitly later if curators want it.
3. **Score thresholds vs raw values** — store the raw
   `AnnotSV_ranking_score`. UI does any bucketing.
4. **Versioning** — two explicit fields, `annotsv_code` (binary) and
   `annotsv_bundle` (annotations directory). Upgrades are rare; an
   explicit mismatch is preferable to a single collapsed version
   string. Both default to `null` so a VG version bump on a deployment
   that hasn't installed AnnotSV does not change
   `VariantAnnotationVersion` and does not trigger reannotation —
   AnnotSV is **strictly opt-in via deployment settings**.
5. **Transcript provider** — pass `-tx` matching the project's
   `AnnotationConsortium` so the provider aligns with VEP. Versions
   inside the AnnotSV bundle are independent of the VEP cache and will
   drift; per-gene rows therefore key off `Gene`, not transcript.

## Future stages (not in this iteration)

ClassifyCNV is dropped from this list — AnnotSV's `ACMG_class` /
`AnnotSV_ranking_score` are produced by ClassifyCNV's criteria
internally, so this iteration already delivers it.

Once the AnnotSV stage exists, additional external tools drop in by
repeating steps 2/5/7:
- **CADD-SV** — single-score column.
- **SVAFotate** — alternate population AF source (CoLoRSdb, etc.) keyed
  by overlap.
- **DeepSVP** — `del_z / dup_z / cnv_z`.

Each is a small, independent stage on the SV `AnnotationRun`, gated by
its own `ANNOTATION_<TOOL>_ENABLED` setting.
