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

Out of scope for this iteration (parked for follow-ups):
CADD-SV, SVAFotate, ClassifyCNV, DeepSVP. The integration shape we
choose for AnnotSV should make adding any of these as an additional
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
ANNOTATION_ANNOTSV_BIN = "/opt/AnnotSV/bin/AnnotSV"          # the TCL entry point
ANNOTATION_ANNOTSV_ANNOTATIONS_DIR = "/data/AnnotSV/Annotations_Human"
ANNOTATION_ANNOTSV_GENOME_BUILD = {                          # AnnotSV's build flag
    "GRCh37": "GRCh37",
    "GRCh38": "GRCh38",
}
ANNOTATION_ANNOTSV_EXTRA_ARGS: list[str] = []                # e.g. ["-SVminSize", "50"]
ANNOTATION_ANNOTSV_TIMEOUT_SECONDS = 60 * 60                 # SV volume is low; cap anyway
```

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
                        genome_build: GenomeBuild) -> list[str]:
    build_arg = settings.ANNOTATION_ANNOTSV_GENOME_BUILD[genome_build.name]
    cmd = [
        settings.ANNOTATION_ANNOTSV_BIN,
        "-SVinputFile", vcf_filename,
        "-outputDir", output_dir,
        "-genomeBuild", build_arg,
        "-annotationsDir", settings.ANNOTATION_ANNOTSV_ANNOTATIONS_DIR,
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

### 4. New annotation version field

`VariantAnnotationVersion` already gates VEP/columns versions. Add an
AnnotSV version pin so a bundle upgrade triggers reannotation the same
way a VEP version bump does:

```python
# annotation/models/models.py — VariantAnnotationVersion
annotsv_version = models.TextField(null=True, blank=True)
annotsv_annotations_version = models.TextField(null=True, blank=True)
```

`vep_check_command_line_version_match` has a counterpart
`annotsv_check_command_line_version_match` that runs `AnnotSV --version`
and compares with `VariantAnnotationVersion.annotsv_version`.

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

### 6. Per-gene SV overlap rows (AnnotSV "split" lines)

AnnotSV emits one **full** line per SV, then N **split** lines, one per
overlapped gene, each carrying gene-scoped morbidity / OMIM /
inheritance / exons-spanned / frameshift columns. These don't slot well
into `VariantTranscriptAnnotation` (which is per-transcript, not
per-gene-overlap-of-an-SV) — recommend a new model:

```python
# annotation/models/models.py
class VariantSVGeneOverlap(models.Model):
    version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    annotation_run = models.ForeignKey(AnnotationRun, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, null=True, on_delete=SET_NULL)
    gene_symbol = models.TextField(null=True, blank=True)

    overlap_type = models.CharField(max_length=16, null=True, blank=True)  # full / partial
    exons_spanned = models.IntegerField(null=True, blank=True)
    frameshift = models.BooleanField(null=True)
    nearest_ss_type = models.CharField(max_length=2, null=True, blank=True)  # 3' / 5'

    omim_id = models.TextField(null=True, blank=True)
    omim_inheritance = models.TextField(null=True, blank=True)
    omim_morbid = models.BooleanField(null=True)

    class Meta:
        indexes = [models.Index(fields=["variant", "version"])]
```

Sharded the same way other annotation tables are
(`annotation_version_querysets.get_queryset_for_annotation_version`).
`AnnotationRun.delete_related_objects` (`models.py:947`) gets the new
class added to its loop:

```python
for klass in [VariantAnnotation, VariantTranscriptAnnotation,
              VariantGeneOverlap, VariantSVGeneOverlap]:
    qs = get_queryset_for_annotation_version(klass, annotation_version)
    qs.filter(annotation_run=self).delete()
```

### 7. TSV ingestion

New module `annotation/vcf_files/bulk_annotsv_tsv_inserter.py`. Mirrors
the design of `BulkVEPVCFAnnotationInserter` but reads a TSV with
pandas (or `csv.DictReader`) and bulk-inserts via existing
`upsert`/`bulk_create` paths:

```python
import pandas as pd
from annotation.models import VariantAnnotation, VariantSVGeneOverlap


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

SPLIT_COLUMN_MAP = {
    "Gene_name": "gene_symbol",
    "Overlapped_tx_length": None,     # derived
    "Tx_start": None,
    "Tx_end": None,
    "Location": "overlap_type",
    "Frameshift": "frameshift",
    "Exon_count": None,
    "Dist_nearest_SS": None,
    "Nearest_SS_type": "nearest_ss_type",
    "OMIM_ID": "omim_id",
    "OMIM_inheritance": "omim_inheritance",
    "OMIM_morbid": "omim_morbid",
}


def import_annotsv_tsv(annotation_run):
    if not annotation_run.annotsv_tsv_filename:
        return
    df = pd.read_csv(annotation_run.annotsv_tsv_filename, sep="\t",
                     dtype=str, na_values=["", "NA", "."])

    # AnnotSV ID format encodes the input VCF ID; we wrote IDs as
    # variant_id during dump (see _unannotated_variants_to_vcf), so
    # we can join straight back to Variant.pk.
    full = df[df["Annotation_mode"] == "full"]
    splits = df[df["Annotation_mode"] == "split"]

    _bulk_update_full(annotation_run, full)
    _bulk_create_splits(annotation_run, splits)
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
- Per-gene rows from `VariantSVGeneOverlap`.
- Region columns (RE_gene, repeats, SegDup, ENCODE blacklist) in a
  collapsible "Genomic context" section.

`annotation/templates/annotation/variantannotation_detail.html` (or
the SV-specific partial — confirm location during implementation) is
the touchpoint. No new datatable views required for the first cut.

## Migration & rollout

1. **Schema** — one Django migration adds AnnotSV fields to
   `AnnotationRun`, `VariantAnnotationVersion`, and `VariantAnnotation`,
   plus the new `VariantSVGeneOverlap` partitioned table.
2. **Settings** — `ANNOTATION_ANNOTSV_ENABLED=False` in default,
   opt-in per env. Shariant / SAP can enable when the bundle is
   installed; vgtest stays off until verified.
3. **Backfill** — bumping `VariantAnnotationVersion.annotsv_version`
   will trigger reannotation the standard way. SV volume is low
   (orders of magnitude fewer than SNVs), so a full backfill is cheap.
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

## Open questions

1. **Per-gene overlaps storage** — new `VariantSVGeneOverlap` model
   (recommended above) vs reusing `VariantGeneOverlap` (which exists
   today for SNV gene proximity). The latter may be ambiguous — SV
   "overlap" semantics differ from SNV `--overlap`-style proximity.
   New model is cleaner; happy to revisit.
2. **ACMG class surfacing in classification app** — should AnnotSV's
   ACMG class auto-populate any classification evidence keys, or stay
   purely in `annotation`? My instinct: annotation only; the
   classification app pulls from it explicitly when curators want it.
3. **Score thresholds vs raw values** — store the raw
   `AnnotSV_ranking_score` and let the UI bucket, rather than baking
   thresholds into the schema.
4. **Bundle versioning** — AnnotSV's annotation directory has its own
   release cadence. Capturing the bundle version (e.g. via a stamp
   file in the annotations dir) is worth doing in
   `annotsv_annotations_version` so reannotation triggers correctly
   when sources update without the binary changing.

## Future stages (not in this iteration)

Once the AnnotSV stage exists, additional stages drop in by repeating
steps 2/5/7:
- **CADD-SV** — single-score column.
- **SVAFotate** — alternate population AF source (CoLoRSdb, etc.) keyed
  by overlap.
- **ClassifyCNV** — second ACMG/AMP scorer for CNVs specifically; store
  alongside `annotsv_acmg_class` so curators can compare.
- **DeepSVP** — `del_z / dup_z / cnv_z`.

Each is a small, independent stage on the SV `AnnotationRun`, gated by
its own `ANNOTATION_<TOOL>_ENABLED` setting.
