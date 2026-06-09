# Issue #1568 — External annotation runs: store/dump VCFs for reuse

**Status:** Design / plan (no code yet)
**Issue:** https://github.com/SACGF/variantgrid/issues/1568 — milestone *SA Path VG 4*
**Author of issue:** davmlaw

---

## 1. Goal

Decouple the heavy VEP step from VariantGrid so annotation compute can run on separate
(non-VM) machines, **and** reuse the resulting annotated VCFs to avoid running VEP twice for the
same deploy.

Driving scenario (from the issue): upgrade `sapath test` (a clone of prod, needs a full
annotation rerun) → then upgrade prod. We want prod to reuse the annotated VCFs produced
during the test run instead of re-running VEP.

### Scope (explicit)

**In scope — reuse only between a database and its own clone, for identical annotation runs.**
The only supported reuse path is: a DB is cloned, both copies need the *same* annotation version, and
annotated VCFs produced on one are imported into the other for **annotation runs that are the same**
(same annotation-version identity, and variant ranges that still line up because one DB is a clone of
the other). Reuse across unrelated databases, or across diverged regions where variant ids no longer
line up, is **not** supported — the §6a id-alignment check exists precisely to reject those cases by
failing the import.

**Out of scope (for this plan):** a general "external compute for any DB" feature; coordinate-based
re-resolution of variants so divergent DBs can share annotation; automated cross-DB run matching that
tolerates id drift. The operator manually selects which annotated files to reuse (those predating any
divergence); VG's job is to make that safe, not to discover the mapping.

### Trigger: management command, not a settings flag

Kicking off externally-managed annotation runs is driven entirely by **one management command**, not by
a global setting. The *same* command both **dumps** (write VCFs + metadata + Snakemake bundle, park runs
in a waiting state) and **imports** (re-import annotated VCFs), selected by different params/subcommands.
This keeps external annotation an explicit operator action on a chosen build/version, rather than a mode
that silently changes the normal scheduler's behaviour.

Issue checklist, restated:

1. Self-describing dump filenames (server, version, run, range-lock min/max as a double check).
2. Runs can be dumped, then parked in an **"external annotation"** waiting state — until told where annotated files are.
3. Operator copies dump files to another machine, annotates, copies results back.
4. The dump command also generates a **Snakemake** workflow using the current VEP command, with a YAML config for VEP paths (on the compute server) + the dump directory.
5. The same command, in import mode, takes annotated VCF filenames, matches them back to annotation runs, does upload-only, sets states.
6. **Clone reuse** validated by the §6a metadata + variant-id-alignment sanity check.

**Decided (this plan):** run identity keys on **variant coordinates + annotation-version
identity** (genome build, VEP version, columns_version, consortium, gene-annotation release), and the
hard guarantee comes from the §6a per-record variant-id-alignment check.

---

## 2. Current architecture (what we build on)

The pipeline already splits cleanly into **dump → VEP → import**, each phase guarded by a field check.

| Concern | Location | Notes |
|---|---|---|
| Orchestrator task | `annotation/tasks/annotate_variants.py:51` `annotate_variants()` | Dumps only if `vcf_annotated_filename is None` (`:76`); imports only if it is set (`:79`). **This guard is the natural seam for "external" mode.** |
| Dump + local VEP | `annotate_variants.py:105` `dump_and_annotate_variants()` | Writes input VCF, builds + runs VEP, sets `vcf_annotated_filename`. |
| Dump writer | `annotate_variants.py:219` `_unannotated_variants_to_vcf()` / `write_qs_to_vcf()` | Writes variants in range to VCF, `variant_id` in INFO. |
| Upload-only retry | `annotate_variants.py:185` `annotation_run_retry(upload_only=True)` | **Already imports an existing annotated VCF** — the seed for the import command. |
| Importer | `annotation/vcf_files/import_vcf_annotations.py:18` | Reads annotated VCF, bulk-inserts; runs `vep_check_annotated_file_version_match()` first. |
| **Annotated-file version self-check** | `annotation/vep_annotation.py:494` `vep_check_annotated_file_version_match()` | Reads `##VEP=...` header from the annotated VCF (`get_vep_version_from_vcf()` `:452`) and asserts it matches the `VariantAnnotationVersion`. **Key cross-DB safety primitive.** |
| VEP command builder | `annotation/vep_annotation.py:100` `get_vep_command()` | Returns the full arg list — reuse for Snakemake generation. |
| Status (derived) | `annotation/models/models.py:1066` `get_status()` | Status computed from timestamp fields in `save()`; no manual writes. |
| Dump filename | `models.py:1104` `get_dump_filename()` | `dump_{pk}_{type}.vcf` — **PK-based; this is what breaks cross-DB reuse.** |
| Range lock | `models.py:952` `AnnotationRangeLock` | `min_variant`/`max_variant` FKs + `count`. Subdivided on crash (`annotation_scheduler_task.py`). |
| Scheduler | `annotation/tasks/annotation_scheduler_task.py:16` | Creates range locks + `AnnotationRun`s, calls `annotate_variants.apply_async`. |
| Status enum | `annotation/models/models_enums.py:30` `AnnotationStatus` | `CREATED, DELETING, DUMP_STARTED, DUMP_COMPLETED, ANNOTATION_STARTED, ANNOTATION_COMPLETED, UPLOAD_STARTED, FINISHED, ERROR`. |
| Dump dir setting | `variantgrid/settings/components/annotation_settings.py:207` | `ANNOTATION_VCF_DUMP_DIR`. |
| Site name | `default_settings.py:424` | `SITE_NAME` (used as the "server" token). |

---

## 3. The matching problem (clone-only)

Reuse is supported **only between a database and its own clone** (§1 scope). A clone starts with
identical PKs, so for runs the operator chooses (those predating any divergence) the origin
`Variant.id`s line up 1:1 with the target DB. Matching therefore has two jobs:

1. **Identity** — the annotated file must be for the *same annotation version*: genome build, VEP
   version, `vep_cache`, `columns_version`, `annotation_consortium`, `gene_annotation_release` (plus
   the data-source versions on `VariantAnnotationVersion`). Largely free: the annotated VCF's `##VEP=`
   header encodes these and `vep_check_annotated_file_version_match()` already enforces it; the sidecar
   meta records them too.
2. **Run/range** — match the file to a local `AnnotationRun` of the right `pipeline_type` in
   `EXTERNAL_DUMP_COMPLETED`. On a clone the meta's `annotation_run_pk` is the primary handle; the
   variant-coordinate range of `min_variant`/`max_variant` (string repr, build-relative) is the
   cross-check.

Neither of those *proves* the variant ids line up — divergence (new variants imported, range-lock
subdivision after a crash; `annotation_scheduler_task.py`) can break alignment without changing the
version. That guarantee is the job of the **§6a per-record variant-id-alignment check**, which is the
real safety net and is always run before import.

---

## 4. Design

### 4.1 Self-describing filenames + sidecar metadata

`get_dump_filename()` becomes self-describing. Proposed stem:

```
{site}__{build}__{consortium}__vep{vep}__cv{columns_version}__gar{gene_annotation_release}__run{pk}__{type}.vcf
```

- `site` = slugified `settings.SITE_NAME` (the "server" token).
- `run{pk}` stays for human/debug uniqueness on the *origin* server; it is **not** used for
  cross-DB matching.
- The annotated counterpart keeps the existing `.vep_annotated_{build}.vcf.gz` suffix.

Filenames are for humans and globbing; **matching is driven by a sidecar JSON**, written next to
the dump at dump time (`dump_{stem}.meta.json`):

```json
{
  "schema": 1,
  "site_name": "sapath_test",
  "annotation_run_pk": 1234,
  "pipeline_type": "S",
  "genome_build": "GRCh38",
  "annotation_consortium": "R",
  "variant_annotation_version": {
    "pk": 42, "vep": 110, "vep_cache": "...", "columns_version": 3,
    "gene_annotation_release": 17, "distance": 5000,
    "gnomad": "...", "dbnsfp": "...", "cosmic": "...", "...": "..."
  },
  "range": {
    "min_variant": "1:12345 G>A",
    "max_variant": "22:998877 T>C",
    "count": 24000
  },
  "dump_count": 24000
}
```

The `variant_annotation_version` block is the existing version-identity fields (reuse the same
kwargs that `vep_dict_to_variant_annotation_version_kwargs()` builds). The `range` strings are the
`Variant` string repr — this is the issue's requested "variant string representation … as a sanity
check."

Add helpers on `AnnotationRun`: `get_dump_metadata() -> dict`, `write_dump_metadata()`, and a
classmethod `parse_dump_metadata(path)`.

### 4.2 The `annotation_external` management command (single entry point)

One command drives everything, selected by params. No global settings flag; external annotation is an
explicit operator action on a chosen build/version.

```
manage.py annotation_external --dump   --genome-build GRCh38 [--pipeline-type S] --output-dir DIR
manage.py annotation_external --import --genome-build GRCh38 --input-dir DIR [--dry-run]
```

**`--dump` mode:**
- Selects the `AnnotationRun`s to externally annotate for the given build/version/pipeline_type
  (those not yet finished; create range locks + runs the same way the scheduler does if they don't
  exist yet — reuse `annotation_scheduler_task` helpers).
- For each run: write the dump VCF (reuse `_unannotated_variants_to_vcf()` / `write_qs_to_vcf()`),
  write the sidecar `.meta.json` (§4.1), set `external=True` and **stop before VEP** → status
  `EXTERNAL_DUMP_COMPLETED`.
- Also emit the Snakemake bundle into `--output-dir` (§4.3).

**`--import` mode:** §4.4.

**New per-run state.** Add a boolean marker field, not a global mode:
- `AnnotationRun.external = models.BooleanField(default=False)`, set at dump time by the command.
- New `AnnotationStatus` member `EXTERNAL_DUMP_COMPLETED = 'e', "Awaiting external annotation"`;
  surface it in `get_summary_state()` as "Awaiting external annotation".
- `get_status()` gains a branch: if `external` and `dump_end` is set and `vcf_annotated_filename is
  None` → `EXTERNAL_DUMP_COMPLETED` (rather than getting stuck at `DUMP_COMPLETED`). It is **not** a
  completed state — it is genuinely waiting on the operator.
- The normal scheduler/`annotate_variants` must **skip `external=True` runs** so it never auto-runs
  VEP on a run the operator is managing externally (guard in `annotation_scheduler_task` and at the
  top of `annotate_variants`).
- Migration: one new field (`external`) + the choices update (no-op DB change for the enum value).

`annotate_variants()` already guards import on `vcf_annotated_filename` being set, so once `--import`
fills `vcf_annotated_filename` the normal import path runs unchanged.

### 4.3 Snakemake generation (part of `--dump`)

The `--dump` run also writes, into `--output-dir` (the operator copies this whole dir to the compute box):
- the dump VCFs + sidecar `.meta.json`,
- a `config.yaml` with VEP paths (`ANNOTATION_VEP_*` settings) **as overridable placeholders**, plus
  the dump dir, output dir, fork count, buffer sizes,
- a `Snakefile` whose rule body is the real VEP command, built from `get_vep_command()` with
  `{input}`/`{output}` substituted and every server path read from `config.yaml` (so the compute box
  can have different VEP install paths than the VM). One rule per `pipeline_type`; SV runs also need
  the AnnotSV step (`run_annotsv`, `annotsv_annotation.py:36`) — second rule, or VEP-only first (§7).

The Snakefile discovers work by globbing `*.meta.json` so it is self-contained on the compute box.
Reusing `get_vep_command()` verbatim keeps the external run byte-identical to the in-VM run (and the
`##VEP=` header check passes on import).

### 4.4 `--import` mode (re-import annotated VCFs)

Inputs: `--input-dir` of annotated VCFs + their sidecar metas.

For each annotated VCF:

1. Read sidecar meta → version identity + coordinate range; read `##VEP=` header and confirm it matches
   the meta (defence in depth).
2. **Match to a local `AnnotationRun`** (this DB) — clone, so ids should align:
   - Find the local `VariantAnnotationVersion` whose identity equals the meta's
     (build/consortium/vep/columns_version/gene_annotation_release/data versions).
   - Find the local `AnnotationRun`(s) of the right `pipeline_type`, in `EXTERNAL_DUMP_COMPLETED`,
     whose range matches the file (the meta's `annotation_run_pk` is the primary key on a clone; the
     coordinate range is the cross-check).
3. **Run the §6a variant-id-alignment check (always on).** Abort the whole import on any mismatch.
4. For each verified run: set `vcf_annotated_filename` (copy/symlink into `ANNOTATION_VCF_DUMP_DIR`) and
   call the existing **upload-only** path (`annotation_run_retry(upload_only=True)` /
   `import_vcf_annotations()`), which already runs `vep_check_annotated_file_version_match()` and sets
   `upload_*` timestamps → status `FINISHED`.

Reuses `import_vcf_annotations()` unchanged; the §6a pre-flight (§6a) is what makes trusting the
origin-DB `variant_id` safe.

`--dry-run` lists matches and any unmatched files/runs without importing.

---

## 5. Reuse flow end-to-end (the sapath scenario)

1. **test VM:** `annotation_external --dump --genome-build GRCh38 --output-dir DIR` → every selected run
   is `EXTERNAL_DUMP_COMPLETED` with VCF + meta on disk, plus the Snakemake bundle in `DIR`.
2. Operator copies `DIR` to the compute box.
3. Compute box: edit `config.yaml` VEP paths, `snakemake` → annotated VCFs.
4. Copy annotated VCFs (+ metas) back to test VM; `annotation_external --import --input-dir DIR` →
   §6a check passes → test runs `FINISHED`.
5. **Keep the annotated VCFs + metas.**
6. **prod** (clone needing the same annotation version): the operator copies back only the annotated
   files whose id range predates any divergence and runs `annotation_external --import` on prod. The
   §6a check guarantees every imported run's ids line up; anything that doesn't is rejected and falls
   back to a normal VEP run.

---

## 6. Files to add / change (summary)

**New:**
- `annotation/management/commands/annotation_external.py` — single command, `--dump` / `--import` modes
  (dump writes VCFs + meta + Snakemake bundle; import matches + verifies + upload-only).
- `annotation/external_annotation.py` — shared helpers: metadata schema (build/parse/validate),
  run selection/matching, `verify_annotated_vcf_variant_ids()` (§6a), Snakemake/config emission.
- Snakefile/config templates (under `annotation/templates/` or generated inline).
- Tests: `annotation/tests/test_external_annotation.py` (meta round-trip, §6a pass/fail cases,
  import-only happy path on fixture VCFs).

**Changed:**
- `annotation/models/models_enums.py` — add `EXTERNAL_DUMP_COMPLETED`; update summary map.
- `annotation/models/models.py` — `AnnotationRun.external` field; `get_status()` branch;
  `get_dump_filename()` self-describing stem; `get_dump_metadata()/write_dump_metadata()/parse`.
- `annotation/tasks/annotate_variants.py` — skip `external=True` runs in `annotate_variants`; factor
  the dump step so `--dump` can call it and stop before VEP.
- `annotation/tasks/annotation_scheduler_task.py` — skip `external=True` runs.
- One migration (the `external` field + choices state).

---

## 6a. Variant-ID alignment safety check (the hard fail)

**Operator strategy (decided):** the `variant_id` baked into a dump's INFO is the *origin* DB's
`Variant.id`. Rather than re-resolving variants by coordinate, the operator will **manually copy
back only the annotated files whose ID range predates the test/prod divergence (the "split")**, where
origin PKs still line up 1:1 with the target DB. The job of VariantGrid is to **fail loudly and
automatically if the operator gets this wrong** — never silently import misaligned annotation.

This is cheap and bulletproof because every annotated VCF record already carries *both* identifiers:

- `variant_id = v.INFO["variant_id"]` — the origin DB's `Variant.id`
  (`bulk_vep_vcf_annotation_inserter.py:715`)
- the true coordinate `CHROM/POS/REF/ALT(/SVLEN)` of that record (`:708`)

and a local `Variant` exposes its coordinate via `Variant.coordinate`
(`snpdb/models/models_variant.py:660`).

**Check:** for every record, the **local** `Variant` whose `pk == variant_id` must have a coordinate
equal to the record's coordinate. On a pre-split clone all PKs align → every row matches → pass. If a
post-split file (or a file from the wrong DB) is supplied, at least one PK maps to a different (or
missing) local variant → **fail before any rows are inserted**.

**Implementation — `verify_annotated_vcf_variant_ids(annotation_run)`** (new, in
`annotation/external_annotation.py`), run as a **pre-flight in the import command before
`import_vcf_annotations()`**:

1. Stream the annotated VCF with cyvcf2, reading only `INFO["variant_id"]` and raw
   `CHROM/POS/REF/ALT/SVLEN` (skip CSQ/HGVS work — fast).
2. In batches (`SQL_BATCH_INSERT_SIZE`): `Variant.objects.filter(pk__in=ids).values("id",
   "locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq", "svlen")` → compare each
   to the file's raw coordinate.
3. Fail on the **first** mismatch or missing id with an actionable message, e.g.:
   `"AnnotationRun N: variant_id 12345 in <file> is 1:999 G>A but local variant 12345 is 7:42 T>C —
   this annotated VCF was produced against a different/diverged database (id past the split?). Aborting,
   no data imported."`
4. Count must also reconcile (rows verified == expected non-skipped count) so a truncated file is caught.

Notes:
- Compare on **raw** stored fields (pre-`as_external_explicit`) so the check matches exactly how the
  dump wrote them; prove exact round-trip on a same-DB dump in tests first, including the reference-alt
  special value and symbolic/SV representation.
- This stacks with the existing `vep_check_annotated_file_version_match()` (which guards *version*
  identity from the `##VEP=` header) — together they reject both "wrong version" and "wrong/diverged DB".
- Pre-flight (separate read pass) is preferred over inline-during-insert so we fail **before** writing
  any partial data; it's an extra cheap read (no HGVS/insert work) vs. a per-batch PK-IN query.
- Always on for the external-import command (not opt-in). Offer `--allow-id-mismatch` only as an
  explicit, loudly-logged escape hatch if ever needed; default is hard fail.

---

## 7. Open questions / risks (resolve before coding)

1. **Range non-alignment under subdivision.** A test file's `[min,max]` may partially overlap a
   prod run's range. Decide: only reuse on full containment (skip partial → normal VEP), or support
   splitting an annotated file across runs. v1: full-containment only, log skips loudly
   ("no silent caps"). The §6a per-row id-alignment check makes this safe regardless — a partially
   overlapping file imported by mistake fails on the first divergent id.
2. **SV pipeline / AnnotSV** in the external workflow — include in v1 or document as VEP-only first?
3. **Security/trust** — importing externally-produced VCFs runs the existing version checks plus the
   §6a id-alignment check; confirm that is sufficient (validates VEP/data versions + variant identity,
   not human provenance).
4. **Cleanup** — kept annotated VCFs accumulate; add retention guidance / a prune command later.

---

## 8. Suggested build order

1. Metadata schema + self-describing filenames + `get_dump_metadata()` (+ tests). *No behaviour change.*
2. `EXTERNAL_DUMP_COMPLETED` state + `AnnotationRun.external` field + scheduler/`annotate_variants` skip
   external runs (+ migration).
3. `annotation_external --dump`: select runs, write VCF + meta, park in `EXTERNAL_DUMP_COMPLETED`.
4. **`verify_annotated_vcf_variant_ids()` id-alignment check (§6a) + tests** — the safety net the whole
   manual-copy workflow relies on. Prove it passes on a same-DB dump and fails on a deliberately
   misaligned one.
5. `annotation_external --import`: same-DB round-trip — §6a pre-flight → existing upload-only path.
6. Snakemake bundle emission inside `--dump`.
7. Clone-reuse (prod) scenario end-to-end, with the §6a check guarding every import.
