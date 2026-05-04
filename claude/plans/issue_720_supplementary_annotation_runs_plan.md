# Issue #720 — Split annotation tools into independent `AnnotationRun`s

Generalises the existing `AnnotationRun` so that VEP, AnnotSV, and SpliceAI (and future tools) each get their own run per range lock, with explicit dependencies. Lets us:

- Run SpliceAI as a CLI tool to fill gaps the precalculated dbNSFP/SpliceAI cache misses (the original ask in #720).
- Roll AnnotSV independently of VEP (no full re-annotation when the AnnotSV bundle bumps).
- Isolate failures: a failed SpliceAI/AnnotSV run no longer threatens VEP data.

There is currently **no production AnnotSV deployment** (dev machine only), so the AnnotSV migration is greenfield — we don't need backwards-compat shims for `annotsv_*` columns on `AnnotationRun` / `VariantAnnotationVersion`.

---

## Key model changes

### 1. `AnnotationRun` gets `variant_type`, `tool`, `depends_on`

`pipeline_type` is renamed to `variant_type` (it describes the variant scope, not a pipeline). The pipeline a run executes is now its `tool`.

```python
# annotation/models/models_enums.py

class VariantType(models.TextChoices):
    """ Which kind of variant this AnnotationRun is scoped to. """
    SHORT = "S", "Short Variant"
    STRUCTURAL_VARIANT = "C", "Structural Variant"


class AnnotationTool(models.TextChoices):
    """ Which annotator runs in this AnnotationRun. """
    VEP = "V", "VEP"
    SPLICEAI = "P", "SpliceAI"
    ANNOTSV = "A", "AnnotSV"
```

`VariantAnnotationPipelineType` → delete. The choice values `"S"` and `"C"` are preserved on the new `VariantType` so DB data round-trips without rewriting rows.

```python
# annotation/models/models.py — AnnotationRun

class AnnotationRun(TimeStampedModel):
    status = models.CharField(max_length=1, choices=AnnotationStatus.choices, default=AnnotationStatus.CREATED)
    annotation_range_lock = models.ForeignKey(AnnotationRangeLock, null=True, on_delete=CASCADE)
    variant_type = models.CharField(max_length=1, choices=VariantType.choices, default=VariantType.SHORT)
    tool = models.CharField(max_length=1, choices=AnnotationTool.choices, default=AnnotationTool.VEP)
    # Dependency: this run won't start until depends_on is FINISHED.
    depends_on = models.ForeignKey('self', null=True, blank=True, on_delete=CASCADE,
                                   related_name='dependents')
    # ... existing fields (status, task_id, dump_*, annotation_*, upload_*, pipeline_*) ...

    class Meta:
        unique_together = ('annotation_range_lock', 'variant_type', 'tool')
```

Fields to **delete** from `AnnotationRun`:
- `annotsv_tsv_filename`
- `annotsv_error`
- `annotsv_imported`

Provenance for AnnotSV moves into the AnnotSV run's existing fields (`pipeline_command`, `pipeline_stdout`, `pipeline_stderr`, `error_exception`, `vcf_annotated_filename` — repurposed to hold the TSV path), which is symmetric with how VEP runs work.

### 2. `VariantAnnotationVersion` loses `annotsv_*`

```python
# annotation/models/models.py — VariantAnnotationVersion
# DELETE these fields:
#   annotsv_code = models.TextField(null=True, blank=True)
#   annotsv_bundle = models.TextField(null=True, blank=True)
```

AnnotSV's version is no longer pinned on the VEP version. It moves to `SupplementaryAnnotationVersion` (below).

### 3. New model `SupplementaryAnnotationVersion`

```python
# annotation/models/models.py

class SupplementaryAnnotationVersion(TimeStampedModel):
    """ Pins the version of a non-VEP annotation tool. Independent of
        VariantAnnotationVersion so a tool can roll without forcing a VEP rerun.

        AnnotationRuns for tools other than VEP carry an FK to one of these. """
    tool = models.CharField(max_length=1, choices=AnnotationTool.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    code_version = models.TextField()        # e.g. "1.3.1" for SpliceAI, "3.5.8" for AnnotSV
    data_version = models.TextField(null=True, blank=True)  # AnnotSV bundle, SpliceAI ref-fasta md5
    params_hash = models.CharField(max_length=64, null=True, blank=True)  # hash of CLI params, optional
    status = models.CharField(max_length=10, choices=VariantAnnotationVersion.Status.choices,
                              default=VariantAnnotationVersion.Status.NEW)

    class Meta:
        unique_together = ('tool', 'genome_build', 'code_version', 'data_version', 'params_hash')
        constraints = [
            models.UniqueConstraint(
                fields=["tool", "genome_build"],
                condition=Q(status="ACTIVE"),
                name="one_active_supp_version_per_tool_per_build",
            ),
        ]

    @staticmethod
    def latest_active(tool, genome_build) -> Optional['SupplementaryAnnotationVersion']:
        return SupplementaryAnnotationVersion.objects.filter(
            tool=tool, genome_build=genome_build, status='ACTIVE'
        ).order_by('-created').first()
```

`AnnotationRun` gets a nullable FK:

```python
supplementary_version = models.ForeignKey(SupplementaryAnnotationVersion,
                                          null=True, blank=True, on_delete=PROTECT)
```

Required when `tool != VEP`, NULL when `tool == VEP`. Enforce in `clean()` / DB constraint.

---

## Tool runner interface

Replace the inline `if pipeline_type == STANDARD ... else ...` branches in `dump_and_annotate_variants` with a registry.

```python
# annotation/runners/base.py  (new)

class AnnotationToolRunner(abc.ABC):
    tool: AnnotationTool
    variant_type: VariantType  # which variant scope this runner handles

    @abc.abstractmethod
    def select_variants_q(self, annotation_version) -> Q:
        """ Extra Q filter on top of variants-in-range. e.g. SpliceAI: indels
            with missing spliceai_max_ds. AnnotSV: no extra filter (all SVs).
            VEP: variantannotation__isnull=True (current behaviour). """

    @abc.abstractmethod
    def run(self, annotation_run: AnnotationRun, vcf_in: str) -> str:
        """ Run the tool. Return path to output file (vcf or tsv). Should set
            annotation_run.pipeline_command / pipeline_stdout / pipeline_stderr. """

    @abc.abstractmethod
    def apply(self, annotation_run: AnnotationRun, output_path: str) -> int:
        """ Import results into VariantAnnotation. Returns rows affected.
            VEP runner: existing BulkVEPVCFAnnotationInserter (insert).
            SpliceAI/AnnotSV: bulk UPDATE on existing VariantAnnotation rows. """


# annotation/runners/__init__.py
RUNNERS: dict[tuple[VariantType, AnnotationTool], AnnotationToolRunner] = {
    (VariantType.SHORT, AnnotationTool.VEP): VEPShortRunner(),
    (VariantType.STRUCTURAL_VARIANT, AnnotationTool.VEP): VEPSVRunner(),
    (VariantType.SHORT, AnnotationTool.SPLICEAI): SpliceAIRunner(),
    (VariantType.STRUCTURAL_VARIANT, AnnotationTool.ANNOTSV): AnnotSVRunner(),
}
```

`annotate_variants` task body becomes:

```python
@celery.shared_task
def annotate_variants(annotation_run_id):
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)

    if annotation_run.depends_on_id:
        dep = annotation_run.depends_on
        if dep.status != AnnotationStatus.FINISHED:
            # Re-queue with countdown; scheduler will also pick it up.
            annotate_variants.apply_async((annotation_run_id,), countdown=60)
            return

    # ... existing task_id lock acquisition ...

    runner = RUNNERS[(annotation_run.variant_type, annotation_run.tool)]
    # dump_and_run_and_apply replaces dump_and_annotate_variants + import_vcf_annotations
    if annotation_run.vcf_annotated_filename is None:
        vcf_in = _dump_variants_for_runner(annotation_run, runner)
        if annotation_run.dump_count:
            output = runner.run(annotation_run, vcf_in)
            annotation_run.vcf_annotated_filename = output
            annotation_run.annotation_end = timezone.now()
            annotation_run.save()
    if annotation_run.vcf_annotated_filename:
        runner.apply(annotation_run, annotation_run.vcf_annotated_filename)
```

`_dump_variants_for_runner` calls `get_variants_qs_for_annotation` and ANDs `runner.select_variants_q(...)`. For VEP runners this matches today's behaviour; for SpliceAI/AnnotSV this naturally selects only the gaps because `depends_on` already ensured VEP populated `VariantAnnotation` first.

---

## Per-runner notes

### VEP runners (`VEPShortRunner`, `VEPSVRunner`)

Wrap the existing `get_vep_command` + `run_vep` + `BulkVEPVCFAnnotationInserter` flow. `select_variants_q` returns `Q(variantannotation__isnull=True)` (current behaviour from `get_variants_qs_for_annotation(annotated=False)`).

Strip the AnnotSV block from `dump_and_annotate_variants` (lines 154-175) — that becomes its own runner.

### `SpliceAIRunner`

```python
class SpliceAIRunner(AnnotationToolRunner):
    tool = AnnotationTool.SPLICEAI
    variant_type = VariantType.SHORT

    def select_variants_q(self, annotation_version) -> Q:
        # Per the issue: variants the precalculated cache missed.
        # cadd_raw_rankscore is NULL when dbNSFP wasn't applied; SpliceAI snv/indel
        # cache covers the same scope as dbNSFP for our purposes.
        # Optionally also bound by max indel length (settings.ANNOTATION_SPLICEAI_MAX_INDEL_LEN)
        # since SpliceAI silently skips above 2 * --distance (default 50bp).
        return Q(variantannotation__cadd_raw_rankscore__isnull=True) & \
               Q(variantannotation__spliceai_max_ds__isnull=True)

    def run(self, annotation_run, vcf_in) -> str:
        sv = annotation_run.supplementary_version
        vcf_out = vcf_in.replace(".vcf", ".spliceai.vcf")
        cmd = [
            settings.ANNOTATION_SPLICEAI_BIN,  # default "spliceai"
            "-I", vcf_in, "-O", vcf_out,
            "-R", VEPConfig(annotation_run.genome_build)["fasta"],
            "-A", "grch37" if annotation_run.genome_build.is_grch37 else "grch38",
        ]
        if settings.ANNOTATION_SPLICEAI_DISTANCE:
            cmd.extend(["-D", str(settings.ANNOTATION_SPLICEAI_DISTANCE)])
        annotation_run.pipeline_command = " ".join(cmd)
        annotation_run.annotation_start = timezone.now()
        annotation_run.save()
        rc, stdout, stderr = execute_cmd(cmd)
        annotation_run.pipeline_stdout = (stdout or "")[:1_000_000]
        annotation_run.pipeline_stderr = (stderr or "")[:1_000_000]
        if rc != 0:
            raise RuntimeError(f"SpliceAI returned {rc}")
        return vcf_out

    def apply(self, annotation_run, output_path) -> int:
        # Parse SpliceAI INFO field (DS_AG/DS_AL/DS_DG/DS_DL/DP_AG/...) and
        # bulk-UPDATE the existing VariantAnnotation row for this variant in
        # this version partition. Update spliceai_max_ds + spliceai_gene_symbol
        # at the same time. Use VariantAnnotation.objects.bulk_update for batches.
        ...
```

Selection-Q rationale: `cadd_raw_rankscore IS NULL` is the canary the issue identified — every dbNSFP record has it, so its absence means dbNSFP (and the SpliceAI snv/indel cache) didn't cover this variant. We additionally require the existing SpliceAI fields are NULL so we don't overwrite cache-derived scores.

### `AnnotSVRunner`

Wraps `run_annotsv` + `import_annotsv_tsv` (existing functions). `select_variants_q` is `Q()` (no extra filter — every SV gets AnnotSV). The runner sets `annotation_run.vcf_annotated_filename` to the TSV path; `apply` calls `import_annotsv_tsv(annotation_run)` (signature already takes an `AnnotationRun`).

Delete:
- `annotation/management/commands/annotsv_run.py` lines that read `annotsv_tsv_filename`/`annotsv_error` — replace with reading the AnnotSV `AnnotationRun`.
- `annotsv_check_command_line_version_match` callsite in `dump_and_annotate_variants`. Move check to `AnnotSVRunner.run()` start, comparing `annotation_run.supplementary_version.code_version` against `annotsv -version`.
- The AnnotSV try/except block in `dump_and_annotate_variants`.

---

## Scheduler changes

`_handle_range_lock` creates one run per `(variant_type, tool)` enabled in settings, sets `depends_on`, and queues the VEP runs immediately (the dependents will be queued by VEP completion or picked up by the scheduler's next pass).

```python
# annotation/tasks/annotation_scheduler_task.py

def _create_runs_for_range_lock(range_lock):
    """ Creates VEP run + dependent supplementary runs in one transaction. """
    genome_build = range_lock.version.genome_build
    variant_type = range_lock.variant_type  # see Q below — need to add this to AnnotationRangeLock OR derive
    with transaction.atomic():
        vep_run, vep_created = AnnotationRun.objects.get_or_create(
            annotation_range_lock=range_lock,
            variant_type=variant_type,
            tool=AnnotationTool.VEP,
        )
        for tool in settings.ANNOTATION_TOOLS_ENABLED.get(variant_type, []):
            if tool == AnnotationTool.VEP:
                continue
            sv = SupplementaryAnnotationVersion.latest_active(tool, genome_build)
            if sv is None:
                continue  # tool enabled but no active version => skip
            AnnotationRun.objects.get_or_create(
                annotation_range_lock=range_lock,
                variant_type=variant_type,
                tool=tool,
                defaults={'depends_on': vep_run, 'supplementary_version': sv},
            )
    if vep_created:
        annotate_variants.apply_async((vep_run.pk,))


def _runnable_qs(variant_annotation_version):
    """ AnnotationRuns whose dependency is satisfied (or absent) and which haven't
        finished. The scheduler can pick these up between fresh range locks. """
    return AnnotationRun.objects.filter(
        annotation_range_lock__version=variant_annotation_version,
    ).exclude(status__in=AnnotationStatus.get_completed_states()).filter(
        Q(depends_on__isnull=True) | Q(depends_on__status=AnnotationStatus.FINISHED)
    ).filter(task_id__isnull=True)
```

The scheduler additionally walks `_runnable_qs` and queues anything pending whose dependency just finished — this is the replacement for `annotation_run_complete_signal` driving fan-out. The signal stays for what it's actually for (cache invalidation downstream), it just no longer creates new runs.

A note on `AnnotationRangeLock`: today range locks aren't typed by variant_type — there's one set of locks and each lock gets a STANDARD run + a STRUCTURAL_VARIANT run. That's preserved; `_create_runs_for_range_lock` is called once per (range_lock, variant_type) pair, mirroring today's `for pipeline_type in VariantAnnotationPipelineType` loop.

---

## Settings

```python
# variantgrid/settings/components/annotation_settings.py

# Per-variant-type list of tools that should be scheduled. VEP is implicit.
ANNOTATION_TOOLS_ENABLED = {
    "S": [AnnotationTool.VEP],                      # default: VEP only on short variants
    "C": [AnnotationTool.VEP],                      # default: VEP only on SVs
}
# Deployments opt in by overriding, e.g.:
#   ANNOTATION_TOOLS_ENABLED["S"].append(AnnotationTool.SPLICEAI)
#   ANNOTATION_TOOLS_ENABLED["C"].append(AnnotationTool.ANNOTSV)

ANNOTATION_SPLICEAI_BIN = "spliceai"
ANNOTATION_SPLICEAI_DISTANCE = 50           # SpliceAI -D
ANNOTATION_SPLICEAI_MAX_INDEL_LEN = None    # None = no cap
```

Delete: `ANNOTATION_ANNOTSV_ENABLED` (replaced by `ANNOTATION_TOOLS_ENABLED["C"]` membership). Keep all the other `ANNOTATION_ANNOTSV_*` settings — the AnnotSV runner still uses them.

A new Celery queue for tools that may run long / on GPU:

```python
# variantgrid/settings/components/celery_settings.py
# Add 'supplementary_annotation_workers' to CELERY_TASK_QUEUES
```

`annotate_variants` routes to the right queue based on `tool`:

```python
queue = 'supplementary_annotation_workers' if annotation_run.tool != AnnotationTool.VEP else 'annotation_workers'
annotate_variants.apply_async((annotation_run.pk,), queue=queue)
```

---

## Migration

Single migration — no data backfill needed beyond column adds because:
- Existing rows have `pipeline_type` `"S"`/`"C"` → renamed to `variant_type`, same DB values.
- `tool` defaults to `AnnotationTool.VEP` on existing rows (current rows are all VEP).
- `annotsv_*` columns on `AnnotationRun` and `VariantAnnotationVersion` are dropped (no production AnnotSV data).
- `unique_together` updates from `(annotation_range_lock, pipeline_type)` to `(annotation_range_lock, variant_type, tool)` — non-conflicting because all existing rows have `tool=VEP`.

```python
# annotation/migrations/0XYZ_split_annotation_tools.py
operations = [
    migrations.RenameField('AnnotationRun', 'pipeline_type', 'variant_type'),
    migrations.AddField('AnnotationRun', 'tool', models.CharField(max_length=1,
        choices=AnnotationTool.choices, default=AnnotationTool.VEP)),
    migrations.AddField('AnnotationRun', 'depends_on',
        models.ForeignKey('self', null=True, blank=True, on_delete=CASCADE,
                          related_name='dependents')),
    migrations.AddField('AnnotationRun', 'supplementary_version',
        models.ForeignKey('SupplementaryAnnotationVersion', null=True,
                          blank=True, on_delete=PROTECT)),
    migrations.RemoveField('AnnotationRun', 'annotsv_tsv_filename'),
    migrations.RemoveField('AnnotationRun', 'annotsv_error'),
    migrations.RemoveField('AnnotationRun', 'annotsv_imported'),
    migrations.RemoveField('VariantAnnotationVersion', 'annotsv_code'),
    migrations.RemoveField('VariantAnnotationVersion', 'annotsv_bundle'),
    migrations.AlterUniqueTogether('AnnotationRun',
        unique_together={('annotation_range_lock', 'variant_type', 'tool')}),
    migrations.CreateModel('SupplementaryAnnotationVersion', fields=[...]),
]
```

`ColumnVEPField.pipeline_type` (annotation/migrations/0080) — leave as-is. That field describes which VEP runs a column applies to and is still pipeline-typed in the VEP sense; it's unrelated to `AnnotationRun.variant_type`. (See open question 1.)

---

## Callsites to rename `pipeline_type` → `variant_type`

Single sweep replacing:
- `AnnotationRun.pipeline_type` → `AnnotationRun.variant_type`
- `get_variants_qs_for_annotation(pipeline_type=…)` → `variant_type=…`
- `VariantAnnotationPipelineType` → `VariantType`
- `settings.ANNOTATION_VEP_BUFFER_SIZE` keys (still keyed by `"S"`/`"C"`, just rename the comment/lookup variable)

Files (from `grep -rn "pipeline_type\|VariantAnnotationPipelineType"`):
- annotation/vep_columns.py
- annotation/views.py
- annotation/admin.py
- annotation/grids.py
- annotation/annotation_version_querysets.py
- annotation/vep_annotation.py
- annotation/models/models_enums.py
- annotation/models/models.py
- annotation/vcf_files/import_vcf_annotations.py
- annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py
- annotation/tasks/annotate_variants.py
- annotation/tasks/annotation_scheduler_task.py
- annotation/management/commands/vep_run.py
- annotation/management/commands/annotsv_run.py
- annotation/management/commands/fix_annotation_sv_c_hgvs.py
- variantgrid/deployment_validation/vep_columns_check.py
- annotation/tests/test_vep_columns.py

Plus templates / admin display strings that show the field label.

`ColumnVEPField.pipeline_type` is a separate column on a different model — leave its name alone (it's the variant scope a VEP column applies to, conceptually it should also become `variant_type`, but a wider rename — see open question 1).

---

## Tests to add

- `test_annotation_run_dependency` — creating an SV+VEP run + SV+ANNOTSV run wires `depends_on` correctly; the dependent task no-ops while parent is unfinished.
- `test_runner_registry` — every `(variant_type, tool)` combo enabled in settings has a registered runner; missing combos raise at scheduler time.
- `test_spliceai_select_q` — selection Q targets indels missing both `cadd_raw_rankscore` and `spliceai_max_ds`; ignores SNVs/variants the cache covered.
- `test_spliceai_apply_updates_existing_row` — running SpliceAI updates the existing `VariantAnnotation` row's `spliceai_*` fields without disturbing other columns.
- `test_annotsv_failure_does_not_fail_vep_run` — AnnotSV runner failure is recorded on the AnnotSV run only; SV+VEP run remains FINISHED.
- `test_supplementary_version_roll` — promoting a new `SupplementaryAnnotationVersion` to ACTIVE causes the next scheduler pass to create new SpliceAI/AnnotSV runs against the new version (existing FINISHED runs against the old version are untouched).

---

## Implementation order

1. Add `VariantType` + `AnnotationTool` enums; rename `pipeline_type` → `variant_type` across the codebase. Migration #1. Tests should still pass.
2. Add `tool`, `depends_on`, `supplementary_version` fields + `SupplementaryAnnotationVersion` model. Migration #2. Existing runs default to `tool=VEP`, no behavioural change.
3. Introduce `AnnotationToolRunner` interface; refactor existing VEP path into `VEPShortRunner` / `VEPSVRunner`. Tests should still pass — pure refactor.
4. Strip AnnotSV out of `dump_and_annotate_variants`; introduce `AnnotSVRunner`; remove `annotsv_*` fields from `AnnotationRun` and `VariantAnnotationVersion`. Migration #3. Update `annotsv_run` management command to operate on the AnnotSV `AnnotationRun`.
5. Add scheduler dependency handling; add `_runnable_qs` walk.
6. Add `SpliceAIRunner` + selection Q + parser/apply.
7. Add `supplementary_annotation_workers` Celery queue.
8. Tests as above.

---

## Open questions for review

1. **`ColumnVEPField.pipeline_type`** — same name, different model, conceptually the same axis (which variant scope this VEP column applies to). Rename it too in step 1, or leave alone? Renaming touches more migrations and the `vep_columns.py` filtering logic.

2. **SpliceAI per-variant provenance table** — earlier I proposed `SupplementaryAnnotationApplied(variant, supp_version, …)` so we can answer "has SpliceAI v1.3.1 run on variant X". With the current design the answer is implied by `VariantAnnotation.spliceai_max_ds IS NOT NULL` for a given partition, plus the AnnotationRun history for the range. Is that good enough, or do you want the explicit per-variant table? It costs ~1 row per supp-annotated variant.

3. **VAV roll behaviour** — when a new `VariantAnnotationVersion` is created, the new partition starts empty, so SpliceAI/AnnotSV must re-run for every range. That's symmetric with VEP. Worth a "copy forward supplementary scores when supp version unchanged" optimisation in v1, or defer? My lean: defer.

4. **SpliceAI `apply()` bulk strategy** — `bulk_update` on potentially millions of rows is slow. Alternative: write a small VCF and use the existing `BulkVEPVCFAnnotationInserter` path with `update_or_create` semantics. Preference?

5. **`AnnotationStatus.FINISHED` when `dump_count == 0`** — today a VEP run with zero variants in scope finishes immediately. For a SpliceAI runner, "no gap variants in this range" should also short-circuit to FINISHED. Confirm this is the intended semantics (so dependents don't wait on a no-op).

6. **Celery queue name** — `supplementary_annotation_workers` OK, or prefer something like `gpu_annotation_workers` to make the GPU intent explicit? (SpliceAI benefits from GPU; AnnotSV doesn't.)

7. **AnnotSV TSV path field** — repurpose `vcf_annotated_filename` to hold the TSV path (cheapest), or add a generic `output_filename` field and migrate? My lean: repurpose, since the field is already a generic "where did the tool write its output".

8. **Order of step 1 vs scope of this branch** — the rename touches a lot of files. Land it as a separate PR, or include in the same branch? My lean: separate PR — pure rename, easy review.
