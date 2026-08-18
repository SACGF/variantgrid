# Issue #720 — Annotation pipeline types: split AnnotSV out of the VEP SV run

**Supersedes the May 2026 version of this plan.** That plan predates #2667 (dispatcher + leases), #1646
(count lane + merge), #1649 (split VEP / import lanes), #1654 (in-place retry), #1658 (per-task paths +
subprocess abort), #1660 (derived path ownership), #1670 (single file-cleanup owner) and #1701 (truncation
checks), and predates `GENE_LEVEL` becoming a pipeline type. Its two-column `variant_type` × `tool` rename
is not the shape we want any more — see below.

---

## 1. What changed, and why the shape of the answer changed with it

The May plan's central claim was that `pipeline_type` was misnamed: it described *which variants*, not
*which pipeline*, so it should be renamed `variant_type` and a second `tool` column added.

`GENE_LEVEL` landing made that claim false. `VariantAnnotationPipelineType` now holds three values that are
each **a tool applied to a class of variant**:

| pipeline type | variant class | tool |
|---|---|---|
| `STANDARD` | short variants | VEP |
| `STRUCTURAL_VARIANT` | symbolic variants | VEP |
| `GENE_LEVEL` | gene fusions | local computation (`annotate_gene_level_run`) |

That is exactly what the word "pipeline" means. The column name is now correct, and the axis we were going
to add is already there — it is just implicit, spread across `if pipeline_type == …` branches in
`dump_and_annotate_variants`, `get_vep_command`, `import_vcf_annotations` and `_reset_run_for_redispatch`.

**So: keep the single column. Add `ANNOTSV` as a pipeline type, and make the per-type metadata explicit in
a registry.** No rename sweep, no new `AnnotationRun` axis, no change to `UploadedVCFPipelineMaxVariant`,
`pipeline_type_variant_q`'s callers, templates, admin or grids.

## 2. The problem this actually solves

From the issue (2026-08-18):

> Problem: SV contains both VEP and AnnotSV - if you enable AnnotSV later, you cannot backfill them.

Today AnnotSV is a stage bolted inside the SV VEP run (`dump_and_annotate_variants`
`annotation/tasks/annotate_variants.py:535-559`, import at
`annotation/vcf_files/import_vcf_annotations.py:80-88`). Three consequences:

1. **No backfill.** Enabling AnnotSV only affects SV runs that have not run yet. Existing FINISHED SV runs
   are never revisited, and `manage.py annotsv_run` cannot help — it needs `vcf_dump_filename`, which
   #1670 deletes on successful import.
2. **AnnotSV only ever sees small SVs.** It is fed the *VEP* dump, which is filtered by
   `settings.ANNOTATION_VEP_SV_MAX_SIZE` (`annotate_variants.py:650-654`). The large SVs are exactly the
   ones AnnotSV's ACMG ranking is for.
3. **Rolling the AnnotSV bundle forces a full VEP re-annotation**, because the version is pinned on
   `VariantAnnotationVersion.annotsv_code` / `annotsv_bundle`.

Splitting AnnotSV into its own pipeline type fixes all three *structurally*, because the scheduler already
backfills missing runs. `_handle_variant_annotation_version`
(`annotation/tasks/annotation_scheduler_task.py:215-221`) does this today:

```python
for pipeline_type in VariantAnnotationPipelineType:
    # Look for missing any AnnotationRun for this lock
    for range_lock in arl_qs.exclude(annotationrun__pipeline_type=pipeline_type):
        _handle_range_lock(range_lock, pipeline_type)
```

Add `ANNOTSV` to the enum and turn it on, and the next scheduler pass creates an `ANNOTSV` run for every
existing range lock on the ACTIVE version. An AnnotSV run dumps its own VCF from the database, so it does
not care that the VEP dump is long gone. **Backfill needs no new machinery — it is what the scheduler
already does.**

Operationally that is cheap: most range locks contain no SVs at all, so `count_annotation_run`
(`annotation_scheduler_task.py:161`) finishes them at `count == 0` on `db_workers` without ever reaching
the dispatcher.

---

## 3. Design

### 3.1 The pipeline type

```python
# annotation/models/models_enums.py

class VariantAnnotationPipelineType(models.TextChoices):
    """ An annotation pipeline is a tool applied to a class of variant. Most are VEP over a variant class;
        GENE_LEVEL and ANNOTSV run other tools. What each one selects, depends on, and runs is declared in
        annotation.pipelines - this enum is only the stored key. """
    STANDARD = "S", "Standard Short Variant"
    STRUCTURAL_VARIANT = "C", "Structural Variant"
    # Never reaches VEP - annotation is computed locally from the gene identity in the alt
    GENE_LEVEL = "G", "Gene Level"
    # Runs AnnotSV over the same variants as STRUCTURAL_VARIANT, after that run has committed its rows
    ANNOTSV = "A", "AnnotSV"
```

`pipeline_type_variant_q` gains `ANNOTSV` alongside the SV branch — same variants, different tool:

```python
# annotation/annotation_version_querysets.py

    if pipeline_type == VariantAnnotationPipelineType.STANDARD:
        return ~q_sv & ~q_gene_level
    elif pipeline_type in (VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                           VariantAnnotationPipelineType.ANNOTSV):
        return q_sv & ~q_gene_level
    elif pipeline_type == VariantAnnotationPipelineType.GENE_LEVEL:
        return q_gene_level
```

The docstring there already anticipates this: *"a type's predicate may overlap another's"*.

### 3.2 Registry: `annotation/pipelines/`

New package. The per-type facts that are currently `if` branches become data.

```python
# annotation/pipelines/base.py

@dataclass(frozen=True)
class PipelineDef:
    """ Everything the scheduler, dispatcher and upload pipeline need to know about one pipeline type,
        without knowing what its tool is. """
    pipeline_type: VariantAnnotationPipelineType
    runner: 'AnnotationPipelineRunner'
    # Must be FINISHED on the same AnnotationRangeLock before this may launch. Supplementary pipelines
    # update the VariantAnnotation rows the VEP pipeline wrote, so they queue behind it.
    depends_on: Optional[VariantAnnotationPipelineType] = None
    # Whether an unfinished run of this type holds up a VCF import (UploadedVCF.is_fully_annotated).
    # Supplementary annotation is additive - the VCF is usable without it. See #1656.
    blocks_vcf_import: bool = True
    # Settings flag gating whether the scheduler creates runs of this type at all.
    enabled_setting: Optional[str] = None

    @property
    def enabled(self) -> bool:
        return self.enabled_setting is None or bool(getattr(settings, self.enabled_setting))


class AnnotationPipelineRunner(abc.ABC):
    """ One pipeline's tool: which variants it takes, how it runs, how its output is imported.

        The two methods map to the two dispatcher lanes (#1649): annotate() executes on annotation_workers
        and leaves the run at ANNOTATION_COMPLETED with its output on disk; import_results() executes on
        db_workers and takes it to FINISHED. """

    pipeline_type: VariantAnnotationPipelineType
    # VEP pipelines only want variants with no VariantAnnotation row yet. Supplementary pipelines set this
    # False and take every variant of their class in range: they UPDATE rows a VEP pipeline wrote, and
    # depends_on is what guarantees those rows exist.
    selects_unannotated: bool = True

    def get_variants_qs(self, annotation_run) -> QuerySet[Variant]:
        """ Variants in this run's range this pipeline is responsible for. """
        range_lock = annotation_run.annotation_range_lock
        annotation_version = range_lock.version.get_any_annotation_version()
        return get_variants_qs_for_annotation(
            annotation_version,
            pipeline_type=self.pipeline_type,
            min_variant_id=range_lock.min_variant_id,
            max_variant_id=range_lock.max_variant_id,
            annotated=not self.selects_unannotated,
        )

    @abc.abstractmethod
    def annotate(self, annotation_run, lease_heartbeat=None):
        """ Dump this run's variants and run the tool over them. """

    @abc.abstractmethod
    def import_results(self, annotation_run):
        """ Load the tool's output into the DB and set upload_end. """

    @abc.abstractmethod
    def get_output_paths(self, annotation_run, dump_filename) -> list[str]:
        """ Every file one attempt writes, derived from that attempt's dump stem (#1660). """

    @abc.abstractmethod
    def tool_finished(self, annotation_run, dump_filename) -> bool:
        """ Whether the tool ran to completion AND its output is still on disk. Drives the reclaim's
            keep-or-scrub decision. """

    @abc.abstractmethod
    def record_resume_state(self, annotation_run, dump_filename):
        """ Record what the import lane needs to resume upload-only off this attempt's finished output.
            The run's filename fields are only written at the runner's final save, so a run reclaimed
            before then has usable output its row cannot name. """
```

```python
# annotation/pipelines/__init__.py

PIPELINES: dict[VariantAnnotationPipelineType, PipelineDef] = {p.pipeline_type: p for p in [
    PipelineDef(VariantAnnotationPipelineType.STANDARD,
                VEPRunner(VariantAnnotationPipelineType.STANDARD)),
    PipelineDef(VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                VEPRunner(VariantAnnotationPipelineType.STRUCTURAL_VARIANT)),
    PipelineDef(VariantAnnotationPipelineType.GENE_LEVEL, GeneLevelRunner()),
    PipelineDef(VariantAnnotationPipelineType.ANNOTSV, AnnotSVRunner(),
                depends_on=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                blocks_vcf_import=False,
                enabled_setting="ANNOTATION_ANNOTSV_ENABLED"),
]}


def get_pipeline(pipeline_type) -> PipelineDef:
    return PIPELINES[VariantAnnotationPipelineType(pipeline_type)]


def enabled_pipeline_types() -> list[VariantAnnotationPipelineType]:
    return [pt for pt, p in PIPELINES.items() if p.enabled]


def blocking_pipeline_types() -> list[VariantAnnotationPipelineType]:
    """ Types a VCF import waits on. @see UploadedVCF.is_fully_annotated """
    return [pt for pt, p in PIPELINES.items() if p.blocks_vcf_import]


def vep_pipeline_types() -> list[VariantAnnotationPipelineType]:
    """ Types that actually invoke VEP - the only ones get_vep_command means anything for. """
    return [pt for pt, p in PIPELINES.items() if isinstance(p.runner, VEPRunner)]
```

**Import layering.** The runners need the dump/path helpers that currently live in
`annotation/tasks/annotate_variants.py`, and that module needs the registry — a cycle. Break it by moving
the file-and-dump primitives down into a new `annotation/annotation_run_files.py`, which imports only
models and querysets:

- `dump_variants`, `write_qs_to_vcf`, `get_annotated_filename`, `get_annotsv_dir`,
  `remove_run_output_files`, `conservation_sidecar_filename` re-export.

This is a move worth making on its own: `annotation/signals/annotation_run_cleanup.py` — the module whose
docstring says it is the single owner of on-disk cleanup — currently reaches up into
`annotation.tasks.annotate_variants` for its path derivation.

Resulting layers, no cycles:

```
annotation/annotation_run_files.py        paths + dumping        (models, querysets)
annotation/pipelines/{base,vep,annotsv,gene_level}.py            (+ annotation_run_files)
annotation/pipelines/__init__.py          the registry
annotation/tasks/annotate_variants.py     celery tasks           (+ pipelines)
annotation/signals/annotation_run_cleanup.py                     (+ pipelines)
```

`get_run_output_paths` moves onto the runner (`get_output_paths`) — with one tool per run there is one
output set per run, and only the runner knows its shape.

### 3.3 The task bodies collapse

`annotate_variants` (`annotation/tasks/annotate_variants.py:266-273`) currently branches on pipeline type:

```python
        with AnnotationRunLeaseHeartbeat(annotation_run, my_task_id) as lease_heartbeat:
            if annotation_run.pipeline_type == VariantAnnotationPipelineType.GENE_LEVEL:
                annotate_gene_level_run(annotation_run)
            elif annotation_run.vcf_annotated_filename is None:
                dump_and_annotate_variants(annotation_run, lease_heartbeat=lease_heartbeat)
```

becomes:

```python
        with AnnotationRunLeaseHeartbeat(annotation_run, my_task_id) as lease_heartbeat:
            if annotation_run.vcf_annotated_filename is None:
                get_pipeline(annotation_run.pipeline_type).runner.annotate(
                    annotation_run, lease_heartbeat=lease_heartbeat)
```

and `import_annotation_run` (`:351-352`):

```python
        with AnnotationRunLeaseHeartbeat(annotation_run, my_task_id):
            get_pipeline(annotation_run.pipeline_type).runner.import_results(annotation_run)
```

Everything else in both tasks — the `task_id` lock, `save_if_owner`, the reclaim classification, the
`finally` lease release, `_trigger_dispatch` — is untouched. Supplementary runs inherit lease renewal,
reclaim, attempt caps, per-task paths and abortability for free, which is the main reason to make them
first-class runs rather than another bolted-on stage.

### 3.4 Output filename

`vcf_annotated_filename` is reused as "this run's tool output", holding AnnotSV's TSV. That is deliberate,
not laziness: the field is load-bearing in three places that we *want* AnnotSV to inherit —

- `_dispatchable_runs_qs` resume lane (`annotation_scheduler_task.py:393-401`)
- `AnnotationRun.is_upload_resumable` (`models.py:1284-1297`)
- `_reset_run_for_redispatch`'s keep-the-expensive-output decision (`annotation_scheduler_task.py:546-561`)

so a reclaimed AnnotSV run resumes at the import step instead of re-running AnnotSV. `annotsv_tsv_filename`
is deleted.

A `RenameField` to `output_filename` is a reasonable tidy-up (the field now holds TSVs) but is orthogonal
and mechanical — do it after, or not at all.

### 3.5 Dependency gating

An AnnotSV run updates rows the SV VEP run creates, so it must not launch until that run is FINISHED. Two
places need the gate.

**Dispatcher** — `_dispatchable_runs_qs` (`annotation_scheduler_task.py:393`):

```python
def _dependency_satisfied_q() -> Q:
    """ #720: a supplementary run may only launch once the run it depends on - over the same range lock -
        is FINISHED, since it updates the VariantAnnotation rows that run wrote. Expressed as EXISTS
        subqueries rather than a join to the sibling run, so the independent-pipeline branches of the OR
        don't get filtered by an INNER JOIN they never wanted. """
    independent = [pt for pt, p in PIPELINES.items() if p.depends_on is None]
    q = Q(pipeline_type__in=independent)
    for pipeline_type, pipeline_def in PIPELINES.items():
        if pipeline_def.depends_on is None:
            continue
        dep_finished = AnnotationRun.objects.filter(
            annotation_range_lock_id=OuterRef("annotation_range_lock_id"),
            pipeline_type=pipeline_def.depends_on,
            status=AnnotationStatus.FINISHED,
        )
        q |= Q(pipeline_type=pipeline_type) & Q(Exists(dep_finished))
    return q


def _dispatchable_runs_qs(vav: VariantAnnotationVersion, now):
    created = Q(status=AnnotationStatus.CREATED)
    resume = Q(status=AnnotationStatus.ANNOTATION_COMPLETED, vcf_annotated_filename__isnull=False)
    return AnnotationRun.objects.filter(
        annotation_range_lock__version=vav,
        external=False,
        task_id__isnull=True,
    ).filter(created | resume).filter(
        Q(lease_expires__isnull=True) | Q(lease_expires__lt=now)
    ).filter(_dependency_satisfied_q())
```

**Count lane — deliberately NOT gated.** Worth stating, because the opposite looks right at first.
`count_annotation_run` (`annotation_scheduler_task.py:161`) finishes a run at `count == 0`, and it is
tempting to think an AnnotSV run would count zero before VEP has written any rows. It doesn't:
`get_variants_qs_for_annotation(annotated=True)` *drops* the "no annotation row yet" filter rather than
inverting it, so a supplementary run counts every SV in its range whether annotated or not.

That makes counting dependency-independent, which is what keeps enabling AnnotSV cheap: an SV-free range
lock — the large majority — has its AnnotSV run finished at `count == 0` on `db_workers` immediately,
without waiting on anything and without ever reaching a worker.

The count is taken through the runner so it is exactly what the dump will be, including each pipeline's
own selection (VEP's SV size cap, a supplementary pipeline's lack of one):

```python
    count = get_runner(annotation_run.pipeline_type).get_variants_qs(annotation_run).count()
```

**A blocked run is not a failure.** If the SV VEP run ERRORs, its AnnotSV run sits CREATED indefinitely. It
holds no worker slot (`_in_flight_runs_qs` treats CREATED as settled), and it blocks no VCF import
(§3.6). When the VEP run is retried and finishes, the AnnotSV run is counted and dispatched on the next
cycle. Nothing to add.

### 3.6 Keeping supplementary runs out of the VCF-import watermark

`UploadedVCF.is_fully_annotated` (`upload/models/models.py:492-505`) walks one
`UploadedVCFPipelineMaxVariant` row per pipeline type present in the VCF and asks
`get_lowest_unannotated_variant_id(pipeline_type=…)`. Rows are created by iterating the whole enum
(`upload/vcf/abstract_bulk_vcf_processor.py:86`). Left alone, adding `ANNOTSV` would make every SV-bearing
VCF wait for AnnotSV before its samples become available — precisely the failure #1656 fixed for SVs.

One change, in two places:

```python
# upload/vcf/abstract_bulk_vcf_processor.py
            for pipeline_type in blocking_pipeline_types():
                data = base_qs.filter(pipeline_type_variant_q(pipeline_type)).aggregate(m=Max("pk"))

# upload/management/commands/one_off_backfill_uploaded_vcf_pipeline_max_variant.py - same substitution
```

`annotation_run_complete_signal` already carries `pipeline_type=`, and
`upload/signals/signal_handlers.py:12` filters pending VCFs on it, so an AnnotSV completion matches no
`UploadedVCFPipelineMaxVariant` row and no-ops. No change needed there.

### 3.7 Supplementary runs own no `VariantAnnotation` rows

The `VariantAnnotation.annotation_run` FK keeps pointing at the VEP run that created the row. An AnnotSV
run only ever UPDATEs. Three useful consequences:

- `AnnotationRun.delete_related_objects()` on an AnnotSV run deletes nothing, so `reset_for_retry`,
  `_reset_run_for_redispatch` and deleting the run are all safe and cannot take VEP's data with them.
- Re-running AnnotSV over a range is idempotent.
- Provenance is the range lock: *"has AnnotSV covered variant X?"* is *"is there a FINISHED `ANNOTSV`
  `AnnotationRun` on the range lock containing X, at which pipeline version?"* No per-variant table.

  This holds because AnnotSV *sweeps* - one run takes every SV in its range. It is not a property of
  supplementary pipelines in general: a pipeline run on a chosen handful of variants at a time covers a
  range lock only partially, so it needs per-variant records instead. See #1755.

This does mean `import_annotsv_tsv` must stop scoping by `annotation_run` — it currently does
`VariantAnnotation.objects.filter(annotation_run=annotation_run)`
(`annotation/vcf_files/bulk_annotsv_tsv_inserter.py:160`), which will match zero rows once the AnnotSV run
is its own run. Scope by annotation version and the lock's variant range instead:

```python
    range_lock = annotation_run.annotation_range_lock
    qs = VariantAnnotation.objects.filter(version=annotation_run.variant_annotation_version,
                                          variant_id__gte=range_lock.min_variant_id,
                                          variant_id__lte=range_lock.max_variant_id)
```

The `version` filter is what keeps the update off other annotation versions' rows — it is load-bearing,
not cosmetic, and it is what the `annotation_run` filter was providing implicitly.

Deliberately *not* `get_queryset_for_annotation_version`: that adds the SQL transformer which rewrites
joins to the version's partition table, and this query reads `VariantAnnotation` directly rather than
joining through `Variant`. Inheritance means a query on the parent reaches the partitions anyway (which is
what the original `annotation_run` filter relied on), so the explicit `version` filter is both correct and
the smaller change.

### 3.8 Tool version pinning

`VariantAnnotationVersion.annotsv_code` / `annotsv_bundle` are deleted — pinning AnnotSV on the VEP version
is the thing that makes a bundle bump cost a full re-annotation. Replaced by one small model:

```python
# annotation/models/models.py

class AnnotationPipelineVersion(TimeStampedModel):
    """ Version of a non-VEP pipeline's tool + data bundle, recorded on each run that used it. Separate
        from VariantAnnotationVersion so a tool can roll without forcing a VEP re-annotation - which is
        the point of #720. VEP pipelines don't use this: their version IS the VariantAnnotationVersion. """
    pipeline_type = models.CharField(max_length=1, choices=VariantAnnotationPipelineType.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    code_version = models.TextField()                       # eg AnnotSV "3.5.8"
    data_version = models.TextField(null=True, blank=True)  # eg AnnotSV bundle version

    class Meta:
        unique_together = ("pipeline_type", "genome_build", "code_version", "data_version")

    def __str__(self):
        data = f"/{self.data_version}" if self.data_version else ""
        return f"{self.get_pipeline_type_display()} {self.code_version}{data} ({self.genome_build})"
```

```python
# AnnotationRun
    # What ran, for non-VEP pipelines. Null for VEP pipelines (VariantAnnotationVersion covers those).
    pipeline_version = models.ForeignKey(AnnotationPipelineVersion, null=True, blank=True,
                                         on_delete=PROTECT)
```

The runner resolves the current version at `annotate()` time (`AnnotSV -version` plus
`settings.ANNOTATION_ANNOTSV_BUNDLE_VERSION`, `get_or_create`) and records it. This **replaces**
`annotsv_check_command_line_version_match` rather than reimplementing it: recording what actually ran is
more useful than refusing to run, because it makes a mid-backfill upgrade visible and re-runnable instead
of turning every remaining run into an error.

Rolling a version is then: upgrade the binary/bundle, then reset the stale runs.

```python
# annotation/management/commands/annotation_pipeline_rerun.py
#   --pipeline-type A --genome-build GRCh38 [--stale] [--limit N]
#
# --stale selects FINISHED runs whose pipeline_version is not the current one, and calls
# reset_for_retry() on each (safe - a supplementary run owns no annotation rows, §3.7), then kicks the
# dispatcher. Also the tool for "AnnotSV was off, now it's on, re-run the ones that already finished
# empty". Replaces manage.py annotsv_run, which is deleted.
```

Surface `pipeline_version` on the annotation run detail page and the runs grid, so "which AnnotSV produced
this" is answerable without the shell.

### 3.9 Abortability and the AnnotSV timeout

From the #720 comment on the #1658 work: abortability belongs in the runner contract, not per-tool.
`run_annotsv` currently calls `subprocess.run` directly (`annotation/annotsv_annotation.py:54`), so a
reclaimed worker's AnnotSV keeps running to `ANNOTATION_ANNOTSV_TIMEOUT_SECONDS` while the new attempt
works the same range.

Route AnnotSV through `execute_cmd(process_callback=lease_heartbeat.set_process)` like VEP. The blocker
noted in that comment is that `execute_cmd` has a bare `communicate()` with no timeout. Add one:

```python
# library/utils/os_utils.py

def execute_cmd(cmd: list, **kwargs) -> CmdOutput:
    process_callback = kwargs.pop("process_callback", None)
    # timeout (#720): kill + drain rather than leaving an orphan holding the pipes, matching what
    # subprocess.run does. Callers get the partial output and a non-zero code, so the existing
    # `return_code != 0` handling covers a timeout with no extra branch.
    timeout = kwargs.pop("timeout", None)
    ...
    try:
        std_out, std_err = pipes.communicate(timeout=timeout)
        return_code = pipes.returncode
    except subprocess.TimeoutExpired:
        pipes.kill()
        std_out, std_err = pipes.communicate()
        return_code = pipes.returncode or -1
        logging.warning("Command timed out after %ss: %s", timeout, cmd)
```

`run_annotsv` then loses its own `TimeoutExpired` handling and takes `process_callback` / `timeout`
through. Every current `execute_cmd` caller is unaffected (`timeout=None` is `communicate()`'s default).

### 3.10 AnnotSV runner

```python
# annotation/pipelines/annotsv.py

class AnnotSVRunner(AnnotationPipelineRunner):
    pipeline_type = VariantAnnotationPipelineType.ANNOTSV
    selects_unannotated = False  # updates the rows the SV VEP run wrote

    def get_variants_qs(self, annotation_run):
        # Deliberately no ANNOTATION_VEP_SV_MAX_SIZE filter. That cap exists because VEP fills the logs
        # with 'too long to annotate'; AnnotSV handles large SVs, and ranking them is what it is for.
        # Feeding AnnotSV the VEP dump is why it has never seen them.
        return super().get_variants_qs(annotation_run)

    def annotate(self, annotation_run, lease_heartbeat=None):
        annotation_run.pipeline_version = get_or_create_annotsv_version(annotation_run.genome_build)
        dump_count = dump_variants(annotation_run, runner=self, task_token=annotation_run.task_id)
        if not dump_count:
            annotation_run.annotated_count = 0
            annotation_run.annotation_end = timezone.now()
            annotation_run.save()
            return

        annotsv_dir = get_annotsv_dir(annotation_run)
        cmd = get_annotsv_command(annotation_run.vcf_dump_filename, annotsv_dir,
                                  annotation_run.genome_build, annotation_run.annotation_consortium)
        annotation_run.annotation_start = timezone.now()
        annotation_run.pipeline_command = " ".join(cmd)
        annotation_run.save()

        os.makedirs(annotsv_dir, exist_ok=True)
        process_callback = lease_heartbeat.set_process if lease_heartbeat else None
        try:
            return_code, std_out, std_err = execute_cmd(
                cmd, process_callback=process_callback,
                timeout=settings.ANNOTATION_ANNOTSV_TIMEOUT_SECONDS)
        finally:
            if lease_heartbeat:
                lease_heartbeat.set_process(None)

        max_output = 1_000_000
        annotation_run.pipeline_stdout = std_out[:max_output] if std_out else std_out
        annotation_run.pipeline_stderr = std_err[:max_output] if std_err else std_err

        tsv = get_annotsv_tsv_filename(annotation_run.vcf_dump_filename, annotsv_dir)
        if return_code != 0 or not os.path.exists(tsv):
            # Its own run now, so failing it is isolated - the SV VEP run's data is already committed.
            # As VEP: persist the output only while we still own the run (#1658).
            if not annotation_run.save_if_owner(annotation_run.task_id):
                logging.warning("AnnotationRun %s was reclaimed - discarding AnnotSV output from task %s",
                                annotation_run.pk, annotation_run.task_id)
            missing = "" if os.path.exists(tsv) else f" (expected TSV not found: {tsv})"
            raise RuntimeError(f"AnnotSV returned {return_code}{missing}")

        annotation_run.vcf_annotated_filename = tsv
        annotation_run.annotation_end = timezone.now()
        annotation_run.save()

    def import_results(self, annotation_run):
        annotation_run.upload_start = timezone.now()
        annotation_run.upload_attempts += 1
        annotation_run.save()
        annotation_run.annotated_count = import_annotsv_tsv(annotation_run)
        annotation_run.upload_end = timezone.now()
        annotation_run.save()

    def get_output_paths(self, annotation_run, dump_filename) -> list[str]:
        return [dump_filename, get_annotsv_tsv_filename(dump_filename, get_annotsv_dir(annotation_run))]
```

The best-effort swallow goes away. An AnnotSV failure now fails an AnnotSV run — visible on the runs page,
retryable from the UI, subject to the attempt cap — while the SV VEP run stays FINISHED with its data
imported. That is what "isolate failures" means here, and it is strictly better than a `annotsv_error`
text field nothing looks at.

`get_annotsv_dir` keeps being keyed on the run, not the task, for the reason recorded on the issue: the TSV
inside is named from the per-task dump stem, so concurrent attempts write side by side, while a per-task
directory would be nameable only by the attempt that created it.

### 3.11 VEP runner

`VEPRunner` is `dump_and_annotate_variants` minus the AnnotSV block, parameterised by pipeline type:

```python
# annotation/pipelines/vep.py

class VEPRunner(AnnotationPipelineRunner):
    def __init__(self, pipeline_type):
        self.pipeline_type = pipeline_type

    def get_variants_qs(self, annotation_run):
        qs = super().get_variants_qs(annotation_run)
        if self.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
            if settings.ANNOTATION_VEP_SV_MAX_SIZE:
                # VEP skips variants above a certain size and fills the logs with 'too long to annotate'
                q_not_too_long = Q(svlen__isnull=True) | Q(abs_svlen__lte=settings.ANNOTATION_VEP_SV_MAX_SIZE)
                qs = qs.annotate(abs_svlen=Abs("svlen")).filter(q_not_too_long)
        return qs

    def annotate(self, annotation_run, lease_heartbeat=None):
        # existing dump_and_annotate_variants body, minus lines 534-559 (the AnnotSV stage)
        ...

    def import_results(self, annotation_run):
        import_vcf_annotations(annotation_run)   # minus its AnnotSV block, lines 79-88
```

The SV-conservation sidecar stage (`annotate_variants.py:562-580`) stays inside `VEPRunner.annotate` for
now — it writes into `_max` columns that `BulkVEPVCFAnnotationInserter` merges during the VEP import, so it
is genuinely part of that import, not a separate pipeline. It is a good candidate for its own pipeline type
later, on exactly this machinery.

### 3.12 Loose ends in existing code

- `AnnotationRun.PIPELINE_TYPE_DESC` (`models.py:1330`) maps only STANDARD and STRUCTURAL_VARIANT — add
  `"gene_level"` and `"annotsv"`, so dump stems stop falling back to the bare enum character.
- `annotation/views.py:334` builds a VEP command for **every** pipeline type to display on the annotation
  versions page — already wrong for `GENE_LEVEL`. Iterate `vep_pipeline_types()`.
- `_reset_run_for_redispatch` (`annotation_scheduler_task.py:507-580`) derives VEP-shaped paths from the
  dump stem. Route through `get_pipeline(run.pipeline_type).runner`: `get_output_paths` for the discard,
  and a runner predicate for "did the tool finish?" (VEP: the `--skipped_variants_file` exists per #1710;
  AnnotSV: the TSV exists). Drop the `annotsv_imported = False` line.
- `annotation/signals/annotation_run_cleanup.py` `get_all_run_output_paths` — replace its hardcoded VEP +
  AnnotSV path set with `runner.get_output_paths(...)` plus the persisted `vcf_annotated_filename`, and
  keep the `remove_annotsv_dir` behaviour keyed on the AnnotSV pipeline type.
- `manage.py annotsv_run` is deleted (§3.8) — the normal run-retry path replaces it.

---

## 4. Migrations

Greenfield: no deployment carries AnnotSV data, so no data migration.

```python
operations = [
    migrations.CreateModel(name="AnnotationPipelineVersion", fields=[...]),
    migrations.AddField("annotationrun", "pipeline_version",
                        models.ForeignKey(null=True, blank=True, on_delete=PROTECT,
                                          to="annotation.annotationpipelineversion")),
    migrations.AlterField("annotationrun", "pipeline_type", ...),          # new choice
    migrations.AlterField("uploadedvcfpipelinemaxvariant", "pipeline_type", ...),  # new choice
    migrations.RemoveField("annotationrun", "annotsv_tsv_filename"),
    migrations.RemoveField("annotationrun", "annotsv_error"),
    migrations.RemoveField("annotationrun", "annotsv_imported"),
    migrations.RemoveField("variantannotationversion", "annotsv_code"),
    migrations.RemoveField("variantannotationversion", "annotsv_bundle"),
]
```

`unique_together ("annotation_range_lock", "pipeline_type")` is unchanged and already correct for the new
type. No column on `AnnotationRun` changes meaning, so existing rows need nothing.

## 5. Settings

`ANNOTATION_ANNOTSV_ENABLED` keeps its name and gains a slightly wider meaning — it now gates whether the
scheduler creates `ANNOTSV` runs at all:

```python
# variantgrid/settings/components/annotation_settings.py

# Schedules an ANNOTSV AnnotationRun per range lock, dependent on that lock's STRUCTURAL_VARIANT run.
# Turning this on for an existing deployment backfills: the scheduler creates ANNOTSV runs for every
# existing range lock, and the count lane finishes the (majority) SV-free ones without dispatching.
ANNOTATION_ANNOTSV_ENABLED = False
```

`ANNOTATION_ANNOTSV_CODE_VERSION` is not added — the runner reads the binary's own `-version`.
`ANNOTATION_ANNOTSV_BUNDLE_VERSION` stays and becomes `AnnotationPipelineVersion.data_version`.
`ANNOTATION_VEP_BUFFER_SIZE` stays keyed by pipeline type; `ANNOTSV` simply has no entry.

## 6. Tests

Worth keeping (each covers logic we wrote, with a branch that is easy to get wrong later):

- `test_annotsv_run_waits_for_vep_run` — an `ANNOTSV` run on a lock whose `STRUCTURAL_VARIANT` run is not
  FINISHED is absent from `_dispatchable_runs_qs`, and appears once it is.
- `test_sv_free_range_finishes_annotsv_run_without_dispatching` — the cheap-backfill property of §3.5:
  an SV-free range lock finishes its AnnotSV run on the count lane, with no dependency wait.
- `test_annotsv_run_does_not_block_vcf_import` — a VCF with SVs is fully annotated once its
  `STRUCTURAL_VARIANT` run finishes, with the `ANNOTSV` run still pending (#1656's guarantee, extended).
- `test_enabling_annotsv_creates_runs_for_existing_locks` — the backfill claim in §2, driven through
  `_handle_variant_annotation_version`.
- `test_annotsv_failure_isolated` — an AnnotSV run erroring leaves the `STRUCTURAL_VARIANT` run FINISHED
  with its rows intact.
- `test_annotsv_does_not_inherit_veps_sv_size_cap` — the `abs_svlen` annotation VEP adds solely to apply
  `ANNOTATION_VEP_SV_MAX_SIZE` is absent from the AnnotSV queryset.
- `test_annotsv_runs_not_created_while_disabled` — the other half of the enable/backfill switch.
- `test_execute_cmd_timeout` — kills, drains, and returns non-zero with partial output.
- `test_annotsv_run_owns_no_annotation_rows` — §3.7's safety condition: the AnnotSV run owns no
  `VariantAnnotation` rows, so resetting or deleting it cannot take the VEP run's data with it.
- `test_reuses_row_for_same_code_and_bundle` / `test_upgraded_binary_gets_its_own_version` — §3.8's
  record-don't-enforce behaviour.

Extend the existing `upload/tests/test_pipeline_max_variant.py` rather than duplicating it — it already
covers `pipeline_type_variant_q` and the per-type watermark, and only needs the blocking-types distinction
added.

## 7. Implementation order

Each step is independently reviewable and leaves the tree green. **All six landed** — see §9 for the
places the implementation departed from this plan.

1. **`execute_cmd` timeout** (`library/utils/os_utils.py` + test). Self-contained, no annotation changes.
2. **Extract `annotation/annotation_run_files.py`** — pure move of the path/dump helpers out of
   `annotation/tasks/annotate_variants.py`, fixing the `annotation_run_cleanup` → tasks import. No
   behaviour change.
3. **Introduce `annotation/pipelines/`** with `VEPRunner` and `GeneLevelRunner` wrapping today's code, and
   route both celery tasks through the registry. AnnotSV still inline in `VEPRunner`. Pure refactor —
   existing tests are the check.
4. **Add the `ANNOTSV` pipeline type**: enum, `pipeline_type_variant_q`, `AnnotSVRunner`, registry entry,
   dependency gating in the dispatcher and count lane, `blocking_pipeline_types()` in the upload app,
   `AnnotationPipelineVersion`, migration; strip the AnnotSV stage from `VEPRunner` /
   `import_vcf_annotations`; rescope `import_annotsv_tsv`; delete `manage.py annotsv_run`.
5. **`annotation_pipeline_rerun` command** + `pipeline_version` on the runs grid and detail page.
6. **Tests** per §6.

Step 3 is the one worth pausing on for review — it touches the hot path of every annotation run without
changing any behaviour, so it should go in on its own and be seen to be a no-op.

---

## 8. The other issues on #720

**In scope, folded in:**

- The #1658 follow-up note on the issue — abortability in the runner contract, and the `execute_cmd`
  timeout that AnnotSV needed before it could use `process_callback`. §3.9.
- The #1660 follow-up note — `get_annotsv_dir` stays keyed on the run, and the reasoning is preserved in
  the code that moves. §3.10.

**Out of scope, deliberately:**

- **SpliceAI (#720's original ask) — its own issue, #1755, and a different shape from this.** SpliceAI
  runs are sparse and manual: a handful of specially chosen variants, with the *same* range lock possibly
  touched again months later. That breaks the range-lock-as-provenance contract in §3.7 outright — a
  FINISHED run over a lock would not mean SpliceAI had covered the lock. So SpliceAI needs per-variant
  coverage records, and a manual run is a poor fit for `AnnotationRun` at all (which is coupled to a range
  lock through `variant_annotation_version`, `genome_build`, `annotation_consortium` and
  `_trigger_dispatch`, and which #1654 deliberately keeps from ever being committed rangeless).

  The sharpest consequence: `spliceai_max_ds IS NULL` cannot distinguish "never ran" from "ran, no splice
  effect", and SpliceAI legitimately returns nothing for variants outside genes or over `2 * -D` long.
  Without a per-variant record the expensive tool is re-run forever on variants already known to score
  nothing.

  A batch mode sweeping range locks *would* fit this framework, but only if the gap set is small — the
  benchmarks on #720 are ~9.5s CPU per variant, and a Quadro P600 was no faster than 8 CPU threads. That
  needs counting before it is worth designing. Either way none of it should hold up the AnnotSV split.

- **#1675 (backfill a single `VariantAnnotation` column).** Different problem — re-running VEP over
  *already-annotated* variants to fix one column, e.g. the COSMIC field renames in #1673. It shares one
  seam with this work: `selects_unannotated = False`. But it also needs a cut-down VEP command line and an
  UPDATE-a-column-subset import mode, neither of which AnnotSV requires (AnnotSV's import is already a
  targeted `bulk_update`). Doing it as a pipeline type on this registry afterwards is a good fit; doing it
  now would double the size of the change for no gain to either issue.

- **#1653 (short/long annotation queues).** Enabled by this, not part of it: once a pipeline is a registry
  entry, a `queue` field on `PipelineDef` routes it, and `_lease_and_launch_run` picks the task by lane
  already. Worth revisiting after step 3 lands.

- **#1679 (latch VEP version mismatch on the VAV).** Adjacent — this plan replaces the AnnotSV version
  *check* with version *recording* (§3.8), which is the opposite move to #1679's. They do not conflict:
  the VEP mismatch is a deployment misconfiguration that should stop scheduling, whereas an AnnotSV
  upgrade mid-backfill is a normal event that should be recorded and re-run. Keep them separate.

- **#1699 (mappability).** Not a supplementary run at all — it is a VEP `--custom` bigWig track, so it
  belongs in `vep_columns.py` and the VEP command, exactly as the issue describes. Listing it here only to
  say it needs nothing from this work.

---

## 9. Where the implementation departed from this plan

Three things only became clear while building it:

1. **The count-lane "trap" was not real** (§3.5). `annotated=True` drops the unannotated filter rather
   than inverting it, so a supplementary run's count never depended on the VEP run having imported. The
   dependency gate was written, then removed: it delayed the SV-free-range fast path for no correctness
   benefit. Dependency gating lives only in `_dispatchable_runs_qs`, where it is genuinely required.

2. **Scoping the AnnotSV import by `version`, not the partition transformer** (§3.7).
   `get_queryset_for_annotation_version` rewrites joins to the version's partition table, which is right
   for querysets that join through `Variant` but returns nothing for a direct `VariantAnnotation` read in
   the test environment (rows written by the ORM land in the parent table). Filtering `version=` directly
   is both correct and closer to what the code it replaces was doing.

3. **The runner contract needed `tool_finished` + `record_resume_state`.** `_reset_run_for_redispatch`
   made two VEP-shaped decisions the plan had not accounted for: "did the tool finish?" (VEP: the
   `--skipped_variants_file` exists, per #1710; AnnotSV: the TSV exists) and "what does the run need to
   resume upload-only?" (VEP: annotated VCF + skipped-variants file; AnnotSV: the TSV). Both are per-tool,
   so both became runner methods rather than branches in the scheduler.

Also delivered but not in the plan: `pipeline_version` on the annotation-runs grid and the run detail
page, and the `views.py` VEP-command loop restricted to `vep_pipeline_types()` (it was already wrong for
`GENE_LEVEL`).
