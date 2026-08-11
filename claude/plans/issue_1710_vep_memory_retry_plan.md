# Issue #1710 — VEP memory failures: contain the blast radius, size the buffer, retry once

AnnotationRun 34818 (GRCh38, standard, 6,250 variants) was OOM-killed three times running VEP at the
default buffer size of 500, peaking at 5.3 GB and 7.3 GB RSS. Because `celeryd_annotation_workers.service`
has no memory cap, both kills were *global* OOMs — the kernel weighed every process on the box — and the
2026-08-11 01:18 kill stalled four other concurrent runs badly enough that `reclaim_stalled_annotation_runs`
swept them up mid-VEP, marked them ANNOTATION_COMPLETED off a half-written `.vcf.gz`, and sent them to the
import lane where they failed with `TruncatedVEPOutputError`. One bad run cost five.

This plan has three independent halves, in descending order of value:

1. **Contain** — a `MemoryMax=` on the annotation worker unit, so a memory-hungry VEP kills only itself.
2. **Size** — the buffer size becomes a per-range value rather than a per-pipeline-type constant, chosen
   from the load actually measured at dump time.
3. **Retry** — a run whose VEP output came back short is redispatched once with the buffer at the floor.

**Halves 2 and 3 may be unnecessary. Read the next section first.**

Related: #1708 (dispatching long delins to the SV pipeline by length) reduces load in the megabase case but
is not the fix here; #1701 supplies the truncation detection this plan retries off.

---

## First, the one-line version — done

The variants in run 34818 were ~1 kb, which (see "What actually costs the memory") rules out the
per-variant explanations and leaves bins-per-buffer — driven by how *sparse* the range is — as the driver.
That is precisely what `--buffer_size` bounds, so the standard pipeline's default is now 100:

```python
_VARIANT_ANNOTATION_PIPELINE_STANDARD: 100,
```

Applied ahead of the rest of this plan because it is one constant, does not change VEP's output, and does
not mint a new `VariantAnnotationVersion` — buffer size is not among the kwargs
`vep_dict_to_variant_annotation_version_kwargs` builds, so there is no re-annotation and reverting is a
one-line edit.

The SV pipeline stays at 250 **for now** — but it is not obviously right either, and probably for stronger
reasons. See "The SV lane has the same question open" below.

### Why a smaller buffer looks close to free

`get_all_features_by_InputBuffer` (`AnnotationSource.pm:109-145`) ends each buffer with
`$self->clean_cache($regions)`, which keeps **only the current buffer's regions** and evicts the rest. Since
`write_qs_to_vcf` dumps position-sorted within each contig, the regions two consecutive buffers share stay
cached across the boundary, and a region is generally entered once per run whatever the buffer size. Total
disk region loads are therefore roughly the number of distinct regions in the dump — **independent of
buffer size**. That is the mechanism behind the "runtime was flat across every size measured; output
byte-identical" note already in `ANNOTATION_VEP_BUFFER_SIZE`.

### So why is the default 500?

Not because lower was expensive — because lower was *pointless* for the data it was measured on. The
recorded benchmark is `500 -> 871MB, 250 -> 858MB` on 12.5k small variants, hence "flattens at 500 - below
that fixed startup dominates". Those variants were dense, so 500 consecutive ones sat in a handful of bins
and bins-per-buffer was never the dominant term. A sparse range inverts that: 500 consecutive variants
scattered across a build can touch several hundred distinct bins, and 500 -> 100 cuts that close to
five-fold.

### The one measurement still worth taking

Re-dump 34818's range and run VEP over it at buffer 100 under `/usr/bin/time -v`. One number: peak RSS.

Not a matrix — the existing benchmark already covers runtime ("flat across every size measured; output
byte-identical") and the `clean_cache` mechanism explains why. What this run answers is whether the buffer
was the driver *at all*. If peak RSS is still multi-GB, the transcript-count × allele-length term is the
real cause (see "One candidate not yet ruled out"), `--transcript_filter` is the lever, and lowering the
default fixed nothing — which is worth knowing before the ticket is closed.

Two things to note while doing it:

- **Mind `--fork`.** `Runner.pm:465` computes `maxForkSize = int($buffer_size / (2 * $fork_number)) || 1`,
  so the buffer is divided among forks and needs to stay comfortably above twice the fork count. VG runs
  `ANNOTATION_VEP_FORK = 1`, so this is a settings comment rather than a blocker — but a deployment that
  raises fork must raise the buffer with it.
- **Per-buffer overhead × 5, on a realistically-sized run.** This is the one real cost of the change, and
  34818 understates it: at 6,250 variants it is 13 buffers vs 63, but standard runs measured on the dev box
  average ~40,000 variants (max 50,001), which is ~80 buffers vs ~400. Benchmark a full-sized dump as well
  as 34818's, or the multiplier that matters goes unmeasured.

### The SV lane: closed, the setting is inert

SV runs are sparser than standard ones, not denser — an `AnnotationRangeLock` is sized for 5,000–25,000
*total* variants, of which only a handful are SVs — so the concern was that bins-per-buffer there is driven
by the distance between variants **and** each variant's own span. And the 250 was tuned on data with the
same flaw as the 500: "1000 DELs (median 49kb) -> worst-case 247 regions in one buffer" at buffer 1000 is
roughly 4 SVs per bin, i.e. a *clustered* set.

Measured instead of argued (dev box, 2026-08-12), and the answer is that none of it matters:

| SV runs (`pipeline_type='C'`) with a dump_count | 173 |
|---|---|
| max dump_count | **24** |
| average | 1.66 |
| dumped zero variants | 147 |
| runs at or above the 250 buffer | **0** |

No SV run is within an order of magnitude of the buffer size, so `--buffer_size` never binds on that lane —
every SV run is already a single buffer, and 250 vs 100 is a distinction without a difference. **Leave the
SV setting at 250.**

The structural reason generalises past this box: exceeding 250 SVs in one run needs SVs to be >1–5% of all
variants in a range lock. Worth confirming on the deployment that hit #1710, but the same query answers it
in seconds.

Two points kept because they still bound any future change there:

- **The re-load cost is real for SVs.** The standard-pipeline argument (position-sorted variants never
  revisit an evicted bin) does not carry over: a long SV loads bins 10-60, a short one in the next buffer
  keeps only bin 12 and `clean_cache` evicts the rest, and a later overlapping SV re-fetches them. This is
  what the existing "250 keeps enough per buffer to avoid re-loading shared regions" note protects.
- **There is a floor buffer size cannot cross.** A 5 Mb deletion is six bins on its own at any buffer size.
  Below that floor the only lever is `--max_sv_size`, which is version-affecting (decision 7).

### What this leaves

Phase 1 still lands — containment is independent of what caused the pressure. Phase 2's per-range column is
still worth having as a **manual operator lever** (a range known to be bad, set by hand). Phase 3's
automatic retry has lost its purpose: with the default already at the floor there is no rung to step down
to, and a run that still dies should fail visibly. Phase 4 goes entirely — predicting a buffer size only
pays when the default is too high to be safe everywhere, which is what lowering it fixed.

Phases 3 and 4 are kept below as a record of the design and of why they are not being built.

---

## Background: two shapes of the same failure

`-9` is not a sufficient signal, and neither is the exit code in general.

| what dies | VEP exit code | where it surfaces |
|---|---|---|
| `perl` is the OOM victim | `-9` | `RuntimeError: VEP returned -9` at the end of the VEP stage |
| the `gzip` child dies | **0** | nothing at VEP time — only at import, via the #1701 checks |

The second shape is Ensembl/ensembl-vep#2033 (open, no fix): with `--compress_output` VEP pipes to an
external compressor and discards `close()`'s result, so it exits 0 having written a truncated file.
`check_all_dumped_records_accounted_for` (`annotation/vcf_files/import_vcf_annotations.py:142`) documents
exactly this.

`-9` is also ambiguous in the other direction: `AnnotationRunLeaseHeartbeat._abort_process`
(`annotation/tasks/annotate_variants.py:138`) calls `process.kill()` deliberately, so a lease abort and an
OOM kill are indistinguishable by return code alone. The heartbeat knows which it did, so the fix is to
ask it.

So detection keys on **the output being short**, not on why:

- `TruncatedVEPOutputError` from the import lane — covers both shapes, and is the gate.
- A negative return code that our own heartbeat did not cause — a fast path for the first shape only,
  saving a wasted dispatch to the import lane.

Per the #1710 comment, the flag and the retry cap are therefore named for **truncation**, not for OOM:
`TruncatedVEPOutputError` says the output was short, not that memory caused it. Lowering the buffer is
cheap and harmless whatever the cause, so acting on the broader signal is the honest trade.

---

## What actually costs the memory

Read out of the VEP 115 source installed at `/data/annotation/VEP/vep_code/115`. **The variants in run 34818
were around 1 kb**, which rules out most of the per-variant story and points back at the buffer.

### What 1 kb rules out

At 1 kb, none of VEP's span-proportional work distinguishes these variants from a SNV:

- **Transcript cache bins.** `get_regions_from_coords` (`AnnotationSource.pm:190`) loads every 1 Mb bin from
  `start - up_down_size` to `end + up_down_size`. At 1 kb that is one bin, occasionally two — the same as a
  SNV. It only becomes a per-variant problem in the megabase range.
- **Plugin data fetches.** `dbNSFP.pm:337` asks for `$vf->{start} - 1 .. $vf->{end}`, but
  `BaseVepTabixPlugin` has `$DEFAULT_EXPAND_RIGHT = 1e6`, so the region actually fetched and cached is
  **1 Mb wide whatever the variant's size**. A 1 kb indel and a SNV pull the same block. (That expansion is
  deliberate — it is the region cache that stops sequential nearby variants re-querying — but it means
  plugin memory is a function of position density, not variant length.)
- **`--check_existing`.** `VariationTabix.pm:229` queries `chr:start-1-end+1` per variant: a 1 kb window.

So there is no single pathological variant here, and `--max_sv_size` would not have caught these at any
plausible setting.

### What is left: bins per buffer, driven by sparsity

The mechanism that survives is the one the ticket started with. `get_all_regions_by_InputBuffer`
(`AnnotationSource.pm:219`) computes bins over the **whole current buffer**, and holds the transcript cache
for all of them at once (dropping the rest in `AnnotationSource::clean_cache`). The count is therefore
"how many distinct 1 Mb bins do 500 consecutive variants touch" — which is about **how spread out the
variants are**, not how long they are.

That is the difference between run 34818 and the 12.5k-small-variant benchmark recorded in
`ANNOTATION_VEP_BUFFER_SIZE`. A dense VCF puts 500 consecutive variants inside a handful of bins. A sparse
range — 6,250 variants scattered across a build — can put 500 consecutive variants in several hundred
distinct bins, each carrying its own transcript cache, and each also being a 1 Mb plugin block under the
expansion above.

This is good news for the plan: **it is exactly what buffer size controls**, and exactly what Phase 4
measures. Cutting 500 → 100 cuts bins-per-buffer close to five-fold, and the existing benchmark says that
costs nothing in runtime.

### One candidate not yet ruled out

A 1 kb allele is copied per overlapping transcript, and for a frameshift VEP computes downstream peptide
sequence per transcript — so cost scales with **transcript count × allele length**. That is the same shape
as the #1228 BRCA1 blowout (368 RefSeq transcripts for one gene), which `--transcript_filter` already exists
to mitigate. It is a plausible contributor and it is *not* addressed by buffer size. Untested — see the
bisection below.

### Settle it by measurement, not by reading

Cheap and decisive, because the input is known: re-dump 34818's range (or reuse its dump if the per-task
file survives) and run VEP over the same 6,250 variants under `/usr/bin/time -v`, changing one thing at a
time — baseline, `--buffer_size 100`, plugins removed, `--no_check_existing`, a tighter
`--transcript_filter`. Peak RSS deltas name the driver in an afternoon of compute. Record
`dump_regions_spanned` (Phase 4) for the range at the same time: if it is in the hundreds, sparsity is
confirmed and Phase 4's model is calibrated on the case that actually failed rather than on the SV
benchmark.

**Do this before implementing Phase 4** — it decides whether the regions model is the right one, and what
`ANNOTATION_VEP_MAX_REGIONS_PER_BUFFER` should be.

### Why splitting the range still does not help

`get_all_regions_by_InputBuffer` computes bins over the current buffer — a sliding window of N consecutive
variants. Halving an `AnnotationRangeLock` halves the run's total span but leaves every window inside it
identical, so bins-per-buffer is unchanged. This is the mechanism behind the ticket's assertion that
`subdivide_annotation_range_lock` is not a memory fix, and it holds for splitting by genomic region too:
what matters is which variants land in a window together, not how wide a net the run casts. Shrinking the
window is the lever; shrinking the range is not.

### The other flags, and what they cost

- **`--cache_region_size` is not tunable.** `CacheDir.pm:149` prefers the cache's own `info.txt` value over
  the parameter, and the downloaded cache is built in 1 Mb bins.
- **`--max_sv_size` is length-based, not SV-based** — `Parser.pm:490` computes `$len = $vf->{end} -
  $vf->{start}` for every variant and sets `vep_skip`, which `AnnotationSource.pm:238` honours before any
  bin is loaded. Worth knowing, but irrelevant at 1 kb, and version-affecting: `sv_max_size` is one of the
  kwargs `get_or_create_variant_annotation_version_from_current_vep` keys its `get_or_create` on
  (`annotation/vep_config.py:83`), so changing it re-annotates everything.
- **`--check_existing`** could be dropped for a pipeline, at the cost of COSMIC IDs. Small at 1 kb.
- **Plugins have no per-variant size guard**, and their fetches are 1 Mb regardless, so the only lever on
  them is which pipeline runs them (`get_vep_command` adds the plugin block for STANDARD only).
- **`--transcript_filter`** is already used for the BRCA1 blocklist (#1228) and is the lever for the
  transcript-count term, if the bisection implicates it.

---

## Phase 1 — Contain: cap the worker unit

`config/systemd/celeryd_annotation_workers.service` has no memory accounting at all. Adding one turns a
global OOM into a cgroup OOM: VEP dies, postgres and gunicorn are never candidates, and the kill becomes
attributable rather than inferred.

```ini
[Service]
...
# Fraction of installed RAM, resolved by systemd - these units are copied to /lib/systemd/system
# verbatim (config/systemd/README.txt), so an absolute figure would ship the same number to a 16GB
# box and a 256GB one. Leaves room for postgres, gunicorn and page cache: only this unit is capped,
# so a cgroup OOM fires only when the annotation pool alone crosses the line.
MemoryMax=50%
# Kill promptly rather than thrash. The collateral damage in #1710 was four runs *stalled* past their
# 900s lease, not killed - swap thrash is what produces that, and reclaim then swept them up mid-VEP.
MemorySwapMax=0
# A cgroup OOM kill must take the offending VEP only. systemd's default OOMPolicy is `stop`, which -
# with Restart=always above - would stop and restart the whole worker pool, killing every concurrent
# run: exactly the collateral damage this cap exists to prevent.
OOMPolicy=continue
```

`MemoryAccounting=` is implied by setting a limit (and defaults on in modern systemd), so it isn't listed.

Three things to get right:

- **`OOMPolicy=continue` is not optional.** The unit is `Restart=always`, so leaving the default would
  convert "one VEP killed" into "the pool bounced", which is worse than today.
- **A percentage, not a byte count.** systemd takes `MemoryMax=` as a percentage of installed physical
  memory, so the cap scales with the box without any RAM detection of our own - no `/proc/meminfo` read,
  no new dependency, and nothing to template at install time. It is also the same shape postgres is
  already sized in (`shared_buffers` as a fraction of RAM), so the two scale together instead of drifting
  apart per deployment. A box that wants a different figure overrides it in a drop-in
  (`/etc/systemd/system/celeryd_annotation_workers.service.d/override.conf`) rather than editing the
  repo's unit.
- **The cap is per-pool, and the per-VEP share tracks CPUs, not RAM.** The unit runs `celery multi`
  (`Type=forking`), all workers in the one cgroup, so `MemoryMax=` bounds the sum of concurrent VEPs plus
  the workers themselves. `config/celery/celeryd_annotation_workers.env` sets no `--concurrency`, so
  celery's prefork default is the CPU count: the budget one VEP actually gets is `MemoryMax / concurrent
  VEPs`, whose divisor scales with cores while the numerator scales with RAM. On a high-core, modest-RAM
  box that share is thinner than the percentage suggests - which is an argument for setting
  `--concurrency` explicitly there rather than for raising the cap.

Once the cap is in place, `memory.events` becomes usable for attribution:

```python
# annotation/vep_annotation.py

CGROUP_MEMORY_EVENTS = "/sys/fs/cgroup/system.slice/celeryd_annotation_workers.service/memory.events"


def get_cgroup_oom_kill_count() -> Optional[int]:
    """ The annotation worker cgroup's cumulative oom_kill counter, or None where it isn't readable
        (no MemoryMax set, cgroup v1, not running under the unit - eg a dev box or a unit test).

        Only a global-OOM box leaves this at zero while VEP is being killed, which is the state Phase 1
        exists to end - so it is diagnostic, never the gate. """
```

Read it either side of the VEP subprocess in `dump_and_annotate_variants` and log the delta onto the run
(`set_task_log("cgroup_oom_kills", delta)`). Classification never depends on it, so nothing regresses on a
deployment without the cap; it just makes "this was memory" provable in the run's own logs.

---

## Phase 2 — Size: `AnnotationRangeLock.vep_buffer_size`

### Where the override lives, and why not on `AnnotationRun`

`AnnotationRun.reset_for_retry` (`annotation/models/models.py:1239`) overwrites the whole row with a fresh
instance, restoring only `pk`, `annotation_range_lock`, `pipeline_type` and `created`:

```python
fresh = AnnotationRun(annotation_range_lock=self.annotation_range_lock, pipeline_type=self.pipeline_type)
fresh.pk = self.pk
fresh.created = self.created
fresh.save()
```

A `vep_buffer_size` column on `AnnotationRun` would be wiped by the very retry that sets it, and the run
would go out at 500 again — looping until the attempt cap tripped. The same applies to any "already
retried once" counter: `attempt_count` is reset too.

Carving preserved fields into `reset_for_retry` erodes the property that makes it safe — its docstring
states the whole-row overwrite exists so there is "no field-by-field reset to keep in sync as the model
grows". So the override goes on `AnnotationRangeLock` instead:

```python
# annotation/models/models.py - AnnotationRangeLock

class AnnotationRangeLock(models.Model):
    ...
    # #1710: per-range VEP --buffer_size, overriding settings.ANNOTATION_VEP_BUFFER_SIZE. Lives on the
    # lock rather than the run because it describes what is *in the range* (a heavy large-indel load
    # spanning many 1Mb regions), not one attempt at annotating it - and because AnnotationRun.
    # reset_for_retry rewrites the whole run row, so a run-level override would be erased by the retry
    # that set it. Set automatically when VEP's output comes back short, or by hand for a range already
    # known to be bad.
    vep_buffer_size = models.IntegerField(null=True)
```

It also gives the one-retry cap for free, with no counter: **a non-null `vep_buffer_size` means this range
has already had its shot.** Since the automatic step-down goes straight to the floor (below), there is
nowhere further to step, and an operator who set the value by hand has stated their intent — in both cases
a still-truncated run should fail visibly rather than loop.

### Threading it into the command

```python
# annotation/vep_annotation.py

def get_vep_command(vcf_filename, output_filename, genome_build: GenomeBuild, annotation_consortium,
                    pipeline_type: VariantAnnotationPipelineType, compress_output: bool = True,
                    variant_annotation_version=None, buffer_size: Optional[int] = None):
    ...
    if buffer_size := buffer_size or settings.ANNOTATION_VEP_BUFFER_SIZE.get(pipeline_type):
        cmd.extend(["--buffer_size", str(buffer_size)])
```

with the resolution itself on the run, so every caller asks the same question:

```python
# annotation/models/models.py - AnnotationRun

def get_vep_buffer_size(self) -> Optional[int]:
    """ #1710: this run's --buffer_size - the range's override if it carries one, else the pipeline-type
        default from settings. Phase 4 additionally lowers it from the load measured at dump time. """
    if self.annotation_range_lock and self.annotation_range_lock.vep_buffer_size:
        return self.annotation_range_lock.vep_buffer_size
    return settings.ANNOTATION_VEP_BUFFER_SIZE.get(self.pipeline_type)
```

Callsites:

- `dump_and_annotate_variants` (`annotation/tasks/annotate_variants.py:490`) — passes
  `buffer_size=annotation_run.get_vep_buffer_size()`. `pipeline_command` is already persisted on the run,
  so the size used is visible per attempt with no extra field.
- `run_vep` (`annotation/vep_annotation.py:289`) — passes its new kwarg through.
- `annotation/views.py:335` (the VEP command shown on the settings page) and
  `build_snakemake_config` (`annotation/external_annotation.py:405`) both build a *template* with no run in
  hand, and keep the settings default. External runs (#1568) are operator-managed; the operator edits
  `vep_buffer_size` in `config.yaml`.

---

## Phase 3 — Retry: step to the floor, once

### Settings

```python
# variantgrid/settings/components/annotation_settings.py

# #1710: buffer size a run is retried at after VEP's output came back truncated. There is exactly one
# retry, so it goes straight to the floor rather than halving - 500 -> 250 only halves large-indels-per-
# buffer, which from a 5.3GB peak still lands near 3GB with the retry spent. The benchmark above records
# 100 as nearly free: runtime was flat across every size measured and output was byte-identical.
ANNOTATION_VEP_BUFFER_SIZE_TRUNCATION_RETRY = 100
```

### Detection point A — the import lane (the gate)

`check_all_dumped_records_accounted_for` and the cyvcf2 read guard raise `TruncatedVEPOutputError` inside
`import_vcf_annotations`, which `import_annotation_run` catches, records on `error_exception` and re-raises
→ status ERROR. The dispatcher only picks up CREATED / ANNOTATION_COMPLETED runs, so today an ERROR run
sits there until someone retries it by hand. There is no automatic path to hook into — the retry has to be
issued explicitly.

```python
# annotation/tasks/annotate_variants.py - import_annotation_run

    except Exception as e:
        ...
        # #1710: VEP's output was short. Record the range as needing a smaller buffer, so the retry
        # scheduled below re-runs it at the floor instead of reproducing the same failure. A range that
        # already carries an override has had its one retry - leave it ERROR and visible.
        retry_truncated = isinstance(e, TruncatedVEPOutputError) and _lower_buffer_size_for_range(annotation_run)
        ...
    finally:
        ...
        if retry_truncated:
            # After this task's own row updates have landed, so the async reset can't race them.
            annotation_run_retry(annotation_run, upload_only=False)
```

with the ordering deliberate: `annotation_run_retry` fires `reset_annotation_run_for_retry` via
`apply_async`, and that task rewrites the whole row. Calling it from the `except` block would race this
task's own `finally`, whose `save_if_owner(my_task_id)` would then be writing over a row the reset had
already replaced. Issuing it at the end of the `finally` — after the lease is cleared and ownership is
settled — leaves the reset with the row to itself.

The step-down itself is one small helper, and is the only place the one-retry cap is enforced:

```python
# annotation/tasks/annotate_variants.py

def _lower_buffer_size_for_range(annotation_run) -> bool:
    """ #1710: pin this range to the retry buffer size, returning whether it was this call that did so.

        False means the range already carries an override - either a previous truncation's retry or an
        operator's manual setting - so it has had its one attempt at a smaller buffer and should fail out
        visibly rather than cycle the box again. """
```

`annotation_run_retry(upload_only=False)` is the right existing entry point: a truncated VCF has to be
re-dumped and re-annotated, not re-imported, and it already routes through `reset_annotation_run_for_retry`
→ dispatcher → VEP lane (#1649). Because the override lives on the lock, the reset leaves it intact and the
next attempt reads it back through `get_vep_buffer_size()`.

### Detection point B — the VEP stage (fast path)

Where `perl` itself is the victim, the run fails at `return_code != 0` in `dump_and_annotate_variants` and
never reaches the import lane. Same treatment, gated on the heartbeat not being the killer:

```python
# annotation/tasks/annotate_variants.py - AnnotationRunLeaseHeartbeat

    def __init__(self, annotation_run, task_id):
        ...
        # #1710: _abort_process kills VEP with SIGKILL, which is indistinguishable from an OOM kill by
        # return code. Record that it was us, so a lost lease isn't misread as a memory failure and
        # doesn't shrink the buffer for a range that was never memory-hungry.
        self.aborted = threading.Event()
```

set in `_abort_process` alongside the existing `process.kill()`, and read in the `return_code != 0` branch:

```python
        if return_code != 0:
            ...
            if return_code < 0 and not (lease_heartbeat and lease_heartbeat.aborted.is_set()):
                # Killed by a signal we didn't send - VEP was the OOM victim (#1710). Same treatment as a
                # truncated output: the range is pinned to the retry buffer size and redispatched once.
```

This is an optimisation, not a second mechanism — it reaches the same `_lower_buffer_size_for_range` and the
same `annotation_run_retry`, called from `annotate_variants`'s `finally` for the same ordering reason. If it
were dropped entirely, the run would still be retried, one wasted import attempt later.

### Interaction with the attempt cap

`ANNOTATION_MAX_RUN_ATTEMPTS` counts *lease reclaims* (dead workers), bumped in `_lease_and_launch_run`, and
is untouched here. `reset_for_retry` sets `attempt_count` back to 0, so a truncation retry starts its lease
budget fresh — correct, because it is a new attempt at a differently-parameterised run, and the one-retry
cap that bounds it is the non-null `vep_buffer_size` on the lock.

---

## Phase 4 — Predict: choose the buffer size at dump time

The reactive path costs a full wasted VEP run. Where the load is knowable in advance it is strictly better
to size up front, and it is knowable: `dump_and_annotate_variants` dumps *before* it builds the VEP command,
and `write_qs_to_vcf` already selects `locus__position`, `locus__ref__seq`, `alt__seq` and `end`.

Calibrate this phase against the bisection in "What actually costs the memory" above, not against the SV
benchmark. The metric below is the right one if bins-per-buffer is confirmed as the driver; if the
transcript-count term dominates instead, `--transcript_filter` is the lever and this phase is worth less
than it looks.

### The metric

Peak RSS tracks **distinct 1Mb regions spanned per buffer**, not variant count — VEP holds the transcript
cache for every 1Mb region the current buffer spans and drops the rest in `AnnotationSource::clean_cache`.
The benchmark already in `ANNOTATION_VEP_BUFFER_SIZE` is a direct calibration of this: 1000 DELs of median
49 kb gave worst-case regions-per-buffer of 247 / 126 / 67 / 29 at buffer sizes 1000 / 500 / 250 / 100, for
1774 / 1212 / 752 / 505 MB — close to linear at roughly 5–6 MB per region over a ~375 MB floor.

It also explains why a heavy large-indel load blows up a *standard* run: a long explicit delins that dodged
symbolic conversion (#1708) spans many regions while counting as one variant against the buffer.

Because variants in a range lock are contiguous by ID and dumped position-sorted, regions accumulate
monotonically through the dump, so the count is a running total costing one integer:

```python
# snpdb/variants_to_vcf.py - write_qs_to_vcf

def write_qs_to_vcf(vcf_filename, genome_build, qs, info_dict=VARIANT_GRID_INFO_DICT,
                    use_accession=False, values_observer=None) -> int:
    ...
    sorted_values = qs.values("id", chrom_key, "locus__position",
                              "locus__ref__seq", "alt__seq", "end", "svlen")
    if values_observer:
        # Callers that want to measure the dump as it's written wrap the value stream here, rather than
        # re-querying it - #1710 counts the 1Mb regions it spans to size VEP's buffer.
        sorted_values = values_observer(sorted_values)
```

with the counter itself in `annotation/` (a generator wrapper tracking `(chrom, position // 1_000_000)`
through to `end // 1_000_000`, resetting on contig change), and its result stored for diagnosis:

```python
# annotation/models/models.py - AnnotationRun

    # #1710: distinct 1Mb regions the dump spans - what VEP's peak RSS tracks. Recomputed per dump, so
    # reset_for_retry clearing it is correct.
    dump_regions_spanned = models.IntegerField(null=True)
```

### The choice

```python
# variantgrid/settings/components/annotation_settings.py

# #1710: the memory budget one VEP buffer may span, in distinct 1Mb regions of transcript cache. The
# benchmark in ANNOTATION_VEP_BUFFER_SIZE puts peak RSS near 375MB + ~5.5MB per region, so 120 is
# roughly a 1GB ceiling per concurrent VEP. Confirm against a re-dump of a known-heavy range before
# trusting the constant.
ANNOTATION_VEP_MAX_REGIONS_PER_BUFFER = 120
# Sizes the chooser may pick, largest first. Every value here is benchmarked in ANNOTATION_VEP_BUFFER_SIZE.
ANNOTATION_VEP_BUFFER_SIZE_LADDER = [500, 250, 100]
```

After `dump_variants`, pick the largest rung whose estimated regions-per-buffer
(`dump_regions_spanned * rung / dump_count`) stays inside the budget, and take the lower of that and the
pipeline-type default — so the chooser only ever *reduces*, never raises a deployment's configured size.
An explicit `AnnotationRangeLock.vep_buffer_size` still wins over both: it is the operator's or the retry's
last word.

Phase 4 does not remove Phase 3. Prediction handles the region-density load; the reactive path still covers
what density does not explain (transcript-dense loci such as the BRCA1 blowout in #1228, a box under
concurrent pressure) and any range whose first dump predates this code.

---

## Model / migration summary

```python
# annotation/migrations/0166_vep_buffer_size.py
operations = [
    migrations.AddField('AnnotationRangeLock', 'vep_buffer_size', models.IntegerField(null=True)),
    migrations.AddField('AnnotationRun', 'dump_regions_spanned', models.IntegerField(null=True)),
]
```

Both nullable with no backfill: null means "no override" and "not measured", which is exactly the state of
every existing row.

---

## Operator surface

- **Grid** (`annotation/grids.py`, `AnnotationRunColumns`) — add
  `RichColumn(key="annotation_range_lock__vep_buffer_size", label="VEP Buffer", orderable=True,
  client_renderer='TableFormat.number')`, so an overridden range is visible without opening the run. The
  size actually used by a given attempt is already in `pipeline_command`.
- **Manual override** — a `set_annotation_run_buffer_size` POST view in `annotation/views.py` beside
  `subdivide_annotation_run` (`:600`), routed in `annotation/urls.py`, with its form in
  `annotation/templates/annotation/view_annotation_run.html` next to the existing subdivide form. This is
  the "park run 34818 by hand" lever the ticket asks for, and the reason the value is a size rather than a
  boolean. Superuser-only, matching subdivide.

`subdivide_annotation_range_lock` stays reachable only from that same manual action — splitting the range
leaves large-indel density per buffer unchanged, so it is not a memory fix and does not belong on the
automatic path.

---

## Tests

- `test_get_vep_buffer_size_override` — a lock carrying `vep_buffer_size` wins over the settings default;
  a null one falls back to it; the value reaches `get_vep_command`'s `--buffer_size` argument.
- `test_reset_for_retry_preserves_buffer_size` — the whole point of putting the field on the lock:
  `reset_for_retry` returns the run to CREATED with the range's override intact.
- `test_truncated_import_lowers_buffer_and_retries` — a `TruncatedVEPOutputError` from
  `import_vcf_annotations` sets the lock to `ANNOTATION_VEP_BUFFER_SIZE_TRUNCATION_RETRY` and schedules a
  full (not upload-only) retry.
- `test_truncated_import_second_time_fails_out` — the same failure against a lock that already carries an
  override leaves the run ERROR and schedules nothing.
- `test_lease_abort_is_not_treated_as_oom` — a `-9` from the heartbeat's own `_abort_process` leaves the
  lock's buffer size untouched (the regression the `aborted` event exists to prevent).
- `test_dump_regions_spanned` — a dump of variants with known spans, including one crossing a 1Mb boundary
  and a contig change, records the expected region count.
- `test_buffer_size_chosen_from_regions` — a region-dense dump picks a lower rung off the ladder; a sparse
  one keeps the pipeline-type default; an explicit lock override beats both.

The systemd change is not unit-testable — verify on the box with `systemctl show -p MemoryMax -p
MemorySwapMax -p OOMPolicy celeryd_annotation_workers.service`, checking the percentage has resolved to the
bytes expected for that box's RAM, then force a cgroup OOM and confirm the pool survives it and
`memory.events` `oom_kill` increments.

---

## Implementation order

1. **Phase 1, on its own.** `MemoryMax=50%` + `MemorySwapMax=0` + `OOMPolicy=continue` on the unit. Nothing
   else depends on it, and it is what stops one bad run taking four good ones down.
2. **The benchmark** — 34818's dump at buffer 500 vs 100, wall clock and peak RSS, plus the flag bisection.
   This decides how much of the rest is needed, so it comes before any of it.
3. **If 100 is flat: lower the default and stop.** Steps 4–5 stay (the operator lever is cheap and useful);
   steps 6–7 are dropped. See "First, try the one-line version".
4. Migration + `AnnotationRangeLock.vep_buffer_size` + `AnnotationRun.get_vep_buffer_size()` + the
   `get_vep_command` kwarg and its callsites. No behaviour change while every lock is null.
5. Grid column + manual override view/form — makes step 4 usable by hand, and gives run 34818 its immediate
   mitigation.
6. Phase 3 detection and retry: `_lower_buffer_size_for_range`, the import-lane gate, the heartbeat
   `aborted` event and the VEP-stage fast path.
7. `get_cgroup_oom_kill_count` logging, once Phase 1 is deployed and the counter is known to move.
8. Phase 4: `values_observer`, the region counter, `dump_regions_spanned`, the ladder chooser, with
   `ANNOTATION_VEP_MAX_REGIONS_PER_BUFFER` calibrated from step 2 rather than from the SV benchmark. Only
   reached if step 2 shows a lowered default is *not* safe to apply globally.

## Decisions

1. **Override on `AnnotationRangeLock`, not `AnnotationRun`.** Survives `reset_for_retry` untouched,
   describes the range rather than an attempt, and doubles as the one-retry cap.
2. **Named for truncation, not OOM.** The reliable signal is a short output; memory is an inference.
   Attribution via the cgroup counter is diagnostic logging only, so nothing regresses on a deployment that
   hasn't set `MemoryMax` yet.
3. **One retry, straight to 100.** With a single shot, halving spends it for a ~40% cut. The existing
   benchmark records 100 as flat on runtime and byte-identical on output.
4. **`attempt_count` is left alone.** It counts dead-worker reclaims; the truncation cap is the non-null
   override, and the two bound different failure modes.
5. **The chooser only lowers.** A deployment's configured `ANNOTATION_VEP_BUFFER_SIZE` stays the ceiling;
   Phase 4 may pick a smaller rung, never a larger one.
6. **External runs (#1568) keep the settings default.** `build_snakemake_config` builds an operator-edited
   template with no run in hand; `vep_buffer_size` is already a key in `config.yaml`.
7. **`--max_sv_size` stays as it is.** It bounds explicit long variants by length, but the variants here
   were ~1 kb, so no plausible setting would have caught them — and it is one of the
   `VariantAnnotationVersion` `get_or_create` kwargs, so lowering it re-annotates the database. Revisit at
   the next version roll if the megabase case (#1708) turns out to matter separately.
8. **Buffer size is the right lever after all — possibly as a plain default change.** At 1 kb the
   per-variant costs are indistinguishable from a SNV's; what remains is distinct 1 Mb bins per buffer,
   driven by how sparse the range is, which is precisely what `--buffer_size` bounds. `clean_cache` keeping
   only the current buffer's regions, over position-sorted input, is why a smaller buffer does not cost
   re-loads. If the benchmark confirms that, lowering the default beats every per-range mechanism in this
   plan on both simplicity and coverage.
9. **Bisect before calibrating.** The transcript-count × allele-length term is not ruled out and is not
   addressed by buffer size. One afternoon of `/usr/bin/time -v` runs over 34818's own dump settles it,
   and nothing downstream should hardcode a constant until it has.
10. **The memory cap scales with RAM; the buffer budget does not.** `MemoryMax=` is a percentage because
   systemd resolves it for free. Phase 4's `ANNOTATION_VEP_MAX_REGIONS_PER_BUFFER` stays a plain constant,
   even though "read the cgroup's `memory.max`, divide by worker concurrency, convert to a regions budget
   via the linear fit" is derivable: each step of that chain is an inference, and a wrong one silently
   mis-sizes every run rather than failing. The constant is benchmarked, greppable and overridable per
   deployment; revisit it if the fit is confirmed across more than one box.
