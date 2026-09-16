# #1179 / #1005 — VCF preprocess pipe: name the stage that failed and show its stderr

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: draft

[#1179](https://github.com/SACGF/variantgrid/issues/1179) (bcftools fails in the pipe without being caught) and
[#1005](https://github.com/SACGF/variantgrid/issues/1005) (better error handling when a preprocess tool dies). One plan
because the fix is the same code: `upload/vcf/vcf_preprocess.py:run_pipe`.

## The problem

Preprocess is one bash pipe (zcat, manage.py vcf_clean_and_filter, bcftools norm, bcftools view --no-header,
manage.py vcf_clean_alts, split) built by `upload/vcf/vcf_preprocess.py:_build_pipe_commands` and run by
`upload/vcf/vcf_preprocess.py:run_pipe`. Since `b27e3be11` (2026-03, pipefail) a stage that dies does fail the
pipeline, so the "looks like it is still running" half of #1179 is gone. What is still missing is *which* stage failed
and *what it said*:

- `run_pipe` gives bash a single `stderr=PIPE`, so the stderr of all six stages arrives as one string, and the
  `CalledProcessError` it raises carries the whole pipe as `cmd`. `upload/tasks/vcf/import_vcf_step_task.py:ImportVCFStepTask.run`
  formats that with `library/utils/os_utils.py:format_called_process_error` into the step's `error_message`, the
  pipeline's `progress_status` and the Rollbar message: "Command '<300 characters of pipe>' returned non-zero exit
  status 255" followed by every stage's stderr mixed together. On a real file the interesting line is buried under
  thousands of `[W::vcf_parse_format_dict2]` warnings.
- The sub steps (`upload/vcf/vcf_preprocess.py:create_sub_step`: clean/filter, sort, normalize) never get
  `output_text` - the code that fills it is the other branch of `run_pipe` (below), which nobody runs.
  `upload/models/models.py:UploadStep.close_sub_steps` then stamps them all with the parent's status, so on the
  pipeline page three sub steps show ERROR and none says why.

Verified on this box (bcftools 1.20) with the #1179 repro - `upload/test_data/vcf/grch38_brca1.vcf` with its eleven
`##FORMAT` lines removed, run through `bcftools norm --multiallelics=- | bcftools view --no-header` under pipefail:
exit 255, `PIPESTATUS` `0 255 0 0`. `bcftools norm` warns per undeclared FORMAT tag, then
`[E::vcf_format] Invalid BCF, the FORMAT tag id=31 at chr17:43074436 not present in the header` /
`[flush_buffer] Error: cannot write to -` and dies; `bcftools view` exits 0 with empty output (it had already read the
header norm wrote), so the "Failed to read from standard input: unknown file type" line in the issue is from an older
bcftools. The culprit is unambiguous from `PIPESTATUS`; today's message hides it.

### Drift from the issues

- **The `vt` bullets in #1005 are done in bcftools form.** `vt` was removed in `6230ec00f` (#127);
  `upload/vcf/vcf_preprocess.py:get_bcftools_tool_version` records the bcftools/htslib version on every sub step, and
  `variantgrid/deployment_validation/tool_version_checks.py:_check_bcftools_version` (`manage.py deployment_check`,
  "Tool versions") enforces `_REQUIRED_BCFTOOLS_VERSION` 1.20 with the install wiki link as the fix text, alongside
  `check_vcf_split_pipe` for the split stage. Nothing here adds a per-import version gate.
- **"Error: No split VCF records"** (#1005's symptom) no longer exists: removed in `320056301` (#1079). An empty
  post-clean VCF is now `upload/tasks/vcf/import_vcf_step_task.py:ImportVCFStepTask._handle_no_vcf_records` - a warning
  and the data stages skipped - reached from `upload/tasks/vcf/import_vcf_tasks.py:ScheduleMultiFileOutputTasksTask`.
  With the pipe failing on any non-zero stage, an empty split dir can only mean an empty VCF.
- **"Break the pipeline into bits so we can tell what part failed"**: the sub steps already are the bits; they just
  never get their stderr. Splitting the pipe into separately run stages would mean writing the intermediate VCFs to
  disk (a 50M-record file three times over) for no gain once `PIPESTATUS` names the stage. Kept as one pipe.
- **`settings.VCF_IMPORT_PREPROCESS_POPEN_SHELL`** has been `True` since `8d86ef9fe` (2024-05) and is assigned only in
  `variantgrid/settings/components/default_settings.py`. The `False` branch of `run_pipe` (one `Popen` per stage, stderr
  to files, `communicate()` in reverse) is dead, and it is the "looks like it is still running" hang: the parent never
  closes its copy of each stage's stdout read end, so a stage whose downstream has died never gets SIGPIPE and blocks
  forever on a full pipe. It goes.
- `claude/research/upload.md` ("Preprocess: one shell pipe") says the sub steps are "where `ToolVersion` and the
  per-stage stdout/stderr live" - only the version is true today. Fixed below.

## Data

No model changes. `UploadStep.output_text` / `error_message` / `tool_version` already exist on the sub steps; they are
finally written.

```python
@dataclass
class PipeStage:
    """ One stage of the preprocess pipe, as run """
    name: str                       # key in pipe_commands, e.g. UploadStep.NORMALIZE_SUB_STEP
    command: list[str]
    stderr_filename: str
    returncode: Optional[int] = None    # None until the pipe has run
    stderr: str = ""                    # truncated to MAX_STDERR_OUTPUT, head and tail

class PipeStageError(CalledProcessError):
    """ CalledProcessError for the stage that caused the failure (cmd/returncode/stderr are that stage's),
        with the whole pipe's results for the message """
    stage: PipeStage
    stages: list[PipeStage]
```

`PipeStageError` subclasses `CalledProcessError` so nothing upstream changes: `ImportVCFStepTask.run` still catches
it and formats it, and `preprocess_vcf`'s unsorted-marker retry still catches it.

## Code changes

### `upload/vcf/vcf_preprocess.py`

- `_build_pipe_commands`: every stage gets a sub step, not just the three that carry the bcftools version -
  `zcat`/`cat`, `vcf_clean_and_filter`, `sort_vcf`, `normalize`, `remove_header`, `vcf_clean_alts`, `split_vcf`, in
  pipe order. `tool_version` is the bcftools `ToolVersion` for the sort/normalize/remove_header stages and `None` for
  the rest (`upload/vcf/gene_level_vcf_preprocess.py` already does this for its split stage "so a failed split reports
  on the pipeline the way the bcftools stages do"). A truncated gzip (`zcat: unexpected end of file`) is a real failure
  mode that today reads as "the pipe failed", so the reader stage is included.
- `run_pipe(pipe_commands, sub_steps, split_env, upload_pipeline)` - one implementation, bash. The shell string is
  built as today (stages joined with `' | '`, arguments joined with spaces, since `SPLIT_VCF_FILTER` is pre-quoted
  shell) with two additions:
  - each stage gets `2>` + `shlex.quote(stderr_filename)`, where `stderr_filename` is
    `get_import_processing_filename(upload_pipeline.pk, f"stderr_out_{stage_name}.log")` (the name the dead branch
    used; a sort retry overwrites them);
  - after the pipe, `; echo "${PIPESTATUS[*]}" > <quoted status file>` (pipestatus.txt in the same dir). `PIPESTATUS`
    is expanded before the `echo` runs, so it is still the pipe's. The status file, not `$?` or `set -o pipefail`,
    is the failure signal: bash's own exit status is the `echo`'s.
  - `Popen(..., shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE, env=split_env)` as now; bash's own
    stderr (a stage binary missing, a syntax error) stays on that pipe. If the status file is missing afterwards, raise
    a plain `CalledProcessError(p.returncode, shell_command, stderr=bash_stderr)` - the pipe never ran.
  - Parse the status file into one `PipeStage` per stage, read each stderr file through the existing
    `MAX_STDERR_OUTPUT` head+tail truncation (the `[E::...]` line is the last thing bcftools writes, so the tail keeps
    it under any number of warnings), and write `sub_step.output_text` for every stage. Log each non-empty stderr at
    warning level as before.
  - If any returncode is non-zero: pick the culprit with `_failed_stage(stages)`, set that sub step's
    `error_message` to its stderr, and raise `PipeStageError`.
- `_failed_stage(stages) -> PipeStage`: the first stage in pipe order whose exit is non-zero and is not a
  consequence of a stage after it dying: `returncode == 141` (128 + SIGPIPE - `zcat`, `cat`, bcftools) or
  `"BrokenPipeError" in stderr` (Python ignores SIGPIPE, so `manage.py vcf_clean_and_filter` exits 1 with a traceback
  when `bcftools norm` goes away under it). If every non-zero stage is SIGPIPE-shaped, the last one. Docstring says
  why; it is the one piece of logic here that is easy to get wrong later.
- `PipeStageError.__str__`: `Preprocess stage 'normalize' (3 of 6) failed: ` + `CalledProcessError.__str__` (which
  names the stage's command, not the pipe), then one line `exit codes: cat=0 vcf_clean_and_filter=0 normalize=255
  remove_header=0 vcf_clean_alts=0 split_vcf=0`. `format_called_process_error` appends `stderr:` with the stage's
  stderr, so the step `error_message`, `progress_status` and the Rollbar message all read: stage name, command, exit
  code, every stage's code, and the stage's own stderr - nothing else's.
- Delete the `Popen`-per-stage branch. `_reset_for_retry` is unchanged: the stderr and status files are overwritten
  by the next run and the split dir is already cleared.

### `variantgrid/settings/components/default_settings.py`

Remove `VCF_IMPORT_PREPROCESS_POPEN_SHELL`. `vg settings` shows no other file assigns it.

### `upload/vcf/gene_level_vcf_preprocess.py`

Unchanged in behaviour - it already gives its split stage a sub step. Add sub steps for its `cat` and `sed` stages
with `tool_version=None` so the three-stage pipe reports the same way.

### Docs

- `upload/CLAUDE.md`, under the "Preprocess is one shell pipe" pattern: `run_pipe` reads `PIPESTATUS` and a stderr
  file per stage; the failing stage is the first non-zero exit that is not SIGPIPE-shaped (141, or a Python
  `BrokenPipeError`) - those are stages killed because the one after them died. Every stage is a sub step and its
  stderr is its `output_text`.
- `claude/research/upload.md`, "Preprocess: one shell pipe": replace the pipefail sentence with the status-file
  mechanism and make the "per-stage stderr" claim true.
- `scripts/vg map` after the setting removal.

## Tests

`upload/tests/vcf/test_vcf_preprocess.py` already builds the `UploadPipeline`/`UploadStep` and patches
`get_bcftools_tool_version`; its retry tests patch `run_pipe` with `CalledProcessError` side effects and keep passing
(`PipeStageError` is one). `test_reset_for_retry_clears_sub_steps_split_files_and_stats` counts `len(sub_steps)` and
follows the new count.

New class `TestRunPipe` in the same module, running the real `run_pipe` on toy stages - no bcftools, no fasta, the
pipe is `cat` / `sh -c` / `yes` / `python3`:

1. **Failing stage is named, with its own stderr and code.** `{"a": ["cat", vcf], "b": ["sh", "-c", "'echo boom >&2;
   exit 3'"], "c": ["cat"]}` (the `sh -c` body pre-quoted, as `SPLIT_VCF_FILTER` is): raises `PipeStageError` whose
   `stage.name == "b"`, `returncode == 3`, `stderr == "boom"`, `cmd` is stage b's command; `str(e)` contains
   `'b' (2 of 3)` and `a=0 b=3 c=0`; sub step b's `error_message` is `boom`, sub step a's `output_text` is empty.
2. **Upstream SIGPIPE and BrokenPipeError are not blamed.** `yes | sh -c 'exit 3' | cat` gives `141 3 0`;
   `python3 -c "print('x\n' * 10**7)" | sh -c 'exit 3' | cat` gives a Python traceback ending in `BrokenPipeError`
   on stage 1. Both name stage 2. This is the `_failed_stage` rule.
3. **Success writes every sub step's stderr and raises nothing.** `cat vcf | sh -c 'echo note >&2; cat' | cat`:
   no exception, sub step 2 `output_text == "note"`, status file consumed.
4. **The #1179 repro**, `@skipUnless(shutil.which(settings.BCFTOOLS_COMMAND))`: new test data
   upload/test_data/vcf/no_format_header.grch38_brca1.vcf (`upload/test_data/vcf/grch38_brca1.vcf` minus its
   `##FORMAT` lines) through `cat | bcftools norm --multiallelics=- - | bcftools view --no-header - | cat` -
   `bcftools norm` splits multi-allelics without `--fasta-ref`, so no reference is needed. Asserts the culprit is the
   norm stage and its stderr contains `[E::vcf_format]`. Skipped on a box without bcftools; CI has it.

Tests 1-3 cover the logic we write (attribution, the SIGPIPE rule, stderr routing); 4 is the issue's own repro run
against the real tool, and is the one to drop if a bcftools release changes the message. `_build_pipe_commands`
gaining sub steps for every stage is covered by the existing count test.

## Manual verification

On this box (bcftools 1.20, GRCh38 reference present):

1. `python3 manage.py deployment_check` - "Tool versions" still all valid (`check_vcf_split_pipe` is untouched).
2. Ask first, then `python3 manage.py import_vcf upload/test_data/vcf/no_format_header.grch38_brca1.vcf --name
   no_format_header --user <you>` (`upload/management/commands/import_vcf.py:Command`; the pipe fails before any
   variant insert). The pipeline page shows status ERROR; `progress_status` starts
   `Preprocess stage 'normalize' (3 of 6) failed: Command 'bcftools norm ...' returned non-zero exit status 255`,
   the `exit codes:` line, then the `[E::vcf_format]` line. Expanding the `normalize` sub step row shows the same
   in Error Message and the warnings in Output Text; the `vcf_clean_and_filter` sub step's Output Text is its own
   (empty, or a `BrokenPipeError` traceback on a larger file). Rollbar has the same message.
3. Same import with `upload/test_data/vcf/grch38_brca1.vcf`: SUCCESS, every sub step has `output_text` (norm's
   `Lines total/split/realigned/skipped` summary), no `error_message`.
4. Retry the failed pipeline from its page: the stderr files are overwritten, the message is the same.
5. `ls <IMPORT_PROCESSING_DIR>/<pipeline pk>/` shows `stderr_out_<stage>.log` per stage and pipestatus.txt;
   `IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS` removes them with everything else.

## Decisions made

- **Keep one bash pipe and read `PIPESTATUS`** rather than one `Popen` per stage. The per-stage branch already
  existed, was never used, and hangs on exactly the failure being fixed; a status file is a dozen lines and keeps
  `split --filter` and the sort retry as they are.
- **The status file replaces pipefail as the failure signal.** With `PIPESTATUS` in hand every stage's code is known;
  keeping `set -o pipefail` as well would make bash's exit the pipe's and hide whether the `echo` ran. One mechanism.
- **Every stage is a sub step.** The alternative - only bcftools stages, logging the rest - leaves a `zcat` failure
  or a `vcf_clean_alts` traceback out of the page. Four more rows per import is the cost.
- **`close_sub_steps` still stamps all sub steps with the parent's status.** The failing one carries
  `error_message`; the others carry `output_text`. Making the successful stages show SUCCESS on a failed pipeline is
  a nicety that touches `ImportVCFStepTask.run` for no diagnostic gain.
- **Culprit = first non-zero non-SIGPIPE-shaped exit**, and the full code table is in the message regardless, so a
  wrong guess costs a reader one line, not the information.
- **No per-import bcftools version check.** `deployment_check` enforces 1.20 with the wiki link; a stage failing on an
  old bcftools now names the stage and the tool's own message, which is what the wiki bullet in #1005 wanted.

## Definition of done

- `scripts/vg tests --explain` names `upload.tests.vcf.test_vcf_preprocess` and it passes with `--keepdb`; the kept
  tests are the four above.
- `vcf_preprocess.py` module docstring states the status-file mechanism and the culprit rule; `upload/CLAUDE.md`
  gotcha line added; `claude/research/upload.md` sentence corrected; `scripts/vg docs check` passes.
- `scripts/vg map` refreshed (setting removed).
- This plan's `Status:` updated.
