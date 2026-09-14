# Sending finalised TSO 500 case reports to Mocha

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: in progress - Parts A-D landed 2026-09-14; the open questions below, and the TSO 500
case_template / case_fields, are outstanding

Design for [sapath#444](https://github.com/SACGF/variantgrid_sapath/issues/444) (call the Mocha API
for TSO 500 reports / JSON). The reports themselves come from
[#444](https://github.com/SACGF/variantgrid/issues/444) (multi-variant classification and reporting),
whose `CaseReport` is the thing being sent. This is SA Path only, so almost all of it lives in the
`variantgrid_sapath` repo; the public repo gets two hooks.

Mocha is SA Pathology Molecular Oncology's tracking system ([SACGF/mocha](https://github.com/SACGF/mocha)).
VG already talks to it one way - the hourly import of tumour details described in the sapath repo's
*docs/mocha_api.md* - and this adds the other direction: when a scientist finalises a TSO 500 case
report, VG tells Mocha where the PDF and JSON are, and Mocha (after human review on its side) uploads
them to Omico.

---

## What Mocha expects

`POST /api/v2/report-location/` (the `ReportLocationCreateView` in Mocha's *core/api/views.py*, client
and payload in Mocha's *scripts/mocha_api.py*):

```json
{
  "patient_code": "C16716",                  // Omico C-number
  "helix_number": "GM-25-0046812",           // GM-number = Omico PathologyRequest.external_lab_id
  "seq_run": "TSO_25_082",                   // optional, may be blank
  "report_filename": "/tau/.../C16716_report.pdf",
  "report_filename_sha256": "...",
  "json_filename": "/tau/.../C16716_report.json",
  "json_filename_sha256": "..."
}
```

Facts that shape the design:

- **Matching** is `patient_code` + `helix_number` → PathologyRequest → its single SequencingData. One
  match is `201` (`200` on an exact resend); nothing yet is `202` and the row is parked and reconciled
  later; an ambiguous match is `202` with "Needs attention". Mocha never rejects a well-formed post,
  so we never need to re-send on our own account.
- **Idempotency** is the whole tuple, hashes included. A rebuilt document with a new hash is a new
  row on Mocha's side, which is exactly what a VG "New version" is.
- **Mocha reads the files off its own disk** - the Omico upload executor does `open(report_filename)`.
  Both services run on frgeneseq05vg, but VG's `MEDIA_ROOT` is not somewhere the `mocha` service user
  can read, and the JSON is not a file at all (only `CaseReport.json_output`). So VG has to export the
  two documents to a shared directory before it posts.
- **The Omico push is currently disabled** in the Mocha view (a human decides when to send), so
  nothing on the Mocha side needs changing for this plan.
- The response echoes the row: `id`, `status` ("Matched" / "Pending" / "Needs attention"),
  `match_error`, `created`.

---

## Data

### sapath: `MochaReportLocation` (new)

One row per `CaseReport` that has been, or should be, told to Mocha. A new report version is a new
`CaseReport`, so it gets its own row; a resend after an error updates the same row.

```python
class MochaReportLocation(TimeStampedModel):
    case_report = models.OneToOneField(CaseReport, on_delete=models.CASCADE)

    # What was posted - kept as text (MochaSampleExtraction is Mocha's cache, never FK'd into)
    patient_code = models.TextField(blank=True, default="")
    helix_number = models.TextField(blank=True, default="")
    seq_run = models.TextField(blank=True, default="")
    report_filename = models.TextField(blank=True, default="")
    report_filename_sha256 = models.CharField(max_length=64, blank=True, default="")
    json_filename = models.TextField(blank=True, default="")
    json_filename_sha256 = models.CharField(max_length=64, blank=True, default="")

    # Our side of the delivery
    status = models.CharField(max_length=1, choices=MochaReportStatus.choices,
                              default=MochaReportStatus.QUEUED)
    error = models.TextField(blank=True, default="")   # why UNRESOLVED / ERROR, for the card
    attempt_count = models.IntegerField(default=0)
    sent = models.DateTimeField(null=True, blank=True)

    # Mocha's side, from the response
    mocha_id = models.IntegerField(null=True, blank=True)
    mocha_status = models.TextField(blank=True, default="")      # "Matched" / "Pending" / "Needs attention"
    mocha_match_error = models.TextField(blank=True, default="")
    response = models.JSONField(null=True, blank=True)

    class Meta:
        ordering = ["-modified"]


class MochaReportStatus(models.TextChoices):
    QUEUED = "Q", "Queued"          # signal fired, task not yet run
    UNRESOLVED = "U", "Unresolved"  # VG could not work out the identifiers - retry after the Mocha import
    SENT = "S", "Sent"              # Mocha accepted it; mocha_status says whether it matched
    ERROR = "E", "Error"            # HTTP / filesystem failure after retries
```

### sapath: `MochaSampleExtraction.helix_number` (promoted column)

`helix_number` is in `raw_data` today (the model docstring says to promote a field when something
needs it). It becomes a column with a data migration that fills it from `raw_data`, and joins
`API_FIELDS` so the import keeps it current.

```python
    helix_number = models.TextField(null=True, blank=True)
```

### sapath settings

```python
SAPATH_MOCHA_REPORT_DIR = "/tau/data/clinical/tso500/vg_reports"   # exported PDF + JSON, readable by the mocha user
SAPATH_MOCHA_REPORT_TEMPLATES = ["TSO 500"]   # ClassificationReportTemplate names whose reports go to Mocha
```

`SAPATH_MOCHA_API_URL` and `SAPATH_MOCHA_TOKEN` already exist. With no report dir, or a template not
in the list, finalising does nothing Mocha-related - the same "unconfigured means skip" rule the
import follows.

### Core: the delivery row shown on the Reports card

A small data holder in `library/`, in the style of `library/integration_status.py:IntegrationDetail`,
so the public template can show what a deployment-specific app did with a report without knowing
which app:

```python
@dataclass
class CaseReportDelivery:
    label: str                       # "Mocha"
    status: str                      # bootstrap contextual class - success / warning / danger / secondary
    text: str                        # "Sent, matched" / "Unresolved: no Mocha extraction for 25-245-16107"
    timestamp: Optional[datetime] = None
    action_url: Optional[str] = None     # a POST that retries / resends
    action_label: Optional[str] = None   # "Send to Mocha"
```

---

## Part A - core hooks (public repo)

Two additions to `classification/`, both generic.

**A1. A finalise signal.** In `classification/models/classification_report_models.py`, next to the
models the way `classification/models/classification.py` declares its signals:

```python
case_report_finalised_signal = django.dispatch.Signal()  # args: "case_report", "user"
```

Sent from `classification/report/case_report_builder.py:finalise_case_report` once, when the report
goes DRAFT → FINAL, after `stamp_report_onto_classifications` has run so a receiver sees the finished
state. Re-finalising a FINAL report (how a LIS id entered later reaches the records) does not fire it.

**A2. A deliveries signal for the Reports card.**

```python
case_report_deliveries_signal = django.dispatch.Signal()  # args: "case_reports"; returns dict[int, list[CaseReportDelivery]]
```

`analysis/views/views_classify_report.py:_classify_report_context` sends it with the card's
`case_reports` list and merges the responses into `deliveries_by_report_id`.
`analysis/templates/analysis/case_report_card.html` grows a "Delivery" column rendered from it: a
badge, the text, and the action button (which uses the card's existing `case-report-action` POST
handler, so a retry re-renders the card like Finalise does). With no receivers the column stays
empty, so Shariant and variantgrid.com are unchanged.

`scripts/vg map` after adding the signals, so `claude/maps/signals.md` lists them.

## Part B - the client and its mock (sapath)

*sapath/mocha_api.py* gains a `ReportLocation` dataclass (the same seven fields as Mocha's own
client), `send_report_location(report) -> dict` over a `_post` that mirrors `_get` (token header,
JSON body, `raise_for_status`, one-minute timeout).

`MockMochaAPI` lives in the tests (*sapath/tests/mock_mocha_api.py*), a test double with the same
interface - there is only ever the one Mocha, so production takes `MochaAPI.from_settings()` and
nothing selects between them. It holds `extractions: list[dict]` for `sample_extractions()`,
`matched: set[tuple[str, str]]` of (patient_code, helix_number) pairs Mocha would resolve, and
`report_locations: list[dict]` of every post. `send_report_location` reproduces the real semantics
so tests assert on our handling of them: an exact resend returns the existing row with
`created: False`; a pair in `matched` answers "Matched"; anything else "Pending" with a
`match_error` in Mocha's wording. The report task takes the client as a `mocha_api` kwarg
defaulting to `MochaAPI.from_settings()`, the way `import_mocha_sample_extractions` already does,
so a test calls the task function with the mock; the view and receiver paths that go through
`.delay()` patch `MochaAPI.from_settings` to return it.

## Part C - resolving a case to Mocha's identifiers (sapath)

*sapath/mocha_report.py* (new), entry point `resolve_report_identifiers(case_report) -> ReportIdentifiers`
(a dataclass of `patient_code`, `helix_number`, `seq_run`) that raises `UnresolvedReport(message)`
with the sentence the card shows.

1. **Specimen.** `CaseReport.source` by level: SPECIMEN is itself, EXTRACTION → `extraction.specimen`,
   SAMPLE → `sample.extraction.specimen`, PATIENT → the specimens the case's samples reach
   (`patients/sample_grouping.py:get_sample_group`) - one is fine, more than one is unresolved
   ("report covers specimens X and Y - build it from the specimen").
2. **Container accession.** `Specimen.reference_id` is the full Helix accession (`25-245-16107B`);
   `HelixAccession.split_accession_and_suffix` in *sapath/models/sapath_helix.py* gives the container.
3. **Mocha extraction.** `MochaSampleExtraction` on that container gives `patient_code` and the newly
   promoted `helix_number`. None yet is unresolved ("no Mocha extraction for 25-245-16107 - the
   hourly import may not have it yet"); several with differing identifiers is unresolved too.
   `Patient.patient_code` is checked against Mocha's when both are set, and a mismatch is unresolved
   rather than silently sending either.
4. **Run.** `Sample.sequencing_run` for each of the case's samples
   (`snpdb/models/models_vcf.py:Sample.sequencing_run`); when they all name one run its `name` is
   sent, otherwise blank - Mocha matches on the other two and the field is optional there.

## Part D - export, send, record (sapath)

*sapath/tasks/mocha_report_task.py* (new), `send_case_report_to_mocha(case_report_id)` on the
`web_workers` queue like the Helix task, with celery `autoretry_for=(requests.RequestException,)`,
`max_retries=5` and exponential backoff. What it does, in order:

1. Load the `CaseReport` and its `MochaReportLocation` (created QUEUED by the receiver), bump
   `attempt_count`.
2. Resolve identifiers (Part C). Unresolved → status UNRESOLVED with the message, no retry (a retry
   is a human clicking the button after the import catches up).
3. Export. Write `<SAPATH_MOCHA_REPORT_DIR>/<patient_code>_<helix_number>_vg<case_report.pk>.pdf` from
   `pdf_file` and `.json` from `json.dumps(json_output, indent=2)`; sha256 of the bytes written.
   The pk in the name keeps versions apart; the C- and GM-numbers keep it readable for the MO team.
4. Post through `mocha_api.send_report_location(...)`, inside
   `IntegrationActivity.track("sapath-mocha-report-location", name="Mocha report location",
   direction=OUTBOUND)` (`eventlog/models.py:IntegrationActivity`) so Server Status shows the last
   send and last failure without a run log of our own.
5. Record: the seven posted fields, `sent`, `mocha_id`, `mocha_status`, `mocha_match_error`, the
   whole response, status SENT. An exception after the retries → ERROR with the message.

*sapath/signals/mocha_report.py* (new) has the two receivers, connected in *sapath/apps.py*:

- `case_report_finalised_signal` → if the template is in `SAPATH_MOCHA_REPORT_TEMPLATES` and a
  report dir is configured, `get_or_create` the QUEUED row and `.delay()` the task.
- `case_report_deliveries_signal` → one `CaseReportDelivery` per report that has a row: SENT +
  "Matched" is `success`; SENT + "Pending" / "Needs attention" is `warning` with Mocha's
  `match_error`; UNRESOLVED and ERROR are `danger` with ours. Every non-QUEUED row gets the
  "Send to Mocha" action; a FINAL report of a Mocha template with no row at all (finalised before
  this landed) gets it too, so the backlog can be sent by hand.

A `sapath_send_case_report_to_mocha` POST view in *sapath/urls.py* is the action: checks
`case_report.check_can_write`, upserts the row to QUEUED and enqueues the task. The existing
`mocha_integration_status` receiver gains a second `IntegrationStatus` for the new key with the
row count and a link to the `MochaReportLocation` admin.

## Part E - Mocha side

Nothing for this plan. Two things to raise on the Mocha repo when this lands, since they are Mocha's
decisions: the `mocha` service user needs read access to `SAPATH_MOCHA_REPORT_DIR`, and whether to
whitelist that directory when the Omico push is re-enabled.

---

## Tests

sapath, all under `@override_settings(SAPATH_MOCHA_REPORT_DIR=<tmp>, SAPATH_MOCHA_REPORT_TEMPLATES=[...],
CELERY_TASK_ALWAYS_EAGER=True)` with a `MockMochaAPI` handed to the task, building a real report the way
`analysis/tests/test_case_report.py` does and the Helix / Mocha rows the way *sapath/tests/test_import_mocha.py* does:

- Resolution: a specimen case with DNA and RNA arms on one container resolves to the one
  extraction's identifiers; two runs → blank `seq_run`; no Mocha extraction → UNRESOLVED with the
  container in the message; a patient case spanning two specimens → UNRESOLVED.
- Export and send: the two files exist, the hashes match their bytes, the mock received the seven
  fields, and the row is SENT with Mocha's status - once "Matched", once "Pending" with
  `match_error` carried through.
- Gating: a template outside the list, or no report dir, leaves no row and no files.
- Resend: the button's view after UNRESOLVED re-queues and, with the extraction now present, ends
  SENT; a second identical send is `created: False` and does not change the row's `sent`.

core (`analysis/tests/test_case_report.py`): `case_report_finalised_signal` fires once on the first
finalise and not on the second; the card renders a delivery a test receiver returns.

Keep only the tests that cover a branch we wrote (the resolution rules, the status mapping, the
gating); the mock's own fidelity to Mocha is not something to test twice.

---

## Docs

- The sapath repo's *docs/mocha_api.md* gets an "Outbound: report locations" section - the settings,
  the directory and its permissions, the file naming, what each card status means and what to do
  about it.
- `classification/CLAUDE.md`: one line each for the two signals, next to the finalise notes.

## Open questions

1. **JSON shape.** Mocha passes `json_output` straight to Omico's `upload_sequencing_report_json`.
   The issue calls the current JSON "produced by Shashi's script"; nobody has yet compared that
   output against a VG report's JSON. *sapath/tso500_report.py* is where any difference gets fixed,
   but the comparison has to happen before Mocha re-enables the Omico push.
2. **Which template names.** `SAPATH_MOCHA_REPORT_TEMPLATES` needs the real name(s) of the TSO 500
   template on sapathtest / production.
3. **Report dir.** The path above is a placeholder; the MO team should say where they want the
   files, since Mocha's UI shows the paths to them.
