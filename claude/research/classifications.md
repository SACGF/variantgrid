# classification — research notes

Verified against 7c4408c62 on 2026-09-06

The classification app is a lab's record of what a variant means for a condition, and everything Shariant built around
sharing those records between labs: versioned evidence, share levels, allele resolution, per-allele grouping, discordance
detection and resolution, condition matching to ontology terms, and outbound ClinVar submission. This is the long story
behind the rules in `classification/CLAUDE.md`; the vocabulary is in `claude/domain.md#classification-classification`,
and the model, URL, task, signal and command inventories are the generated maps (`claude/maps/models.md#classification`,
`claude/maps/urls.md#classification`, `claude/maps/tasks.md#classification`, `claude/maps/signals.md#first-party`,
`claude/maps/commands.md`). The two readmes, `classification/__classification_readme.md` and
`classification/__discordance_readme.md`, are the original author's summary.

## Flows

### A record arrives

Every write, whatever the door, ends in `classification/models/classification_inserter.py:BulkClassificationInserter.insert`.
The DRF endpoint `classification/views/classification_view.py:ClassificationView.post` feeds it one record (the web form)
or a `records` list with an `import_id` (OmniImporter and `sync`); the same view is mounted at v1, v2 and v3 paths
(`claude/maps/urls.md#classification`), the version only changing the response shape via `api_version`. A file upload
from a lab's curation system lands as `classification/models/uploaded_classifications_unmapped.py:UploadedClassificationsUnmapped`
and `classification/tasks/classification_import_map_and_insert_task.py:ClassificationImportMapInsertTask` shells out to
the external OmniImporter (`settings.CLASSIFICATION_OMNI_IMPORTER_APP_DIR`) to map it to VariantGrid JSON before calling
the same inserter. `csv_classification_inserter` (`claude/maps/commands.md`, #1481) is the command-line
door. Each record is addressed by `classification/models/classification_ref.py:ClassificationRef.init_from_str` -
`org/lab/lab_record_id.version` - and `lab_record_id` is unique per lab, so a re-import of the same file updates rather
than duplicates.

`insert` pops the control keys (`publish`/`share`, `delete`, `delete_reason`, `source`, `editable`, `return_data`) off
the payload and treats what is left as evidence. A new record goes through
`classification/models/classification.py:Classification.create_with_response`; an existing one through
`classification/models/classification.py:Classification.patch_value`. Records that arrive with `source: api` and no
`editable` flag have every patched field marked immutable at API level, which is what stops a curator editing a value in
the web form that the lab's system will overwrite on the next sync (`classification/enums/classification_enums.py:SubmissionSource.can_edit`).
A bulk import wraps the rows in a `classification/models/classification_import_run.py:ClassificationImportRun` keyed by
`username#import_id`, whose `ONGOING` status is the guard the rest of the app checks before doing per-record work.

### Evidence is patched, not assigned

`Classification.patch_value` is the only mutator. It normalises the incoming dict into a
`classification/models/evidence_key.py:VCDataDict` of `VCDataCell`s, uppercases the gene symbol, strips whitespace from
the c.HGVS, sends `classification/models/classification.py:classification_validation_signal` (receivers in
`classification/models/classification_variant_fields_validation.py:validate_variant_fields` and friends attach
per-key validation messages), then diffs each cell against the current evidence and keeps only what changed. A cell whose
existing immutability outranks the submission source is dropped with an `immutable` warning rather than raising. What
survives becomes the `delta` of a new `classification/models/classification.py:ClassificationModification`, unless
`classification/models/classification.py:ClassificationModification.is_edit_appendable` says the last unpublished
modification is by the same user and source within a minute, in which case the delta is merged into it - this is why a
curator typing in the form for a minute makes one version, not thirty. The patch also refreshes the denormalised
`clinical_significance` and `allele_origin_bucket` columns, and on a first save with `requires_auto_population` set,
fills annotation-derived keys through `classification/autopopulate_evidence_keys/autopopulate_evidence_keys.py:classification_auto_populate_fields`.

### Resolving the variant

A classification never links to a Variant directly. `classification/models/classification.py:Classification.ensure_allele_info_with_created`
calls `classification/models/classification_variant_info_models.py:ImportedAlleleInfo.get_or_create` with exactly what
the lab sent (c.HGVS or g.HGVS, transcript, genome build patch version), unique on an md5 of that text because Postgres
cannot index a 3 kb HGVS (#753). The first time an ImportedAlleleInfo is created it derives a `VariantCoordinate` from the
HGVS (or, since #1506, recognises a `BCR::ABL1` gene pair as a gene-level fusion) and records validation; then
`classification/models/variant_resolver.py:VariantResolver.queue_resolve` attaches it to a per-build
`classification/models/classification.py:ClassificationImport` and, when the inserter finishes or 100 are queued, fires
`classification/tasks/classification_import_task.py:process_classification_import_task`.

That task runs `classification/classification_import.py:process_classification_import`: known coordinates are matched in
bulk through `VariantPKLookup`; unknown ones are written to a synthetic VCF and pushed through the ordinary upload
pipeline (`classification/classification_import.py:_classification_upload_pipeline`, gene-level coordinates on their own
pipeline) so that variant insertion stays behind the single `variant_id_single_worker`. After insertion
`classification/tasks/classification_import_process_variants_task.py:ClassificationImportProcessVariantsTask` links the
new Variants back, populates ClinGen allele ids and schedules liftover (`settings.LIFTOVER_CLASSIFICATIONS`). Each step
that produces a Variant ends in `classification/models/classification_variant_info_models.py:ImportedAlleleInfo.set_variant_and_save`,
which fills the per-build `classification/models/classification_variant_info_models.py:ResolvedVariantInfo` caches
(c.HGVS per build, transcript version, gene symbol - what the grids sort on), recomputes validation, sets the status and
sends `classification/models/classification_variant_info_models.py:allele_info_changed_signal`. Receivers of that signal
copy the allele onto every classification sharing the ImportedAlleleInfo, reassign groupings and recalculate clinical
contexts. `classification/signals/classification_liftover.py:liftover_run_complete_handler` does the last pass with
`force_complete=True`, so a build still missing after liftover is recorded as unattainable rather than pending (#1420).

Validation is a versioned row, `classification/models/classification_variant_info_models.py:ImportedAlleleInfoValidation`,
with tags for normalisation diffs, liftover diffs, missing builds and unsupported transcripts; any `E` severity sets
`include=False` and `classification/models/classification.py:Classification.include_based_on_allele_info` keeps the
record out of exports until a human confirms it (`confirmed_by`) or the match is redone.
`ImportedAlleleInfo.hgvs_converter_version` records which resolver and cdot data produced the match (#1321), so a later
HGVS library upgrade can be audited with `classification/models/classification_variant_info_models.py:ImportedAlleleInfo.dirty_check`.

### Publish and share

Nothing outside the owning lab sees an unpublished modification. `classification/models/classification.py:Classification.publish_latest`
hands the last edited modification to `classification/models/classification.py:ClassificationModification.publish`, which
snapshots `published_evidence`, flips `is_last_published` from the previous version, grants the ShareLevel's Guardian
group read permission on the modification, rewrites `Classification.summary` through
`classification/models/evidence_mixin_summary_cache.py:ClassificationSummaryCalculator`, and sends
`classification/models/classification.py:classification_post_publish_signal`. That one signal drives most of the app
(the receiver list is in `claude/maps/signals.md#first-party`): the submitted flag and clinical-context recalculation,
grouping assignment, condition text sync, ClinVar exclusion patterns, the significance-change flag, common-variant
partition moves in snpdb, and gene-count refresh in annotation.

Share levels are `classification/enums/classification_enums.py:ShareLevel` - user, lab, organisation, logged-in users,
public - and only the last two are `is_discordant_level`, i.e. count as "shared" for discordance, grouping visibility
and export. The web form republishes at the record's current level (`classification/views/views.py:create_classification_object`)
and the inserter refuses to go lower, republishing at the current level with a `shared_higher` warning, but the model has
no constraint: the ratchet is convention. `user` is never a publish target - the inserter maps it to "do not publish".
Per-key visibility is separate: `classification/models/evidence_key.py:EvidenceKey.max_share_level` and
`classification/models/classification.py:Classification.get_visible_evidence` blank out keys such as patient identifiers
for anyone whose `lowest_share_level` is below the key's ceiling, even on a public record.

### Withdrawal and deletion

`classification/models/classification.py:Classification.set_withdrawn` is the soft delete: it sets `withdrawn` with a
`classification/enums/classification_enums.py:WithdrawReason`, opens the withdrawn flag and sends
`classification/models/classification.py:classification_withdraw_signal`, which removes the record from its grouping
and recalculates its clinical context. The inserter's `delete: true` withdraws a shared record and hard-deletes an unshared
one (unless `settings.CLASSIFICATION_ALLOW_DELETE` is off, as it is on Shariant); `delete: "withdraw"` always withdraws;
`delete: false` un-withdraws. Since #1481 a re-imported withdrawn record stays withdrawn and returns a `withdrawn` warning
instead of silently resurrecting.

### Grouping

`classification/models/classification_grouping.py:ClassificationGrouping` is the per-lab, per-allele-origin-bucket,
per-share-level roll-up that the listing grids read, under an `AlleleOriginGrouping` and `AlleleGrouping` per allele.
`classification/models/classification_grouping.py:ClassificationGrouping.assign_grouping_for_classification` runs on
publish, withdraw, allele change and condition change; it moves the record's `ClassificationGroupingEntry` and marks the
old and new groupings dirty. `classification/models/classification_grouping.py:ClassificationGrouping.update` rebuilds the
cached latest modification, condition terms, zygosities and pathogenic/somatic difference; it runs immediately for a
grouping's first record and otherwise waits for `classification/signals/classification_hooks_grouping.py:_instant_undirty_check`,
which does nothing while an import run is ongoing and lets `classification_imports_complete` do one sweep instead.
`classification/models/classification_grouping.py:ClassificationGroupingSearchTerm` rows (gene symbol, condition terms,
SCV) are fed by `classification_grouping_search_term_signal` so the grid's search never touches evidence JSON.
`settings.CLASSIFICATION_NEW_GROUPING` (on for Shariant) switches the listing pages to these groupings.

### Discordance

Publishing at a discordant level calls `classification/models/clinical_context_utils.py:update_clinical_context`, which
places the record in the `classification/models/clinical_context_models.py:ClinicalContext` for its allele, allele-origin
bucket and name (`default` unless a user moved it) and calls
`classification/models/clinical_context_models.py:ClinicalContext.recalc_and_save`. The verdict comes from
`classification/models/clinical_context_models.py:DiscordanceStatus.calculate`: each shared, non-withdrawn modification's
clinical significance is mapped to a bucket via `classification/models/evidence_key.py:EvidenceKeyMap.clinical_significance_to_bucket`,
and two buckets among the counted records is `DiscordanceLevel.DISCORDANT`; VUS-A/B/C differences and B-vs-LB are the
"concordant with differences" levels, and a somatic bucket is always `MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED`. If the
context is discordant but open pending-change flags would resolve it, `pending_concordance` is set - a live calculation,
never stored, which is why `DiscordanceReportTriage.is_outstanding` recomputes it.

`recalc_and_save` stores `last_evaluation` and, when no import is ongoing, sends `clinical_context_signal`; with
`settings.DISCORDANCE_ENABLED` (off by default, on for Shariant) the receiver calls
`classification/models/discordance_models.py:DiscordanceReport.update_latest`. A context that turns discordant opens a
`DiscordanceReport` with a `DiscordanceReportClassification` per shared record (original modification, final filled on
close); each later recalc `update`s it, and the moment the context is concordant again it `close`s as `CONCORDANT`.
A report closed as `CONTINUED_DISCORDANCE` reopens only if a new lab joins or the context becomes concordant
(`classification/models/discordance_models.py:DiscordanceReport.should_reopen_continued_discordance`). Flags follow the
report: `classification/models/discordance_models.py:DiscordanceReport.apply_flags_to_context` is the one place the
discordant flags on the context and on each classification are opened and closed. `discordance_change_signal` fans out to
per-lab `DiscordanceNotification` rows (batched into one email per lab per import, #3384) and to
`classification/models/discordance_models.py:ensure_discordance_report_triages_for`, which keeps a
`DiscordanceReportTriage` per actively involved lab (#3486; triage statuses "will amend", "for discussion", "confident").
The discussion itself is a `review` (see `claude/research/review.md`); its outcome view
`classification/views/discordance_report_views.py:action_discordance_report_review` raises the pending-changes flag with
`{"to_clin_sig": ...}` on the classifications a lab agreed to change, and
`classification/signals/classification_hooks_significant_change.py:clinical_significance_change_check` closes it when the
republished value arrives.

### Condition matching

Labs send condition as free text; ClinVar and cross-lab comparison need ontology terms. On publish,
`classification/models/condition_text_matching.py:ConditionTextMatch.sync_condition_text_classification` normalises the
text into a per-lab `classification/models/condition_text_matching.py:ConditionText` and ensures the hierarchy root →
gene symbol → mode of inheritance → classification exists as `ConditionTextMatch` rows. Terms set at any level apply to
every classification below it that has no override; `classification/models/condition_text_matching.py:apply_condition_resolution_to_classifications`
walks the tree on save and `classification/models/condition_text_matching.py:apply_condition_resolution` writes the
result into `Classification.condition_resolution`, raises the condition-resolution flag (history, #881) and sends
`condition_set_signal`. `classification/models/condition_text_matching.py:ConditionTextMatch.attempt_automatch` assigns
only what `classification/models/condition_text_matching.py:ConditionMatchingSuggestion.is_auto_assignable` is certain
of - an embedded id such as `MONDO:0005021` at the root, or a single leaf term with a known gene relationship at gene
level; everything else waits for a curator on the condition matching page. Multiple terms carry a
`classification/models/condition_text_matching.py:MultiCondition` of uncertain or co-occurring.

### ClinVar export

`classification/models/clinvar_export_prepare.py:ClinvarExportPrepare.update_export_records` (the `clinvar_export`
command, or the button on the export page) walks each ClinVarKey's shared, condition-resolved modifications per allele
and lets `classification/models/clinvar_export_prepare.py:ClinVarConsolidatingMerger` decide one
`classification/models/clinvar_export_models.py:ClinVarExport` per key + allele + condition, choosing the most recent
record and merging conditions that are the same or more specific. `classification/models/clinvar_export_convertor.py:ClinVarExportConverter.convert`
renders the submission as `ValidatedJson`: open flags (withdrawn, discordant, internal review, outstanding edits, pending
changes, "don't share with ClinVar"), unconfirmed variant matching and a somatic bucket all become embedded errors, and
an export with errors is `IN_ERROR` and skipped rather than raising. `classification/models/clinvar_export_models.py:ClinVarExportBatch.create_batches`
groups valid `NEW_SUBMISSION`/`CHANGES_PENDING` exports by key, allele-origin bucket and assertion criteria (10,000 per
batch) into `ClinVarExportSubmission` snapshots; `classification/models/clinvar_export_sync.py:ClinVarExportSync.next_request`
then drives one batch through submit → poll → fetch response file against the test or production API
(`settings.CLINVAR_EXPORT`), storing every exchange as a `ClinVarExportRequest` and writing returned SCVs back so the
next submission is an update rather than a novel record. `classification/models/clinvar_export_exclude_utils.py:published`
applies a key's `ClinVarKeyExcludePattern`s at publish time by opening the not-public flag.

## Why it is shaped this way

**Evidence is JSON with an EvidenceKey schema.** Labs disagree on which fields exist, what they are called and which
options they take, and Shariant onboarded them one at a time. A column per field would have meant a migration per lab.
`classification/models/evidence_key.py:EvidenceKey` rows carry type, options, category, order, `max_share_level`,
`mandatory`, `immutable` and `variantgrid_column` (which annotation column auto-populates it), so the form, the validator,
the CSV/REDCap exporters and the ClinVar mapping are all generated from the same rows, and adding a key is a data
migration (`ls classification/migrations | grep ekey`). `classification/models/evidence_key.py:EvidenceKeyOverrides`
layers lab and organisation config (hide a key, change options, enable a namespace such as `somatic:` or `acmg:`) on top,
merged by `classification/models/evidence_key.py:EvidenceKeyOverrides.merge` and cached per lab through
`classification/models/evidence_key.py:EvidenceKeyMap.with_overrides`. Each value is a `classification/models/evidence_mixin.py:VCBlobDict`
of value / note / explain / db_refs / validation / immutable, so a note and its provenance travel with the value.

**Every edit is a modification.** A classification is a clinical opinion that gets reported on and re-examined years
later; auditors ask what a lab said on a date, and discordance reports need to compare the version each lab had published
at detection against what they publish now. Storing deltas (`ClassificationModification.delta`) with a snapshot at publish
(`published_evidence`) gives both cheaply, and the read permission living on the modification rather than the
classification is what makes "you can see what was published, not what is being edited" a property of the data instead
of every view. The cost is that `Classification.evidence` is a denormalised copy that must be kept in step, which is why
the app insists on `patch_value`.

**Buckets come from key metadata, not the enum.** `classification/enums/classification_enums.py:ClinicalSignificance`
only knows B/LB/VUS/LP/P, but deployments add VUS-A/B/C, oncogenic tiers, risk alleles and "other" as options on the
`clinical_significance` key, and which of those count as the same bucket is a policy each deployment sets in the option's
`bucket` attribute. Reading the bucket from `EvidenceKeyMap` means a new option is a data change, not code, and the
`bucket` on `allele_origin` options is what `classification/enums/classification_enums.py:AlleleOriginBucket.bucket_for_allele_origin`
uses to split germline from somatic contexts.

**Allele resolution is shared and asynchronous.** Thousands of records cite the same HGVS; resolving once per distinct
input (`ImportedAlleleInfo`) rather than per classification made re-matching after a transcript or cdot upgrade tractable
(`classification_re_matching`, `fix_variant_matching` in `claude/maps/commands.md`). Matching goes through
the VCF pipeline because that is the only path allowed to create Variants, and liftover is a separate task chain because
it can take hours on a big import - so a classification has a nullable `variant` and `allele` for as long as that takes.

**Imports defer the expensive work.** One file can touch thousands of records, and recalculating discordance, groupings
and common-variant filters per row was both slow and noisy (an allele can flip discordant and back within one import).
`classification/models/classification_import_run.py:ClassificationImportRun.ongoing_imports` gates all three; the
`classification_imports_complete_signal` sweep does the work once, and `ClinicalContext.pending_cause` remembers what to
say in the notification. #1725 dates the guard to 2021-12 ("delay discordances until import completes").

**Discordance is per ClinicalContext, not per allele.** Two labs can legitimately classify the same allele differently
for different conditions (or germline vs somatic), so the unit of comparison is allele + allele-origin bucket + a name a
user can split records into, with `default` as the name nearly every record has.

## History

The models date from the 2020 rewrite (`7a00283d9`, "renamed lots of things, re-did flags"), with condition matching
reworked into a hierarchy that December. ClinVar export arrived in 2021-08 against ClinVar's test API; the import-run
guard in 2021-12; `UploadedClassificationsUnmapped` and the OmniImporter task in 2022-03; `ImportedAlleleInfo` in 2022-11,
replacing per-classification matching flags with `ImportedAlleleInfoValidation` rows in 2023-01. 2023 added the review
app and discordance triage (2023-08), medically-significant prioritisation and bulk discordance emails (#3486, #3384),
and condition resolution history (#881). Somatic support split the significance axis in 2024-02 ("classification" vs
"somatic clinical significance"), and `allele_origin` became mandatory in 2024-11 (variantgrid_private#2926).
`ClassificationGrouping` landed in 2025-02 for the listing grids. 2026 renamed the `institution` share level to
`organisation` (#1472, with `ShareLevel._missing_` accepting the old value), allowed g.HGVS-only imports (#1063), recorded
the HGVS converter version (#1321), made gene fusions first-class variants (#1506), added the external-lab filter (#394),
the CSV inserter (#1481), the `ReclassificationEvent` timeline behind the reclassification analytics page
(`classification/models/classification_reclassification_models.py:ReclassificationEventBuilder`, #1523), and the
denormalised `summary` JSON with its `summary__p_sort_idx` index that the grids sort on.

## Traps

- `Classification.variant` and `.allele` are nullable and are being retired in favour of `allele_object` (which reads
  `allele_info.allele`). Queries such as `ClassificationModification.latest_for_user(allele=...)` still go through
  `classification__variant__in`, so a record whose allele info resolved but whose `variant` column was never copied is
  invisible to them; `fix_allele_info` and `classification_ensure_alleles_and_liftover` repair that.
- `latest_for_user` without `published=True` and without a single `classification` filters on `is_last_edited`, which
  drops records that have an unpublished edit the caller could not see anyway; the FIXME in the method is real.
- `ClinicalContext.classifications_qs` excludes research labs unless `settings.DISCORDANCE_RESEARCH_ENABLED`, so a
  research lab's record can be in a context yet not counted - `lab_count_all` vs `lab_count` on `DiscordanceStatus`.
- `DiscordanceStatus._calculate_rows` silently ignores a clinical significance with no bucket (`has_ignored_clin_sigs`);
  a new option added without a `bucket` attribute makes records vanish from discordance rather than fail.
- `EvidenceKeyMap.instance` and `clinical_significance_to_bucket` are 60-second `timed_cache`s: a data migration that
  edits a key is not seen by a running worker for a minute, and tests that create keys in `setUp` must not rely on the
  cache having been warmed by an earlier test.
- Unknown keys are accepted by default (`settings.CLASSIFICATION_ALLOW_UNKNOWN_KEYS`, kept True so instances can sync from
  each other); `unknown_evidence_key_cleaner` is how they are found later. Shariant turns it off.
- `ImportedAlleleInfo.get_or_create` strips spaces from the HGVS before hashing but nothing else - `NM_000059.3:c.1A>G `
  and `nm_000059.3:c.1A>G` are two rows resolving to one allele.
- `ClassificationImportRun.record_classification_import` reuses an ONGOING run of the same identifier only if it was
  modified within `MAX_IMPORT_AGE` (5 minutes); a client that pauses longer starts a second run, and a run never marked
  `complete` is only closed by the next call's `cleanup`. Until then everything gated on `ongoing_imports` is deferred.
- `ClassificationModification.publish` compares the requested level with the record's current one and returns False when
  nothing changed - so a caller that patched and then asked to publish at the same level with no delta gets no new
  version and no signal; `BulkClassificationInserter.insert` relies on this to avoid empty republishes.
- The web form and API both go through `ClassificationView`, which is under the login-exempt `/classification/api/`
  prefix; keep any new view there a DRF `APIView` (`claude/guides/operations.md#authentication-surface`).
- ClinVar batches are germline only: `ClinVarExportSync.next_request` raises `NOT_SUPPORTED_YET` for any other
  `allele_origin_bucket`, and the converter marks somatic exports as errors first, so they never reach a batch.
- `ClassificationRef.init_from_str` treats a numeric id as a `Classification.pk` and anything else as `org/lab/record`;
  a lab whose `lab_record_id`s are integers must always be addressed with the lab prefix.
- Tests: `classification/tests/models/test_utils.py:ClassificationTestUtils` builds the lab/user pairs;
  `classification/tests/utils/test_urls.py` is the URLTestCase with `LIFTOVER_CLASSIFICATIONS=False`;
  `classification/tests/views/test_query_scaling.py` guards the grid query count. There are no tests for the OmniImporter
  path or ClinVar sync beyond `classification/tests/models/test_clinvar_export_sync.py` and
  `classification/tests/utils/test_clinvar_prepare.py`.
