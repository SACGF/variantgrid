# classification — agent notes
Owns: Classification / ClassificationModification, EvidenceKey + EvidenceKeyMap, ShareLevel, ImportedAlleleInfo, ClinicalContext +
DiscordanceReport, ConditionText matching, ClassificationGrouping, ClinVarExport, the classification API and import pipeline.
Start with:
- models/classification.py — Classification, ClassificationModification, the app's signals (top of file), patch_value / publish.
- models/evidence_key.py — EvidenceKey (schema of the evidence JSON), EvidenceKeyMap (cached lookup + lab overrides), VCDataCell.
- enums/classification_enums.py — ShareLevel, SpecialEKeys, ClinicalSignificance, SubmissionSource, CriteriaEvaluation.
- models/classification_variant_info_models.py — ImportedAlleleInfo, ResolvedVariantInfo (HGVS to Allele, per-build c.HGVS cache).
- models/clinical_context_models.py + models/discordance_models.py — ClinicalContext, DiscordanceStatus, DiscordanceReport.
- models/classification_inserter.py — BulkClassificationInserter, the one write path shared by the API, file imports and sync.
Patterns here:
- Evidence is JSON keyed by EvidenceKey.key; each value is a blob of value/note/explain/db_refs/validation
  (models/evidence_mixin.py:VCBlobDict). Read it with models/evidence_mixin.py:EvidenceMixin.get and name keys via
  enums/classification_enums.py:SpecialEKeys rather than string literals.
- Change evidence only through models/classification.py:Classification.patch_value (or Classification.create_with_response /
  models/classification_inserter.py:BulkClassificationInserter.insert): it validates, honours per-field immutability by
  SubmissionSource, refreshes the cached clinical_significance / allele_origin_bucket and writes the ClassificationModification.
  Never assign into Classification.evidence directly.
- Publish via Classification.publish_latest, which calls models/classification.py:ClassificationModification.publish: flips
  is_last_published, grants the ShareLevel group the read perm and sends classification_post_publish_signal. Only published
  modifications are visible outside the owning lab; view permission lives on ClassificationModification, Classification itself
  carries only the write perm (models/classification.py:Classification.filter_for_user).
- Query what a user may see with models/classification.py:ClassificationModification.latest_for_user (published=True, allele=...,
  shared_only=...). It applies Guardian perms plus withdrawn and allele-origin filters; pass allele= rather than variant=.
- Parse "org/lab/record_id.version" ids with models/classification_ref.py:ClassificationRef.init_from_str; lab_record_id is unique per lab.
- Look keys up through models/evidence_key.py:EvidenceKeyMap.instance (60s timed_cache) and apply lab config with
  EvidenceKeyMap.with_overrides(classification.evidence_key_overrides). EvidenceKey.max_share_level hides a field from users
  below that level (models/classification.py:Classification.get_visible_evidence).
- Link a classification to a variant only via ImportedAlleleInfo: models/classification.py:Classification.ensure_allele_info ->
  models/classification_variant_info_models.py:ImportedAlleleInfo.get_or_create (unique on md5 of the imported HGVS + transcript
  + build patch). allele_info_changed_signal then fans the resolution out to classifications, groupings and clinical contexts.
- Discordance is per ClinicalContext (allele + allele_origin_bucket + name), recalculated on publish / withdraw / delete by
  models/clinical_context_models.py:ClinicalContext.recalc_and_save. Buckets come from the "bucket" attribute on the
  clinical_significance EvidenceKey options (models/evidence_key.py:EvidenceKeyMap.clinical_significance_to_bucket), not from
  the ClinicalSignificance enum, and only records at a ShareLevel.is_discordant_level (logged_in_users, public) count.
- Hook the lifecycle with the signals at the top of models/classification.py (classification_validation_signal,
  classification_post_publish_signal, classification_withdraw_signal, classification_variant_set_signal,
  classification_revalidate_signal, variants_classification_changed_signal). Receivers live in signals/ and are imported by
  apps.py:ClassificationConfig.ready; annotation/apps.py and snpdb/signals/common_variants_classification_changed.py subscribe
  from outside the app.
- A lab's `classification_config` can give an evidence key a `default_value`, which a record created from the web
  form starts that key with (`models/classification.py:Classification.create_with_response` with
  `populate_with_defaults=True` - the API, file imports and sync never see one). `$user` and `$today` resolve at create
  time through `models/evidence_key.py:resolve_default_value`, so a lab that makes `curated_by` / `curation_date`
  mandatory can say what they start as.
- The case report (#444) is one Django template per lab rendered server side: report/case_report_context.py builds
  ReportVariant / ReportContext (ordering, amp_tier, kinds, measures), report/renderers.py turns one HTML into the PDF
  (xhtml2pdf) and DOCX (html2docx), and the JSON comes off the same context. ClassificationReport.context() builds its
  rows from the same ReportVariant, so `record` / `classifications` / `gene_groups` mean one thing.
- One run is a CaseReport, built / rebuilt / finalised through report/case_report_builder.py. `Rebuild documents`
  re-renders over `context_snapshot` so a template fix never moves the numbers a report has quoted; `New version`
  builds a fresh context from each pinned classification's current published version and marks the old report
  SUPERSEDED. Finalising stamps report_date, variant_reported and report_id onto each pinned classification, publishes
  it and re-points the CaseReportClassification at the new published version - all through patch_value, so
  re-finalising or re-entering the same LIS id writes nothing.
- Going DRAFT -> FINAL sends models/classification_report_models.py:case_report_finalised_signal (once - re-finalising
  does not), which is how a deployment specific app files the report somewhere else; what it made of it comes back to
  the Reports card through case_report_deliveries_signal as library/case_report_delivery.py:CaseReportDelivery rows
  (models/classification_report_models.py:get_case_report_deliveries collects them). SA Path answers both to send
  finalised TSO 500 reports to Mocha.
- A report's JSON is the app that owns that template's shape answering case_report_json_signal, else the canonical
  context dump (report/renderers.py:render_json). A JSON another system parses is an interface, so the app that has
  to keep it in step with that system writes it in Python (SA Path's TSO 500 shape) - there is no JSON template.
- A CaseReport's permissions are its own (models/classification_report_models.py:CaseReport.can_view / can_write):
  its lab's users may build, finalise and rebuild, and anyone who can see the case may read. Every document is served
  by views/views_case_report.py rather than a media URL - MEDIA_ROOT has no permissions of its own.
Gotchas:
- "Every edit makes a ClassificationModification" has one exception: models/classification.py:ClassificationModification.is_edit_appendable
  folds a patch into the previous unpublished modification when it is the same user and source within a minute.
- Share level only ratchets up, but the model does not enforce it: BulkClassificationInserter.insert re-publishes at the current
  level with a "shared_higher" warning when asked for a lower one, and the web form always republishes at the current level
  (views/views.py:create_classification_object). ShareLevel keys are user / lab / organisation / logged_in_users / public;
  there is no "institution".
- Classification.variant and .allele are nullable: a record whose ImportedAlleleInfo failed validation keeps variant=None and
  is excluded from exports (models/classification.py:Classification.include_based_on_allele_info). Filter on it rather than
  assuming every classification has a variant.
- settings DISCORDANCE_ENABLED defaults to False: ClinicalContext.recalc_and_save still stores the status, but DiscordanceReports
  and notifications are only created when it is on. models/discordance_models.py:DiscordanceReport.is_pending_concordance is
  a live calculation, never stored.
- Bulk-import guard: while a ClassificationImportRun is ONGOING, per-classification work (grouping, common-variant filters,
  clinical-context recalc) is skipped and caught up on classification_imports_complete_signal
  (models/classification_import_run.py:ClassificationImportRun.ongoing_imports). An abandoned ONGOING run silently stalls all of it.
- /classification/api/* is in PUBLIC_PATHS, so GlobalLoginRequiredMiddleware skips it: every view under that prefix must be a
  DRF APIView so the IsAuthenticated default applies (views/classification_view.py:ClassificationView is the model). A plain
  Django view there is anonymous.
- EvidenceKeys are rows seeded by data migrations (ls migrations | grep ekey): adding or renaming a key is a migration, and
  EvidenceKeyMap.instance caches for 60s. Tests that need a key create it with EvidenceKey.objects.create in setUp
  (tests/views/test_classification_view.py).
- What copying from a previous classification brings across is EvidenceKey.copy_scope (NONE / ALLELE / GENE) filtered by
  copy_allele_origin (models/classification.py:ClassificationConsensus.consensus_patch). Scope is how far a value travels, so
  a GENE key also copies at allele level; copy_allele_origin=GERMLINE keeps segregation and de novo data out of somatic records.
  Pass copy_scopes=COPY_SCOPES_GENE for a gene-only copy, and models/classification.py:ClassificationConsensus.apply_to to
  write one into a record (fills empty fields only, as SubmissionSource.CONSENSUS).
- The target allele origin bucket is decided before the candidates are listed, never after: nothing germline is ever offered
  as the source for a somatic record. The create page takes it from the user's allele_origin_focus (flippable), the Classify
  & Report dialog from Tag.allele_origin_bucket, the in-form box from the record's own bucket. External labs' records are
  shown for context with no copy control - their evidence was assembled under a config reviewed elsewhere.
- Gene level candidates come from models/classification.py:ClassificationConsensus.gene_consensus_groups - one row per
  distinct set of GENE-scope values, newest representative first, capped at ten. Most records in a gene carry identical gene
  content because they were copied from each other. The pick is always a human's: AMP tiering and therapy content are gene
  *and* tumour type, so the row's "and N other records" spread is the deciding information, not something to automate.
  The three surfaces that offer it (create page, Classify & Report dialog, the form's Gene Content card in
  views/views_gene_consensus.py) all call that one function.
- Classification.clinical_significance, allele_origin_bucket and summary are denormalised caches written by patch_value / publish
  (models/evidence_mixin_summary_cache.py:ClassificationSummaryCalculator); filter and sort on them, never recompute from
  evidence in a grid.
- Condition resolution flows from ConditionTextMatch to classifications via models/condition_text_matching.py:apply_condition_resolution
  (sends condition_set_signal); Classification.condition_resolution is a cache of that match, not the source.
- ClinVar export is one ClinVarExport per ClinVarKey + allele + condition (models/clinvar_export_models.py:ClinVarExport), built by
  models/clinvar_export_convertor.py:ClinVarExportConverter as ValidatedJson: embedded errors keep it out of a batch, no exception is raised.
- Withdrawing is a soft delete (models/classification.py:Classification.set_withdrawn) that sends classification_withdraw_signal;
  the record stays in the lab's lists and can be un-withdrawn. A hard delete goes through pre_delete / post_delete receivers
  that recalc contexts and groupings.
- Finalising leaves a classification with unsubmitted edits alone and names it
  (report/case_report_builder.py:has_unsubmitted_edits): publishing it would push someone's work in progress out under
  the report's name. The finalise response lists them and the tab alerts, rather than silently skipping.
- `variant_reported` says which *kind* of finding a variant was (primary / secondary / incidental), which is the
  curator's call - so finalising only ever settles reported against `not_included`, and a record that already names a
  kind keeps it.
- html2docx has no head: it prints the contents of `<style>` and `<script>` as the Word file's first paragraph, so
  report/renderers.py:render_docx strips them. A template's print CSS is design, not report content.
- A column aligned Results Summary survives html2docx as a `<table>` (a real `<w:tbl>`, equal column widths) and not as a
  `<pre>` (spaces and `<w:br/>`s, but no `<w:rFonts>`, so Word sets it proportional): the TSO 500 template keeps `<pre>`
  because its PDF is read against the legacy document, and report/renderers.py:render_docx opens a `<code>` run at each
  `<pre>` and rewrites html2docx's "Mono" to Courier New - closing that `<code>` instead loses every line break.
- `kind_groups` always carries all five kinds, empty variants and all (the way the tiers are always printed), so a
  template can say "Gene Fusions - None detected". A Results Summary that skips a kind reads as though it wasn't looked for.
  The legacy TSO 500 report is the one exception, omitting Splicing Variants entirely when empty - the template decides.
- An amplification carries two magnitudes and they are different quantities: `copy_number` is the caller's
  absolute count (VCF `CN`) and `fold_change` its ratio against the normal (VCF `SM` / `FC`). Both are
  autopopulated from the sample genotype, routed on `VCF.copy_number_field` by
  library/genomics/vcf_enums.py:VCFConstant.COPY_NUMBER_FIELD_IS_RATIO, so a record may hold either or both.
- What kind of event the case report prints a record as comes off its gene-level alt
  (`classification/report/case_report_context.py:_kind_and_alteration`) - a fusion, a copy number call or a splice
  junction. A SpliceGirl VCF imports as gene-level splice junctions (`upload/tasks/import_splicegirl_vcf_task.py`);
  the coordinate `<DEL>`s it made before #1903 still read as small variants - the kind is the alt's, not the caller's.
  The printed name of the event ("MET exon 14 skipping") is the `splice_label` evidence
  key, autopopulated from the junction's label through `genes/gene_splice.py:display_splice_label` - the panel's own
  wording where a `genes/models/models_splice_event.py:SpliceEvent` names it, and the breakpoints written out
  ("AR GRCh37 X:66905968-66914514") where nothing does, which is the prompt for the scientist to name it.
- A gene-level classification target is named rather than given as HGVS: `models/classification_variant_info_models.py:ImportedAlleleInfo.resolve_gene_level`
  runs the fusion, whole-gene copy number and splice string resolvers through `genes/gene_level_strings.py:resolve_gene_level_string`
  before any HGVS conversion. A gene-level value keeps its spaces - `tidy_hgvs_whitespace` strips them from an HGVS
  only, so `EGFR amplification` reads the same in c.HGVS as in the g.HGVS annotation writes; records imported before
  that hold `EGFRamplification`, which `gene_level_strings_respace` puts right. What takes that path is the *shape* of the
  value (`genes/gene_level_strings.py:looks_gene_level`), so one whose gene turned out to be a typo fails as a
  gene-level record - `gene_level_unresolved`, with the resolver's reason as its message - rather than as a broken
  HGVS. A splice label shape has to be in `genes/gene_splice.py:SPLICE_STRING_PATTERN` to be recognised at all;
  nothing has to be pre-registered beyond that. claude/research/classifications.md has the table of written forms.
- With `VARIANT_GENE_LEVEL_ENABLED` off (Shariant), `genes/gene_level_strings.py:looks_gene_level` is always False, so
  'BCR::ABL1' takes the HGVS path and fails as an HGVS - `cant_resolve_to_variant_coordinate`, never `gene_level_unresolved`.
- A gene-level variant sits on no transcript, so `gene_symbol` is autopopulated from the event's GeneLevelId (a fusion's
  anchor first), which already holds the approved symbol the caller's MYCL1 resolved to
  (`autopopulate_evidence_keys/evidence_from_variant.py:get_evidence_fields_from_gene_level_event`) - the transcript path
  fills that key for everything else.
- A gene-level record's validation must not be read as a c.HGVS submission: `imported_as_c_hgvs` returns False when
  `ImportedAlleleInfo.is_gene_level`, because the named value ('ARV7') lands in `imported_c_hgvs` but sits on no
  transcript. Otherwise `_calculate_validation` tags `transcript_type_not_supported` as "E" and `should_include` keeps
  every gene-level record out of exports forever, and the per-build check demands a c.HGVS `ResolvedVariantInfo`
  deliberately never writes for one.
- A `case_field` can carry `prefill_key`: the build form starts that field from the named evidence key on the case's
  first classification that has one (SA Path's clinical indication), and its `default` otherwise.
- A bool `case_field` can also carry `measure` (a `patients/models_enums.py:MEASURE_CONTEXT_KEYS` value) and `tick_when`:
  the form shows that `SpecimenMeasure` beside the checkbox and starts the tick from the rule - `{"called": true}`,
  `{"call_in": [...]}`, `{"value_below": n}`, or a list of those that ticks when any holds
  (`classification/models/classification_report_models.py:measure_tick`). Under each group the form lists the rule in
  words (`describe_tick_when`) and the measure's `threshold` / `threshold_source`, which the TSO 500 import fills from
  the `TSO500_*_CALL_BANDS` settings - the policy is on the form so a scientist who disagrees can ask for it to change.
  A case with no such measure falls back to `default`, and a draft's own answer wins over both. The keys are hand
  written in admin, so `validate_case_fields` fails the save on an unknown measure or rule.
- A bool `case_field` carries `qc` (a `seqauto/models/models_enums.py:LIBRARY_QC_CONTEXT_KEYS` value) instead of `measure`
  where the flag means "did the caller's library QC pass for that category" - the Amplifications / Variants / Fusions
  flags. Its `tick_when` rules are `{"passed": true|false}` and `{"completed": false}` against the category's
  `seqauto/models/models_seqauto.py:LibraryQC` row (`classification/models/classification_report_models.py:library_qc_tick`;
  `tick_for` picks which of the two a field is judged by), the newest run per category winning
  (`report/case_report_context.py:specimen_library_qc`, which reads the specimen's rows directly - a `LibraryQC` is
  keyed on the run and the pair, and claims the specimen rather than either extraction). The form lists the category's
  metrics against the guidelines DRAGEN's own file quotes, `settings.TSO500_LIBRARY_QC_GUIDELINES` overriding one
  where the lab's own number differs.
- report/__init__.py stays empty on purpose: models/classification_report_models.py imports report/template_validation.py
  for the save-time fixture render, so a package __init__ that reached into classification.models would be a cycle. For
  the same reason template_validation.FIXTURE_CONTEXT is hand written rather than built from ReportContext -
  tests/report/test_case_report_context.py:FixtureContextTest is what keeps the two in step.
Tests:
- tests/models/test_utils.py:ClassificationTestUtils.setUp builds the Org / Lab / user pair (lab_and_user, external_lab_and_user);
  create records with Classification.create(user=..., lab=..., data={SpecialEKeys.C_HGVS: {'value': ...}}, source=SubmissionSource.API).
- tests/utils/test_urls.py is this app's URLTestCase (create_fake_variants, then admin / owner / non-owner passes with 403 checks);
  tests/views/test_query_scaling.py asserts the classification datatable query count stays flat as rows grow.
- URLTestCase already sets CELERY_TASK_ALWAYS_EAGER and LIFTOVER_CLASSIFICATIONS=False, so import tests run inline without liftover.
- tests/utils/data_utils.py:ConditionMock is stale by its own docstring (the ontology mock no longer builds); mock ontology lookups instead.
Deep reference: __classification_readme.md · __discordance_readme.md · claude/research/classifications.md · claude/maps/models.md#classification
