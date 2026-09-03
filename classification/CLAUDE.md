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
Tests:
- tests/models/test_utils.py:ClassificationTestUtils.setUp builds the Org / Lab / user pair (lab_and_user, external_lab_and_user);
  create records with Classification.create(user=..., lab=..., data={SpecialEKeys.C_HGVS: {'value': ...}}, source=SubmissionSource.API).
- tests/utils/test_urls.py is this app's URLTestCase (create_fake_variants, then admin / owner / non-owner passes with 403 checks);
  tests/views/test_query_scaling.py asserts the classification datatable query count stays flat as rows grow.
- URLTestCase already sets CELERY_TASK_ALWAYS_EAGER and LIFTOVER_CLASSIFICATIONS=False, so import tests run inline without liftover.
- tests/utils/data_utils.py:ConditionMock is stale by its own docstring (the ontology mock no longer builds); mock ontology lookups instead.
Deep reference: __classification_readme.md · __discordance_readme.md · claude/research/classifications.md · claude/maps/models.md#classification
