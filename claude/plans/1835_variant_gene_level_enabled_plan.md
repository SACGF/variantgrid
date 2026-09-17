# VARIANT_GENE_LEVEL_ENABLED: one setting that turns gene-level variants off

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-17
Status: in progress
Issue: #1835 (gene-level classifications), #1876 (gene-level search)

A deployment that classifies germline small variants only (Shariant) has no use for fusions, whole-gene
copy number events or splice events, and every surface that accepts one is a way for a stray value to
mint a Variant on the fake contig. `VARIANT_GENE_LEVEL_ENABLED` is the switch, shaped like
`VARIANT_SYMBOLIC_ALT_ENABLED`: a default-on boolean, refused at the coordinate so every write path is
covered, and read by the surfaces that would otherwise offer the feature.

## Data

No model changes. The setting, in `variantgrid/settings/components/default_settings.py` beside the symbolic-alt block:

```python
# Gene-level variants - fusions, whole-gene copy number and splice events stored as Variants on the
# fake contig (@see snpdb.gene_level_variants). Off, a value naming one is refused at the coordinate,
# the search shortcuts and TSO 500 upload types are withheld and the gene-level annotation pipeline
# is never scheduled.
VARIANT_GENE_LEVEL_ENABLED = True
```

`ANNOTATION_GENE_LEVEL_ENABLED` (`variantgrid/settings/components/annotation_settings.py`) is folded into it: its
only reason to exist was "a deployment with no gene-level variants gets an always-empty run per range lock", which
is exactly the case the new setting names. Deployment settings that set it:

| File | Change |
|---|---|
| `variantgrid/settings/env/shariantcommon.py` | `ANNOTATION_GENE_LEVEL_ENABLED = False` becomes `VARIANT_GENE_LEVEL_ENABLED = False` |
| `annotation/fake_annotation.py` | the override dict key renamed the same way |

`../variantgrid_sapath` reads neither name, so nothing changes there.

## Where it gates

| Surface | Code | When off |
|---|---|---|
| Any write of a gene-level coordinate | `snpdb/models/models_variant.py:VariantCoordinate` | a `model_validator` raises `ValueError('Gene-level variants disabled via settings.')` for a coordinate on `GENE_LEVEL_CONTIG_NAME`, the twin of `validate_svlen` |
| Classification import (API, file, web form, sync) | `genes/gene_level_strings.py:looks_gene_level` | returns False, so no value is ever recognised as gene-level: it takes the HGVS path and fails as an HGVS |
| Imported Allele Info page | `classification/templates/classification/imported_allele_info.html` | the "Gene Level Unresolved" filter is rendered only when enabled |
| Search | `snpdb/signals/variant_search.py` | the three receivers pass `enabled=settings.VARIANT_GENE_LEVEL_ENABLED`, so they are never connected and their examples never listed |
| Upload file types | `upload/import_task_factories/import_task_factory.py:ImportTaskFactory` | an `enabled` property, default True; the four gene-level factories return the setting and `get_import_task_factories` skips disabled ones |
| API capabilities | `seqauto/views_rest.py:CapabilitiesView` | `upload_file_types` is derived from the factories, so the TSO 500 and gene-level CNV types drop out with no change here |
| Annotation | `annotation/pipelines/__init__.py` | the `GeneLevelRunner` `PipelineDef` names `VARIANT_GENE_LEVEL_ENABLED` as its `enabled_setting` |
| All Variants page | `snpdb/variant_filters.py:get_all_variant_types` | `VariantType.FUSION` is appended only when enabled, the way `get_symbolic_variant_types` already works |
| Server status | `variantopedia/templates/variantopedia/server_status_settings_detail.html` | a "Gene-level Variants Enabled" row under the Variant Matching card, beside the symbolic one |

### Classification import in detail

`genes/gene_level_strings.py:looks_gene_level` is the one decision every classification write path makes (API
v1-v3, file upload, the create form, Shariant sync all reach it through `ImportedAlleleInfo.get_or_create`,
and `ImportedAlleleInfo.is_gene_level` reads it for a record with no coordinate). With the setting off it
returns False before matching any pattern, and the feature is absent:

- `resolve_gene_level_string` answers `not_applicable`, so `resolve_gene_level` returns False and the value
  goes to the HGVS converter, where 'BCR::ABL1' fails to parse the way it did before gene-level existed. The
  record carries the ordinary `cant_resolve_to_variant_coordinate` tag and the converter's message.
- `is_gene_level` is False, so `imported_as_c_hgvs` and `_calculate_validation` treat it as a c.HGVS
  submission; `gene_level_unresolved` is never written.
- `classification/classification_import.py:_classification_upload_pipeline` needs no change: with no gene-level
  coordinates the gene-level insert pipeline never runs.
- The Imported Allele Info page renders its "Gene Level Unresolved" checkbox only when the setting is on
  (`classification/views/imported_allele_info_view.py` passes it into the context); the query parameter it
  drives stays as it is.

The `VariantCoordinate` validator sits underneath this as the guarantee: a coordinate built any other way
(`from_string` on a `GENE_LEVEL:` string, a VCF the CNV factory would have claimed, a management command) is
refused at construction.

### Search in detail

`search_receiver(enabled=...)` is evaluated at import, the same as `SEARCH_COSMIC_ENABLED`, so the receivers
for `search_variant_gene_fusion`, `search_variant_gene_copy_number` and `search_variant_splice_event` are
simply never connected. A disabled deployment's search help lists none of their examples. "EGFR amplification"
then falls through to the ordinary receivers and returns whatever they make of it.

### Upload factories in detail

`get_import_task_factories` grows one line: `if itf.enabled`. `ImportTaskFactory.enabled` returns True; the
override on `DragenTSO500AllFusionsImportTaskFactory`, `DragenTSO500CombinedVariantOutputImportTaskFactory`,
`GeneLevelCNVImportTaskFactory` and `GeneLevelInsertVariantsOnlyImportFactory` in
`upload/import_task_factories/import_task_factories.py` returns `settings.VARIANT_GENE_LEVEL_ENABLED`. A CNV
VCF whose header declares a `SEGID` field is then claimed by the ordinary VCF factory and imported on its
written coordinates, which is what happened before gene-level CNVs existed.

`seqauto/views_rest.py:API_FEATURES` stays as it is. A client learns the server's gene-level status from
`upload_file_types`, which already reflects the factories.

## Existing data

A deployment that flips the setting off with gene-level Variants already in its database keeps them: rows are
read, displayed, synced (the `remote_gene_level` filter in `sync/shariant/variant_grid_upload.py` is per
destination and stays) and exported as today. The setting stops new ones arriving. Shariant has none, so
this is the germline-only case the setting is for.

## Tests

`override_settings(VARIANT_GENE_LEVEL_ENABLED=False)` covers the runtime reads:

- `snpdb/tests`: `VariantCoordinate.from_string` on a `GENE_LEVEL:` string raises; the same string builds
  with the setting on.
- `genes/tests/test_gene_fusions.py` (or a sibling): `looks_gene_level('BCR::ABL1')` is False under the setting.
- `classification/tests/models/test_imported_allele_info.py`: `get_or_create(imported_c_hgvs='BCR::ABL1', ...)`
  fails as an HGVS - `variant_coordinate` empty, `status` FAILED, `cant_resolve_to_variant_coordinate` tagged,
  `gene_level_unresolved` absent - and no gene-level `UploadedFile` is created.
- `upload/tests`: `get_import_task_factories` returns no gene-level factory, and the capabilities endpoint's
  `upload_file_types` omits them (`seqauto/tests/test_capabilities.py` is where that endpoint is tested).
- `snpdb/tests`: `get_all_variant_types` ends without `FUSION`.
- `annotation/tests/test_gene_level_annotation.py`: `enabled_pipeline_types` drops `GENE_LEVEL` under the
  renamed key.

The search receivers' `enabled` flag is framework behaviour evaluated at import; `variantopedia/tests/test_gene_level_search.py`
keeps testing the enabled path and the flag itself gets no test.

## Docs

- `snpdb/gene_level_variants.py` docstring gains a "Turning it off" paragraph naming the setting and the four
  places it is read - this file is where every `@see` points.
- `claude/guides/operations.md` deployments section: Shariant runs with `VARIANT_GENE_LEVEL_ENABLED = False`.
- `classification/CLAUDE.md`: one line under the gene-level gotchas saying `looks_gene_level` is False on a
  disabled server, so a gene-level value fails as an HGVS there.
- `upload/CLAUDE.md`: one line that `ImportTaskFactory.enabled` is how a settings-gated file type withdraws
  from upload and from the capabilities endpoint.

## Order of work

1. Setting, `ANNOTATION_GENE_LEVEL_ENABLED` rename (settings, `PipelineDef`, `fake_annotation.py`).
2. `VariantCoordinate` validator and its test.
3. `looks_gene_level` gate, the Imported Allele Info filter, and their tests.
4. Search receivers.
5. `ImportTaskFactory.enabled`, the four overrides, `get_import_task_factories`, tests.
6. `get_all_variant_types`, server status row.
7. `snpdb/gene_level_variants.py` docstring and the CLAUDE.md lines; `scripts/vg docs check`.
