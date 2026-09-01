# Trio X-linked analysis: use genetic sex, allow un-hiding template nodes (variantgrid_sapath#439)

Written by Claude Fable 5 (claude-fable-5), 2026-09-01

## Problem

Millennium/Helix enters prenatal cases under the mother's record, so a male fetus imports with
`patient.sex == FEMALE` (the sync code already notes this:
`variantgrid_sapath/sapath/models/sapath_helix_vg_sync.py:45`). The chain that follows:

1. `_xlinked_recessive_errors()` (`analysis/models/nodes/family_inheritance.py:44`) trusts
   `proband_sample.patient.sex` only, and errors on FEMALE: *"X-linked recessive inheritance
   doesn't currently work with female proband"*. Used by both `TrioNode.get_trio_inheritance_errors`
   and `QuadNode` (`quad_node.py:220`).
2. When the auto analysis template runs (`analysis/analysis_templates.py:52-71`), a node erroring
   with `hide_node_and_descendants_upon_template_configuration_error=True` gets `visible=False`,
   along with all its descendants.
3. Analysis views filter on `visible=True` (`analysis/views/views.py:203`,
   `analysis/views/views_json.py:68`), so the hidden X-linked branch is unrecoverable and leaves
   no visible trace — the trio silently lacks X-linked analysis.

## Fix (two parts, independent)

Part A makes the sex check consult genetic sex inferred from the sample's chrX genotypes, so new
trios get the X-linked branch automatically. Part B adds a way to reveal auto-hidden nodes, as a
manual escape hatch for existing analyses and any future auto-hide misfire.

## Part A — genetic sex from chrX genotype stats

Per-sample chrX hom/het counts already exist (`CohortGenotypeStats.x_hom_count` / `x_het_count`,
`snpdb/models/models_cohort_stats.py`), computed during VCF import by `calculate_vcf_stats`
(`annotation/tasks/calculate_sample_stats.py:43`) *before* samples are marked import SUCCESS — so
stats are normally present by the time any trio/template exists. `chrx_sex_guess`
(`models_cohort_stats.py:108`) already infers sex from the hom/het ratio.

### A1. `chrx_sex_guess` returns the enum, with a minimum-count guard

- Change `CohortGenotypeStats.chrx_sex_guess` to return `Sex` (currently returns the display
  label string). Update the one usage, `snpdb/templates/snpdb/data/view_sample.html:219`, to
  render `.label` via a small view/template adjustment.
- Add a minimum-evidence guard: return `Sex.UNKNOWN` unless `x_hom_count + x_het_count >=
  MIN_CHRX_VARIANTS_FOR_SEX_GUESS` (module constant, suggest 30). Small targeted panels have too
  few chrX calls for the ratio to mean anything.

### A2. `Sample.genetic_sex_guess` property

- Move the stats-row lookup from `_get_sample_genotype_stats` (`snpdb/views/views.py:656`) onto
  `Sample` in `snpdb/models/models_vcf.py` (e.g. `get_genotype_stats()` returning
  `Optional[CohortGenotypeStats]`, same resolution: own VCF cohort's CGC, per-sample row,
  `filter_key IS NULL`, `passing_filter=False`, swallowing the same missing-data exceptions).
- Add `genetic_sex_guess` property: `stats.chrx_sex_guess` when stats exist, else `Sex.UNKNOWN`.
- Update `snpdb/views/views.py` to call the model method.

### A3. Effective sex in the X-linked error check

In `_xlinked_recessive_errors` (`family_inheritance.py:44`):

- Effective sex = `proband_sample.genetic_sex_guess`, falling back to `patient.sex` when the
  guess is `UNKNOWN` (no stats row, legacy data, or below the count guard).
- Error on effective sex FEMALE, as now. The message should say which source decided, e.g.
  `"X-linked recessive inheritance doesn't currently work with female proband (genetic sex from
  chrX genotypes)"` vs `"... (recorded patient sex)"` — when a record-keeping error hides the
  node, the argument error stored on the `AnalysisTemplateRunArgument` then explains itself.

Behaviour changes to be aware of (both intended):

- Genetic MALE + recorded FEMALE (this issue): node no longer errors, X-linked branch survives
  templating. QuadNode gets the same fix for free.
- Genetic FEMALE + recorded MALE: now errors where it previously didn't — correct, and the count
  guard keeps low-evidence guesses from triggering it.

### A4 (optional). Pedigree icon rendering

`TrioNode.get_rendering_args` (`trio_node.py:370`) and `QuadNode` (`quad_node.py:274`) draw the
proband square/circle from `patient.sex`. Use the same effective-sex resolution so a genetically
male fetus renders male. Cheap once A2 exists; skip if scope-trimming.

## Part B — reveal auto-hidden template nodes

Hidden nodes still exist in the DB with `visible=False`; nothing in the UI shows they exist.

- New POST view + URL, e.g. `analysis/<analysis_id>/reveal_hidden_nodes/`, permission-checked
  with `analysis.check_can_write(user)`. Sets `visible=True` on the analysis's hidden nodes,
  bumps their versions / marks dirty as `populate_analysis_from_template_run` does in reverse,
  then `reload_analysis_nodes(analysis.pk)` so counts/status recompute. The revealed nodes render
  with their real errors showing, which is the point — the user sees *why* the branch was hidden.
- UI: in the analysis settings details tab (`analysis/templates/analysis/
  analysis_settings_details_tab.html`), when `analysis.analysisnode_set.filter(visible=False)`
  is non-empty, show the hidden-node count with each node's stored template argument error
  (`AnalysisTemplateRunArgument.error`) and a "Show hidden nodes" button posting to the new view.
- Nodes keep their saved x/y appearance, so revealing needs no layout work.

## Tests

- `chrx_sex_guess`: female/male/unknown ratio thresholds and the minimum-count guard.
- `_xlinked_recessive_errors`: genetic MALE overrides recorded FEMALE (no error); genetic
  UNKNOWN falls back to patient sex (error preserved); genetic FEMALE errors.
- Template run integration: trio template containing an X-linked node with
  `hide_node_and_descendants_upon_template_configuration_error=True`; recorded-female/genetic-male
  proband → node stays visible; recorded-female/no-stats proband → node and descendants hidden
  (existing behaviour).
- Reveal view: hidden nodes become visible, and a read-only user gets permission denied.

Audit the tests when done per CLAUDE.md — keep only the ones covering our branching logic.

## Out of scope

- Fixing sex at the source in the Helix sync (`variantgrid_sapath`) — the mother's record is what
  Millennium sends; genetic sex is the reliable signal.
- MOINode / other inheritance modes — only the X-linked recessive proband check reads sex.
