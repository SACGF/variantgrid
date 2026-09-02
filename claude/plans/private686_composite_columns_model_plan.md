# Variant grid composite columns - a model, not a Python dict

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-02. Implemented and revised by Claude Opus 5 (claude-opus-5), 2026-09-02

Follow-on to the variantgrid_private#686 grid work (commits `40e457d0b`..`24abe8466`). Composite cells
become catalogue entries with members recorded in the database, the custom columns arrange page shows
composites instead of their members, a generic renderer covers the many "headline plus detail" groups,
and every collection on every deployment (user and system) is migrated onto the full set of composites
in one release.

---

## 1. Where #686 left things

A composite cell today is an ordinary catalogue column (the "anchor") given a richer client renderer.
`COMPOSITE_COLUMN_ROW_FIELDS` in `snpdb/grid_columns/custom_columns.py` maps each anchor to the partner
columns its renderer reads, and `get_variant_grid_columns` selects those partners out hidden whenever the
anchor is in the collection. Two more Python lists (`VARIANT_COLUMN_ROW_FIELDS`,
`CLASSIFICATIONS_COLUMN_ROW_FIELDS`) do the same for the representative variant and classifications
cells. Migrations 0225-0230 removed the partners from the **system** collections only
(`user__isnull=True`, "All columns" exempt).

Consequences on a deployment:

- Every user-owned collection (clones of the old defaults included) still holds the partners, so those
  users see the SpliceAI cell **and** eight separate SpliceAI columns, the gnomAD AF cell **and** a
  separate Pass/Fail column, and so on.
- The arrange page has no vocabulary for a composite. Dragging "gnomAD AF" silently brings a hidden
  partner; dragging a partner in as well shows it twice; picking only partners never surfaces the cell.
- Membership is a Python constant, so the page, the deployment column check and future composites all
  have to be hand-wired.
- 277 catalogue columns to drag around.

## 2. Decisions (agreed in discussion)

1. **A composite is its own catalogue entry.** No anchor. It is a display-only `VariantGridColumn`
   (the same kind as `classifications` today), with an ordered list of member columns.
2. **Members are hidden from the arrange page by default.** They stay in the catalogue, because filter
   nodes, evidence-key autopopulate, the VEP column map and the annotation version diff all key on them,
   but they are not offered as things to drag. A checkbox reveals them for the rare lab that wants one
   standalone. A member added that way renders as a normal column beside the composite.
3. **Every collection migrates**, user and system, "All columns" included. Exports lose nothing: a
   hidden member is still in the CSV/VCF whenever its composite is present.
4. **Rendering stays in JS; membership and sort menus come from the model.** A generic
   first-member-plus-hover renderer draws most composites; a handful keep bespoke renderers.
5. **All composites land in this release.** The collapse migration touches every collection once.
   Later additions are cheap (§9) but the bulk of the catalogue is grouped now.

**Design rule for what to combine.** Combine when either (a) the members are usually blank together,
so one empty cell lets the eye skip (SpliceAI, Mastermind, denovo-db, MaveDB), or (b) one member is the
headline everyone scans and the rest are context only read after it catches the eye (gnomAD, the score
families). Where two members are read independently side by side, leave them separate. Identifiers
and always-filled gene-level text (pathways, GO, expression, interactions, UniProt text) are a width
problem, not a skim problem, and stay as they are.

## 3. Data

`snpdb/models/models_columns.py` - one new table. `VariantGridColumn` gains no fields: a composite is
an existing-shape row (display-only via `queryset_field=False`, or keyed like `variant`), and
`CustomColumnsCollection` / `CustomColumn` are untouched.

```python
class CompositeColumnMember(models.Model):
    composite = models.ForeignKey(VariantGridColumn, on_delete=CASCADE, related_name="composite_members")
    column = models.OneToOneField(VariantGridColumn, on_delete=CASCADE, related_name="composite_membership")
    sort_order = models.IntegerField()
    in_sort_menu = models.BooleanField(default=True)

    class Meta:
        ordering = ("composite", "sort_order")
```

What the shape commits us to:

- **One composite per column** (`OneToOneField`). Relaxing it later means a plain FK plus
  `unique_together`, with the hidden-field logic deduplicating - cheap if ever needed. §6.3 holds under
  the strict rule today.
- **The headline is the first member** (`sort_order` 0): the default sort key, and what the generic
  renderer prints. Everything after it is hover detail. If a renderer ever needs to tell members
  apart beyond that (a tag drawn inline, a threshold), it's a column added to this table then.
- **Composites nest one level.** A member is a real column; a composite can't be a member. Nesting would
  need the grid code recursing - out of scope, and nothing in §6.3 wants it.
- **Column identity stays the pk** (`grid_column_name`), so composites are addressed by a stable string
  from migrations, overrides and the JS exactly as columns are now.

Code-side helpers on `VariantGridColumn` (`is_composite`, ordered members, headline, derived sort menu,
membership css classes) are implementation detail.

## 4. Grid

### 4.1 `get_variant_grid_columns` (`snpdb/grid_columns/custom_columns.py`)

Delete `COMPOSITE_COLUMN_ROW_FIELDS`, `VARIANT_COLUMN_ROW_FIELDS` and
`CLASSIFICATIONS_COLUMN_ROW_FIELDS`. `CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS` stays - those are
`get_variantgrid_extra_annotate` aliases and become members like any other catalogue column, keyed by
their alias.

For each `CustomColumn` in the collection whose column is a composite:

- Build the visible `RichColumn` with the composite's label/description/width. A display-only composite
  gets `model_field=False`, `queryset_field=False`, `orderable=True` and
  `sort_keys=[first_member.variant_column]` (`RichColumn` refuses an orderable column with neither `key` nor
  `sort_keys`, and the JS sort menu skips non-orderable columns). A keyed composite (`variant`) keeps
  its key.
- `sort_menu` is derived: `[{"label": m.label, "column": m.variant_column} for m in members if m.in_sort_menu]`,
  in member order. The hand-written `sort_menu` literals in `_get_standard_overrides` go.
- `client_renderer_kwargs` gets `{"members": [{"path", "label"}, ...]}` in member order so the generic
  renderer knows what to draw (§4.4). Bespoke renderers ignore it.
- Add each member as a hidden column, in member order, immediately after the composite, unless the
  collection already shows that member standalone (then it is simply visible, once). Apply the same
  columns-version / genome-build filter the visible columns get (`q_columns_this_version`) so a member
  only annotated for another build drops out of the hidden set as well - renderers must cope with an
  absent key, which the existing ones already do because they null-check every field. A composite
  whose members are *all* filtered out drops with them: there is nothing left for the cell to draw,
  and its sort keys and sort menu would name columns the row doesn't carry. The sort keys, the sort
  menu and `client_renderer_kwargs.members` are all built from the surviving members, so the headline
  is always a value the row actually has.

The cohort node's zygosity cell (`CohortNode._get_node_extra_columns`) is built per node from
node-specific count columns that are not in the catalogue. It stays in Python with its own `sort_menu`.
The catalogue composite for the **global** zygosity counts uses the same renderer.

### 4.2 Everything downstream that just works

- Export (`include_in_csv` on hidden members) - unchanged.
- The sort menu JS (`DataTableDefinition.setupSortMenus`) already resolves each entry by column
  `data`, which is the member's path.
- The column summary form (`ColumnSummaryForm.get_summarisable_columns`) offers keyed catalogue
  columns; the members are keyed, so summaries are per member, as now.
- Filter nodes, evidence keys, VEP column map, version diff - all key on member columns, untouched.

### 4.3 Deployment check

`variantgrid/deployment_validation/column_check.py`: replace the hardcoded `SPECIAL_COLUMNS` pseudo
paths with `exclude(queryset_field=False)`, keeping the annotation-alias entries. Add a check that
every composite has at least one member.

### 4.4 Renderers (`variantgrid_formats.js`, overrides in `snpdb/grids.py`)

Overrides are keyed by the composite's `variant_column` (equal to its pk for display-only composites,
the convention `classifications` already uses). Each carries `client_renderer` and `width`.

**Generic** - `VariantGridFormat.composite(data, type, rowData, meta)` reads `renderKwargs.members`:

- The first member's value is the cell, rendered through that member's own client renderer when it
  has one (cosmic id link, MaveDB URN link - resolved server side into the member entry as
  `renderer`), else as text. Choice fields arrive already expanded by the server renderer.
- Every other non-blank member goes in the tooltip as `label: value`. Blank when the whole group is
  blank, which is what makes sparse groups skimmable.
- A `.composite-cell` wrapper (ellipsis on overflow) plus `.composite-empty-headline` for a group
  with hover detail but no headline value - a muted `&middot;` so the cell stays hoverable. The
  bespoke cells keep the per-composite styles #686 gave them.

**Bespoke** (existing, re-keyed): `representativeVariant` (`variant`), `classifications`,
`impactConsequence`, `spliceai`, `maxentscan`, `mastermind`, `aloft`, `predictions`,
`dbZygosityCounts`. **Bespoke, extended**: `gnomadAf` becomes `gnomad`, drawing overall AF, the popmax
AF with its population tag and the filter flag inline, with populations / AC-AN / hom-hemi / FAF / XY on
hover. `gnomadPopmax` is retired into it. `predictions` gains the per-tool calls on hover from its new
detail members.

## 5. Arrange page (`view_custom_columns`)

`snpdb/views/views.py:view_custom_columns` and `snpdb/templates/snpdb/settings/view_custom_columns.html`:

- Query columns with `prefetch_related("composite_members__column")` and annotate whether each is a
  member, so the template does no per-row queries.
- **Available columns** lists non-members and composites. Members render into the same list with a
  `composite-member` class and are hidden by CSS until the "Show columns already inside a composite"
  checkbox (unchecked by default, above the list) toggles `show-members` on `#columns-widget`. No server
  round trip - the POST handler already accepts any catalogue id.
- **A composite card** shows its label, description and a compact member list ("shows: max delta score
  · acceptor gain · ...", from `member_columns` labels).
- **A member card** (in either list) carries a note "shown inside &lt;composite label&gt;".
- `page_help/settings/view_custom_columns_help.html`: a short paragraph on composites and the checkbox.

`global.scss` (+ hand-applied `global.css`): `.composite-column` badge on the card, `.composite-member`
hidden unless `#columns-widget.show-members`.

## 6. Migrations

`snpdb/migrations/0225`-`0230` are on origin and stay as they are. Everything below is new, on top of
`0236_open_targets_l2g_scores_column`.

### 6.1 `0237_compositecolumnmember` - schema

`CompositeColumnMember`.

### 6.2 Collapse helper - reusable by every later composite migration

`library/django_utils/composite_columns.py`, beside `bulk_insert_class_data` (migrations import helpers
that take `apps`, never live models):

```python
def collapse_into_composite(apps, composite_id: str) -> None:
    """ In every collection: replace the run of this composite's members with the composite, at the
        position of the first member (or the composite itself, if already present). A collection with
        none of them is untouched. Bumps version_id on each collection changed (historical models skip
        CustomColumn.save(), which is what normally bumps the node grid definition cache key) """
```

Semantics, per collection (`ordered` = its column ids by `sort_order`):

1. `group = [c for c in ordered if c == composite_id or c in members]`; nothing in group → skip.
2. `position = ordered.index(group[0])`; remove every id in `group`; insert `composite_id` at `position`.
3. Delete `CustomColumn` rows not in the new order, `get_or_create` the composite, renumber
   `sort_order` from 0, `version_id += 1`, save.

Idempotent. This is 0227's loop with the "All columns" exemption removed and the members read from
`CompositeColumnMember`. Members are read from the historical model, so a migration adds its member rows
first and then calls this.

### 6.3 `0238_composite_columns` - data

Steps:

1. **Insert the display-only composite rows** (`bulk_insert_class_data`; `model_field=False`,
   `queryset_field=False`, `variant_column` = pk, `annotation_level` per the table). `variant` and
   `classifications` already exist and just gain members.
2. **Restore each old anchor** to its pre-#686 label, width and description - it is a plain member
   again. Recover the originals from the migration that inserted the column
   (`git log -S"'grid_column_name': 'gnomad_af'" -- snpdb/migrations`, likewise consequence,
   gnomad_popmax_af, maxentscan_percent_diff_ref, mastermind_count_1_cdna, predictions_num_pathogenic,
   total_db_het, aloft_pred). `spliceai_max_ds` was born in 0227 with composite wording - give it a plain
   "Highest of the four SpliceAI delta scores" description and no width. The composite rows take over
   the #686 wording and widths.
3. **Insert `CompositeColumnMember` rows** per the table, first listed member = `sort_order` 0.
4. **Collapse every collection**: `for composite_id in COMPOSITES: collapse_into_composite(apps, composite_id)`.
5. Reverse: delete the new composite rows and all member rows (cascades to `CustomColumn`);
   collections keep whatever members they still had. Lossy - say so in the migration docstring.

The 0225 "leading columns" reorder (variant, classifications, tags, tags_global first) is left alone.

**Composites.** Members are catalogue pks in `sort_order`; the first is the headline (sort key, what
the generic renderer prints), the rest are hover detail. `sort=no` lists members with `in_sort_menu=False`. Renderer G = generic, B = bespoke.

Existing cells, now on the model:

| pk | level | label | members | sort=no | R |
|---|---|---|---|---|---|
| `variant` (exists, keyed `id`) | - | Variant | sorts on its own key; chrom, position, ref, alt, svlen, hgvs_c, hgvs_p, hgvs_g, symbol | all | B |
| `classifications` (exists) | D | Classifications | clinvar_highest_pathogenicity; clinvar_review_status, clinical_significance, conflicting_clinical_significance, clinvar_variation_id, clinvar_allele_id, clinvar_origin, clinvar_preferred_disease_name, clinvar_disease_database_name, clinvar_clinical_sources, drug_response, clinvar_highest_oncogenicity, clinvar_oncogenic_classification, clinvar_oncogenic_review_status, clinvar_somatic_tier, clinvar_somatic_review_status, max_internal_classification, internally_classified, internally_classified_labs, max_internal_somatic_classification, internally_classified_somatic | ids, disease names, sources, labs, conflicting | B |
| `consequence_impact` | V | Consequence | consequence; impact | - | B |
| `gnomad` | V | gnomAD | gnomad_af; gnomad_popmax, gnomad_filtered, gnomad_popmax_af, gnomad_ac, gnomad_an, gnomad_hom_alt, gnomad_hemi_count, gnomad_popmax_ac, gnomad_popmax_an, gnomad_popmax_hom_alt, gnomad_afr_af, gnomad_amr_af, gnomad_asj_af, gnomad_eas_af, gnomad_fin_af, gnomad_mid_af, gnomad_nfe_af, gnomad_oth_af, gnomad_sas_af, gnomad_xy_af, gnomad_xy_ac, gnomad_xy_an, gnomad_non_par, gnomad_faf95, gnomad_faf99, gnomad_fafmax_faf95_max, gnomad_fafmax_faf99_max, gnomad2_liftover_af | gnomad_popmax, gnomad_filtered, gnomad_non_par | B |
| `spliceai` | V | SpliceAI | spliceai_max_ds; spliceai_pred_ds_ag, _ds_al, _ds_dg, _ds_dl, _dp_ag, _dp_al, _dp_dg, _dp_dl, spliceai_gene_symbol | dp_*, gene_symbol | B |
| `maxentscan` | V | MaxEntScan | maxentscan_percent_diff_ref; maxentscan_ref, maxentscan_alt, maxentscan_diff | - | B |
| `mastermind` | V | Mastermind | mastermind_count_1_cdna; mastermind_count_2_cdna_prot, mastermind_count_3_aa_change, mastermind_mmid3 | mmid3 | B |
| `aloft` | V | ALoFT | aloft_pred; aloft_prob_dominant, aloft_prob_recessive, aloft_prob_tolerant, aloft_high_confidence, aloft_ensembl_transcript | transcript | B |
| `predictions` | V | Predictions | predictions_num_pathogenic; predictions_num_benign, sift, polyphen2_hvar_pred_most_damaging, mutation_taster_pred_most_damaging, mutation_assessor_pred_most_damaging, fathmm_pred_most_damaging, metalr_rankscore | the per-tool calls | B |
| `db_zygosity` | D | Database counts | total_db_het; total_db_hom, total_db_ref, total_db_unk | - | B |

New cells:

| pk | level | label | members | sort=no | R |
|---|---|---|---|---|---|
| `alphamissense` | V | AlphaMissense | alphamissense_pred; alphamissense_score, alphamissense_rankscore | - | G |
| `bayesdel` | V | BayesDel | bayesdel_noaf_score; bayesdel_noaf_rankscore | - | G |
| `cadd` | V | CADD | cadd_phred; cadd_raw, cadd_raw_rankscore | - | G |
| `clinpred` | V | ClinPred | clinpred_pred; clinpred_score, clinpred_rankscore | - | G |
| `eve` | V | EVE | eve_class; eve_score | - | G |
| `metarnn` | V | MetaRNN | metarnn_pred; metarnn_score | - | G |
| `primateai` | V | PrimateAI | primateai_pred; primateai_score | - | G |
| `revel` | V | REVEL | revel_score; revel_rankscore | - | G |
| `vest4` | V | VEST4 | vest4_score; vest4_rankscore | - | G |
| `varity` | V | VARITY | varity_r_score; varity_er_score | - | G |
| `mutpred2` | V | MutPred2 | mutpred2_score; mutpred2_top5_mechanisms | mechanisms | G |
| `hipred` | G | HIPred | hipred_prediction; hipred_score | - | G |
| `gene_damage_index` | G | Gene damage index | gene_damage_index_phred; gene_damage_index_score | - | G |
| `gene_indispensability` | G | Gene indispensability | gene_indispensability_pred; gene_indispensability_score | - | G |
| `gnomad_constraint` | G | gnomAD constraint | gnomad_pli; gnomad_oe_lof, gnomad_pnull, gnomad_prec | - | G |
| `gnomad_sv` | V | gnomAD SV overlap | gnomad_sv_overlap_af; gnomad_sv_overlap_percent, gnomad_sv_overlap_name, gnomad_sv_overlap_coords | name, coords | G |
| `conservation` | V | Conservation | phylop_100_way_vertebrate; phylop_30_way_mammalian, phylop_46_way_mammalian, phastcons_100_way_vertebrate, phastcons_30_way_mammalian, phastcons_46_way_mammalian, gerp_pp_rs | - | G |
| `cosmic` | V | COSMIC | cosmic_id (link); cosmic_count, cosmic_legacy_id | ids | G |
| `denovo_db` | V | denovo-db | denovo_db_case_count; denovo_db_control_count, denovo_db_studies, denovo_db_primary_phenotypes, denovo_db_pubmed_ids | text | G |
| `open_targets` | V | Open Targets | open_targets_gwas_l2g_score; open_targets_gwas_diseases, open_targets_gwas_gene_id, open_targets_qtl_biosample, open_targets_qtl_gene_id, open_targets_study_id, open_targets_study_type, open_targets_variant_id | all but headline | G |
| `promoter_ai` | V | PromoterAI | promoter_ai_score; promoter_ai_tss_pos | - | G |
| `mavedb` | V | MaveDB | mavedb_score; mavedb_urn (link) | urn | G |
| `dbscsnv` | V | dbscSNV | dbscsnv_ada_score; dbscsnv_rf_score | - | G |
| `nmd_ptc` | T | NMD / PTC | nmd_escape_status; nmd_escaping_variant, ptc_distance_codons, ptc_last_junction_distance | - | G |
| `protvar` | V | ProtVar | protvar_stability; protvar_int, protvar_pocket | - | G |
| `essential_gene` | G | Essential gene | essential_gene_crispr; essential_gene_crispr2, essential_gene_gene_trap | - | G |
| `annotsv_acmg` | V | AnnotSV ACMG | annotsv_acmg_class; annotsv_acmg_score, annotsv_acmg_criteria | criteria | G |
| `annotsv_benign_af` | V | AnnotSV benign AF | annotsv_b_loss_af_max; annotsv_b_gain_af_max, annotsv_b_ins_af_max, annotsv_b_inv_af_max | - | G |
| `annotsv_splice` | V | AnnotSV splice site | annotsv_dist_nearest_ss; annotsv_nearest_ss_type | type | G |
| `annotsv_omim` | V | AnnotSV OMIM | annotsv_omim_morbid; annotsv_omim_id, annotsv_omim_inheritance, annotsv_omim_phenotype | all but headline | G |
| `annotsv_region` | V | AnnotSV region | annotsv_repeat_type_left; annotsv_repeat_type_right, annotsv_segdup_left, annotsv_segdup_right, annotsv_encode_blacklist_left, annotsv_encode_blacklist_right, annotsv_encode_blacklist_characteristics_left, annotsv_encode_blacklist_characteristics_right | all | G |

Open Targets overlaps the #1822 plan (`1822_open_targets_variant_details_plan.md`), which landed first -
`open_targets_gwas_l2g_scores` (the per-record scores, snpdb migration 0236) joins the composite as a
detail member.

`alphamissense_rankscore` carried `annotation_level` `'A'`, which is not a `ColumnAnnotationLevel` - the
arrange page raises on it. The data migration corrects it to transcript level, matching its siblings.

Left standalone on purpose: af_1kg, af_uk10k, topmed_af (different sources, read side by side with
gnomAD), exon/intron (already one server renderer), domains/interpro, grantham, repeat_masker,
splice_region, variant_class, codons/amino_acids, annotsv_pathogenic_overlaps / exons_spanned /
frameshift / re_gene, tags, Sample, and every identifier / gene-text column.

## 7. Code removal

- `COMPOSITE_COLUMN_ROW_FIELDS`, `VARIANT_COLUMN_ROW_FIELDS`, `CLASSIFICATIONS_COLUMN_ROW_FIELDS` and
  the hidden-field blocks in `get_variant_grid_columns` that read them.
- The nine `sort_menu` literals in `_get_standard_overrides` (`snpdb/grids.py`). Keep each
  `client_renderer` and `width` entry, re-keyed from the anchor path to the composite pk. Standalone
  member renderers (`gnomad_filtered`, `cosmicLink`, `mavedbUrn`, `masterMind`) stay - a member shown
  standalone keeps its link, and the generic renderer reuses them for headlines.
- `VariantGridFormat.gnomadPopmax` (folded into `gnomad`).

## 8. Tests

Update:

- `analysis/tests/test_representative_variant_column.py` - partner assertions read members from the
  model; the sort-menu test checks the derived menu names grid columns.
- `snpdb/tests/test_custom_columns_genome_build.py` - add: a composite whose member is build-specific
  (`gnomad_faf95`, GRCh38 only) drops that member from the hidden set for GRCh37 and still renders.

New:

- `snpdb/tests/test_composite_columns.py`
  - composite in a collection → visible composite RichColumn (no key, sort_keys = first member, derived
    sort_menu honours `in_sort_menu`, `client_renderer_kwargs.members` in member order) followed by its
    members hidden with `include_in_csv`.
  - keyed composite (`variant`) keeps its key and gets its members hidden.
  - member already visible in the collection is emitted once, visible.
  - `collapse_into_composite` (run against `django.apps.apps`): anchor+partners run → composite at the
    run's position; lone partner → composite; composite already present with stray partners → partners
    gone, position kept; collection without any → untouched; `version_id` bumped only when changed.
  - every composite in the catalogue has members and every member is a real column (guards the
    migration table).
- Arrange page (existing view test or `snpdb/tests/test_urls.py`) - lists composites, members carry
  the `composite-member` class.

Audit at the end per CLAUDE.md: keep the collapse cases, the grid column construction and the catalogue
invariant; drop anything that only restates a field declaration.

## 9. Later additions

A new composite is: a mockup if bespoke, a catalogue row plus member rows, an override entry (generic
renderer needs none), and a migration calling `collapse_into_composite`. Candidates deliberately left
out (§6.3 "left standalone") can be revisited with users; conservation is the one most likely to be
argued the other way.

## 10. Files

| File | Change |
|---|---|
| `snpdb/models/models_columns.py` | `CompositeColumnMember`, `VariantGridColumn` helpers, css classes |
| `snpdb/migrations/0237_compositecolumnmember.py` | schema |
| `library/django_utils/composite_columns.py` | `collapse_into_composite(apps, composite_id)` |
| `snpdb/migrations/0238_composite_columns.py` | rows, restores, members, collapse |
| `snpdb/grid_columns/custom_columns.py` | composite → RichColumn + hidden members; lists removed |
| `snpdb/grids.py` | overrides re-keyed to composite pks, sort_menu literals removed, generic renderer wired |
| `variantgrid/static_files/default_static/js/variantgrid_formats.js` | `composite` generic renderer, `gnomad` extended, `predictions` hover |
| `snpdb/views/views.py` | `view_custom_columns` context |
| `snpdb/templates/snpdb/settings/view_custom_columns.html` | composite cards, member toggle |
| `variantgrid/static_files/default_static/page_help/settings/view_custom_columns_help.html` | composites paragraph |
| `variantgrid/static_files/default_static/css/global.scss` + `.css` | card styles |
| `variantgrid/deployment_validation/column_check.py` | display-only exclusion, member check |
| tests per §8 | |
