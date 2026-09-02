# Annotation descriptions page: group by composite, show the cell (variantgrid_private#3905)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-02

## Goal

`annotation/view_annotation_descriptions` answers "what does this grid column mean and where did it
come from". Since #3904 the grid draws 41 composite cells (`gnomad`, `spliceai`, `cadd`, ...) whose
members ride along hidden, and the page has not caught up:

- Every composite falls into the view's "not VEP, not AnnotSV, so VariantGrid calculated it during
  annotation import" bucket. `gnomad` is listed as a VariantGrid-calculated variant-level annotation,
  alphabetically among its own members, with nothing tying them together.
- Rows are headed by `grid_column_name` (`gnomad_af`) while the grid shows labels (`gnomAD AF`).
- The five columns with no annotation level (`chrom`, `position`, `ref`, `alt`, `svlen`) are never
  rendered at all - the template only iterates the levels it knows.
- A composite whose members all drop out for the build (AnnotSV under a build without it, gnomAD SV)
  vanishes from the grid but still lists here.

The page becomes: a **Composite cells** card at the top, one section per composite showing the cell
exactly as the grid draws it (header with its ▾ sort menu, then the cell, hover title included)
alongside its description and a table of the columns it draws; then the existing per-level cards
listing the remaining plain columns. The example cells are rendered by the real client renderers
from hand-written, fictional row data, so they never rot and never leak real data.

## Data

### Example rows - `snpdb/grid_columns/composite_examples.py`

The only new data is fictional and lives in code. One row per composite, keyed by composite pk,
holding the values the **client** row carries, keyed by member `variant_column` path:

```python
# Composite pk -> the grid row the cell draws, as the client receives it: keyed by each member's
# variant_column path. Allele frequencies are unit values - the view formats them the way the grid
# does (@see AF_UNIT_COLUMNS) so the example follows the deployment's percent setting.
COMPOSITE_EXAMPLE_ROWS: dict[str, dict[str, Any]]

# The sample call cell is built per sample by the analysis grid rather than catalogued, so its
# example carries its own column: the row keyed by SAMPLE_EXAMPLE_PREFIX + format column
SAMPLE_EXAMPLE_PREFIX = "sample_1_"
SAMPLE_EXAMPLE_ROW: dict[str, Any]
```

Conventions for the values (the full set is in the appendix):

- **Fictional throughout.** Gene `EXAMPL1`, transcript `NM_012345.6` / `ENST00000123456`, chromosome
  `12`, position `12345678`, sample `NA00001`, lab `Example Lab`, disease `Example syndrome`. Numbers
  are plausible and internally consistent (popmax AF ≥ overall AF; `predictions_num_pathogenic` equals
  the number of damaging tool calls listed; MaxEntScan % diff = 100 × diff / ref; rankscore rises with
  score).
- **Every member has a value**, so the hover title shows the whole group. A headline is never blank.
- **Values are what the client sees.** Where the server renders a choice field, the row holds the
  label the grid shows, not the code. The classification chips are the exception: they key on raw codes
  (`clinvar__highest_pathogenicity` 4, `max_internal_classification` `'3'`, `clinvar__review_status` a
  `ClinVarReviewStatus` code), because that is what the grid row carries for them. Check each bespoke
  renderer in `variantgrid_formats.js` for what it reads.
- **AFs are unit values** and the view formats them; see "AF formatting" below.

No model changes. No migration.

## Page layout

```
VariantGrid Annotations                                      [GRCh38 ▾]
 (page help)  · customise columns · search example

┌ Composite cells ─────────────────────────────────────────────────────────────┐
│ A composite draws several columns in the one cell. The first is the headline  │
│ - what the cell reads as and sorts by - the rest are on hover, and the ▾ in   │
│ the header sorts by any of them. Members are still exported and filterable.  │
│                                                                              │
│  ┌──────────────┐   gnomAD                                    [Variant level]│
│  │ gnomAD     ▾ │   Allele Frequency (0-1) among all gnomAD genotypes ...    │
│  ├──────────────┤   Sorts by gnomAD AF. Hover the cell for the detail.       │
│  │0.012% 0.051% │                                                            │
│  │  NFE   [Pass]│   Label            Column            Role     Level Source │
│  └──────────────┘   gnomAD AF        gnomad_af         headline  V   VEP ... │
│                     gnomAD PopMax AF gnomad_popmax_af  sort      V   VEP ... │
│                     gnomAD PopMax    gnomad_popmax     hover     V   VEP ... │
│                     ...                                                      │
│  ──────────────────────────────────────────────────────────────────────────  │
│  (next composite)                                                            │
└──────────────────────────────────────────────────────────────────────────────┘

┌ Database  Internal database stats ───────────────────────────────────────────┐
│  (plain columns only - members now live in their composite's section)        │
└──────────────────────────────────────────────────────────────────────────────┘
┌ Transcript level ... ┐  ┌ Variant level ... ┐  ┌ HGNC ┐ ┌ UniProt ┐ ┌ ClinVar ┐ ┌ Gene ┐
```

Composite sections are ordered by the composite's annotation level in the page's existing order
(none/variant identity, sample, Database, Transcript, Variant, Gene), then label. The sample call
cell is a section like any other, replacing the hardcoded "Sample Level" card: its member table is
that card's existing hand-written rows (Zygosity, AD, AF, DP, GQ, PL) plus a row for FT (sample
filters), which the cell also draws.

A composite's level and its members' levels do not always agree (`alphamissense` is variant level,
its members transcript level; `cadd` mixes both). Each member row therefore carries its own level
colour box, and the section's level chip is the composite's own.

The member table reuses the per-level tables' columns (Source / Category / Source field / Details)
with three more: label (bold), `grid_column_name` (muted), role. Role is `headline` for
`sort_order == 0`, `sort` for other members with `in_sort_menu`, `hover` for the rest.

## The example cell

Each section holds a one-column table with the grid's own classes, so it picks up every cell style
from `global.scss`'s `table.variantgrid-datatable` block and DataTables' header styling:

```html
<table class="variantgrid-datatable dataTable composite-example" data-column-json="..." data-row-json="...">
  <thead><tr><th class="sorting" style="width: 130px" title="{description}">gnomAD</th></tr></thead>
  <tbody><tr><td></td></tr></tbody>
</table>
```

Page JS (inline `{% block head %}` script, as the arrange page does) fills each one the way
`DataTableDefinition` would, mirroring `datatable_definition.js` lines ~314-325 and
`setupSortMenus`:

```js
const col = table.data("columnJson"), row = table.data("rowJson"), extra = EXAMPLE_EXTRA;
const html = eval(col.render)(row[col.data], "display", row, {extra: extra, kwargs: col.renderKwargs || null});
td.html(html instanceof jQuery ? html.prop("outerHTML") : html);
// header ▾: a .dt-sort-menu with one .dropdown-item per col.sortMenu entry, opened via FloatingPanel;
// picking an entry just closes the menu
// the cell is a picture: every <a> inside it loses its href/target so the fake variant links nowhere
```

`FloatingPanel` (`global.js`), `escapeHtml`, `VariantGridFormat` and DataTables' CSS are all loaded by
`base.html`, so nothing new is included. `inAnalysis()` is false here, so the Variant cell draws no
checkbox, matching the non-analysis grids. Cells render single-line (the default; the two-line-rows
user setting is not applied).

### Column JSON

The client needs the same column definition the grid sends. Build it through the real code path:

1. `snpdb/grid_columns/custom_columns.py`: expose the two private builders as one public function,
   `composite_rich_column(column, members, column_overrides) -> RichColumn`, which applies
   `_catalogue_column_kwargs` then `_composite_column_kwargs` and any override for the composite's
   own path, exactly as `get_variant_grid_columns` does. `get_variant_grid_columns` calls it too.
2. `snpdb/grids.py`: `AbstractVariantGrid._get_standard_overrides` uses nothing from `self` - make
   it a module-level `get_standard_overrides(af_show_in_percent) -> dict`, with the AF column list
   hoisted to a module constant `AF_UNIT_COLUMNS` (see AF formatting). The grid keeps calling it.
3. `snpdb/views/datatable_view.py`: extract the per-column dict built in the definition endpoint
   (`data`, `label`, `render`, `renderKwargs`, `sortMenu`, `width`, `headerTitle`, ...) into
   `rich_column_json(rc, default_column_width) -> JsonObjType`, used by both the endpoint and this view.
4. `snpdb/grids.py`: the body of `AbstractVariantGrid.get_extra` becomes module-level
   `variant_grid_client_extra(genome_build) -> JsonObjType` (`genomeBuild`, `clinvarStars`,
   `genotypeQuality`, `commonGnomadAf`); the grid's `get_extra` returns it. The page passes it to every
   example as `EXAMPLE_EXTRA`.

The sample call column is built the same way the analysis grid builds it - label
`"NA00001 Zygosity"`, `client_renderer` `VariantGridFormat.sampleZygosity`, `samplePrefix`
`sample_1_`, `sort_menu` from the sample sort key labels. Move `SAMPLE_SORT_KEY_LABELS` and
`SAMPLE_COMPOSITE_COLUMNS` out of `VariantGrid._get_grid_genotype_columns` in `analysis/grids.py`
into module constants in `snpdb/grid_columns/grid_sample_columns.py` (which `analysis/grids.py`
already imports from) so the view reads the same labels.

### AF formatting

The grid formats unit AFs server side so CSV and grid agree, and the client renderers only add
markup. The example rows store unit AFs; the view runs each key in `AF_UNIT_COLUMNS` (and the sample
row's `samples_allele_frequency`) through `get_allele_frequency_formatter(source_in_percent=False,
dest_in_percent=settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT)` via a `CellData`, the way
`NodeColumnSummary.get_dataframe` does. `commonGnomadAf` in the extra is in the same units, so the
gnomAD cell's common-variant muting behaves as in the grid.

## View - `annotation/views.py: view_annotation_descriptions`

- Read the catalogue once with `prefetch_related("composite_members__column")` and call
  `VariantGridColumn.annotate_composite_membership`, as the arrange view does.
- Keep `_first_for_build` / the VEP-AnnotSV-VariantGrid row classification for each column. Factor
  the per-column row dict into `_column_row(vgc)` returning `None` for a VEP column not populated in
  this build, so composites and plain columns share it.
- For each composite: members = those whose `_column_row` is not `None`. A composite with none is
  dropped (the grid's rule). Build the section:
  `{column, level, rows: [member row + role], column_json, row_json, sort_by: headline label}`
  where `column_json` comes from `composite_rich_column` on the surviving members and
  `get_standard_overrides(...)`, serialised with `rich_column_json`, and `row_json` is the example row
  after AF formatting. A composite missing from `COMPOSITE_EXAMPLE_ROWS` raises - the test below
  makes that a test failure, not a page error.
- Plain columns (not composite, not a member) go to the per-level dicts as now.
- Context gains `composite_sections`, `sample_section`, `example_extra`.
- The `annotation-source` badge for a composite section reads "Composite"; add the CSS class.

`@cache_page(WEEK_SECS)` stays.

## Template

`annotation/templates/annotation/view_annotation_descriptions.html`:

- Add the Composite cells card above the level cards, with the intro paragraph from the layout above
  and one section per composite (and the sample section first among the sample-level ones). Give
  each section `id="composite-{{ pk }}"` so the arrange page and changelog can deep link.
- Delete the hardcoded Sample Level card; its rows move into the sample section's member table.
- Level cards are unchanged except they now only receive plain columns. Show label bold with
  `grid_column_name` muted beneath, matching the member tables.
- Move the page's `<style>` block additions (`.composite-example` sizing, section divider, role
  chips) into `global.scss` under a `.annotation-descriptions` wrapper, and hand-apply to `global.css`
  per the SCSS rule in CLAUDE.md.

## Page help

`variantgrid/static_files/default_static/page_help/annotation/view_annotation_descriptions_help.html`:
add a short "Composite cells" paragraph (headline / hover / ▾ sort menu, members still exported and
filterable, the example cells are fictional) and link to the arrange page.

## Tests

`snpdb/tests/test_composite_examples.py`:
- Every composite in the catalogue has an example row, and every key in it is the `variant_column`
  of one of that composite's members (or `id` for `variant`). Catches drift when a composite gains or
  loses members.
- Each example's headline member value is non-blank.
- Every `SAMPLE_EXAMPLE_ROW` key is `SAMPLE_EXAMPLE_PREFIX` + a sample format column.

`annotation/tests/test_urls.py` already asserts 200. Add to `annotation/tests/test_views.py` (or the
nearest existing view test module): the rendered page contains a section per composite that has
members for the build, and the composite sections' `column_json` names the members' `sort_menu`
entries - one assertion on `gnomad` is enough.

Audit at the end: keep only tests covering our rules (the drift check, the per-build drop).

## Appendix - the example values

Bespoke-renderer composites, complete (paths abbreviated: `va__` = `variantannotation__`):

**variant** (renderer `representativeVariant`, value `id`)
`id 1 · locus__contig__name "12" · locus__position 12345678 · locus__ref__seq "C" · alt__seq "T" · svlen null ·
va__hgvs_c "NM_012345.6(EXAMPL1):c.1234C>T" · va__hgvs_p "NP_012345.1:p.(Arg412Trp)" ·
va__hgvs_g "NC_000012.12:g.12345678C>T" · va__symbol "EXAMPL1"`

**sample call** (`sampleZygosity`, prefix `sample_1_`)
`samples_zygosity "HET" · samples_allele_frequency 0.48 · samples_allele_depth 23 · samples_read_depth 48 ·
samples_genotype_quality 99 · samples_phred_likelihood 0 · samples_filters "PASS"`

**classifications** (`classifications`)
ClinVar germline LP with 2 stars, internal germline VUS from one lab, ClinVar somatic Tier III, no
internal somatic: `clinvar__highest_pathogenicity 4 · clinvar__review_status <2-star code> ·
clinvar__clinical_significance "Likely pathogenic" · conflicting false · variation_id 123456 ·
allele_id 654321 · origin "germline" · preferred_disease_name "Example syndrome" ·
disease_database_name "MONDO:0000001" · clinical_sources "Example Lab" · drug_response null ·
highest_oncogenicity null · oncogenic_classification null · oncogenic_review_status null ·
somatic_tier "tier_3" · somatic_review_status <1-star code> · max_internal_classification "3" ·
internally_classified "1" · internally_classified_labs "Example Lab" · max_internal_somatic_classification null ·
internally_classified_somatic null`

**consequence_impact** (`impactConsequence`): `consequence "missense_variant" · impact "MODERATE"`

**gnomad** (`gnomad`): `gnomad_af 0.000123 · gnomad_popmax "NFE" · gnomad_filtered false ·
gnomad_popmax_af 0.00051 · ac 31 · an 251432 · hom_alt 0 · hemi_count 0 · popmax_ac 29 · popmax_an 56780 ·
popmax_hom_alt 0 · afr 0.00002 · amr 0.00009 · asj 0 · eas 0 · fin 0.00004 · mid 0 · nfe 0.00051 · oth 0.0001 ·
sas 0.00001 · xy_af 0.00012 · xy_ac 15 · xy_an 124800 · non_par true · faf95 0.00009 · faf99 0.00008 ·
fafmax_faf95_max 0.00038 · fafmax_faf99_max 0.00034 · gnomad2_liftover_af 0.00011`

**spliceai** (`spliceai`): `max_ds 0.62 · ds_ag 0.02 · ds_al 0.62 · ds_dg 0.00 · ds_dl 0.01 ·
dp_ag -12 · dp_al 3 · dp_dg 41 · dp_dl -7 · gene_symbol "EXAMPL1"`

**maxentscan** (`maxentscan`): `percent_diff_ref -38.4 · ref 9.21 · alt 5.67 · diff -3.54`

**mastermind** (`mastermind`): `count_1_cdna 7 · count_2_cdna_prot 9 · count_3_aa_change 14 · mmid3 "EXAMPL1:R412W"`

**aloft** (`aloft`): `aloft_pred "Recessive" · prob_dominant 0.12 · prob_recessive 0.81 · prob_tolerant 0.07 ·
high_confidence true · ensembl_transcript "ENST00000123456"`

**predictions** (`predictions`): `num_pathogenic 4 · num_benign 2 · sift "Deleterious" ·
polyphen2_hvar "Probably damaging" · mutation_taster "Disease causing" · mutation_assessor "Low" ·
fathmm "Tolerated" · metalr_rankscore 0.91` (4 damaging: SIFT, Polyphen2, MutationTaster, MetaLR)

**db_zygosity** (`dbZygosityCounts`, keys `global_variant_zygosity__*`): `het_count 3 · hom_count 1 · ref_count 42 · unk_count 0`

Generic-renderer composites - headline first, then hover detail; the implementer fills any member
not named here in the same spirit:

| Composite | Headline | Detail |
|---|---|---|
| alphamissense | pred "likely_pathogenic" | score 0.87, rankscore 0.92 |
| bayesdel | noaf_score 0.31 | rankscore 0.88 |
| cadd | phred 24.3 | raw 3.91, rankscore 0.93 |
| clinpred | pred "Damaging" | score 0.96, rankscore 0.94 |
| eve | class "Pathogenic" | score 0.81 |
| metarnn | pred "Damaging" | score 0.89 |
| primateai | pred "Damaging" | score 0.83 |
| revel | score 0.78 | rankscore 0.90 |
| vest4 | score 0.86 | rankscore 0.91 |
| varity | r_score 0.74 | er_score 0.69 |
| mutpred2 | score 0.71 | top5_mechanisms "Loss of helix (P = 0.02); Gain of loop (P = 0.04); ..." |
| hipred | prediction "Y" (haploinsufficient) | score 0.72 |
| gene_damage_index | phred 3.2 | score 0.41 |
| gene_indispensability | pred "E" (essential) | score 0.88 |
| gnomad_constraint | pli 0.97 | oe_lof 0.12, pnull 0.01, prec 0.02 |
| gnomad_sv | overlap_af 0.0021 | overlap_percent 94.5, name "gnomAD-SV_v3_DEL_12_00001", coords "12:12300000-12400000" |
| conservation | phylop_100_way 7.85 | phylop_30 1.21, phylop_46 2.03, phastcons_100 1.0, phastcons_30 0.99, phastcons_46 0.98, gerp_pp_rs 5.63 |
| cosmic | cosmic_id "COSV100000001" | count 3, legacy_id "COSM1000001" |
| denovo_db | case_count 2 | control_count 0, studies "Example2021,Example2023", primary_phenotypes "developmental disorder", pubmed_ids "10000001,10000002" |
| open_targets | gwas_l2g_score 0.63 | diseases "Example trait", gwas_gene_id "ENSG00000123456", qtl_biosample "whole blood", qtl_gene_id "ENSG00000123456", study_id "GCST000001", study_type "gwas", variant_id "12_12345678_C_T", l2g_scores "ENSG00000123456:0.63" |
| promoter_ai | score 0.42 | tss_pos -85 |
| mavedb | score -1.24 | urn "urn:mavedb:00000001-a-1" |
| dbscsnv | ada_score 0.97 | rf_score 0.88 |
| nmd_ptc | escape_status "escaping" | escaping_variant true, ptc_distance_codons 12, ptc_last_junction_distance -34 |
| protvar | stability -1.8 | int "yes", pocket "no" |
| essential_gene | crispr "E" | crispr2 "E", gene_trap "N" |
| annotsv_acmg | class 4 | score 0.9, criteria "1A,2C,5F" |
| annotsv_benign_af | b_loss_af_max 0.012 | b_gain 0.004, b_ins 0, b_inv 0 |
| annotsv_omim | omim_morbid "yes" | omim_id "600001", omim_inheritance "AD", omim_phenotype "Example syndrome" |
| annotsv_region | repeat_type_left "AluY" | repeat_type_right "L1", segdup_left "chr12:12000000-12010000", segdup_right null, encode_blacklist_left null, encode_blacklist_right "chr12:12400000-12410000", characteristics_left null, characteristics_right "High Signal Region" |
| annotsv_splice | dist_nearest_ss 42 | nearest_ss_type "5'" |

Where a member is a choice field, look at the field's choices and use the display label the server
renders (the row above shows the intent, e.g. HIPred "Y" may render as "Yes").
