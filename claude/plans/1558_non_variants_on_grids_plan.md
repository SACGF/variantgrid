# #1558 — Showing non-variants (SV / CNV / fusions) on the variant grids

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-07
Status: in progress

[#1558](https://github.com/SACGF/variantgrid/issues/1558): SVs, CNVs and gene fusions are all `Variant` rows
now, but the grids show them as if they were small variants. This is Phase 6 of
`claude/plans/tso500_overall_plan.md` (the #1558 half; #1706's specimen/extraction grids stay there).

The single-vs-multiple grid question is settled: **one grid**. A fusion or CNV is a row like any other so
compound-het, gene lists and every downstream node keep working over one result set - that is what storing
them as `Variant` was for (`snpdb/gene_level_variants.py`). What this plan adds is the presentation that
makes the kind obvious, the per-call data that only these rows have, and the places fusions currently fall
out of.

---

## 1. Where things stand (verified on vg-test2, 2026-09-07)

| Kind | Stored as | Grid today |
|---|---|---|
| Small variant | explicit ref/alt | Representative cell: gene + c.HGVS, or g.HGVS, or coordinate |
| SV | `<DEL>` / `<DUP>` / `<INV>` with `svlen`, VEP class `deletion` / `duplication` / `inversion` | g.HGVS when short (`chr1 g.9770512_9787106dup`), else coordinate + alt + size. No size on the HGVS form; nothing says "structural" |
| CNV | Same `Variant` as an SV, plus the caller's copy number in `CohortGenotype` | Copy number not shown anywhere |
| Fusion | `Variant` on the `GENE_LEVEL` contig, alt `<FUSION:HGNC:nnn>`, `GeneFusion` 1:1 (#1506) | `hgvs_c` = `BCR::ABL1` so the cell shows the pair - indistinguishable from a gene symbol; calls (breakpoints, caller) invisible |

Three shapes of "copy number" exist in the database, which decides §4:

| Source | Where | Example (from `CohortGenotype`) |
|---|---|---|
| DRAGEN germline CNV (`dragen.cnv.vcf`) | `FORMAT/CN` integer, `FORMAT/SM` ratio | `format[0] = {CN: [3], SM: [1.498], GT: ./1}` |
| DRAGEN TSO 500 somatic (`ExampleSample_DNA_2600000001C.cnv.vcf`) | `FORMAT/SM` only, no CN | `SM 1.494` for a `<DUP>` |
| Pisces TSO 500 legacy (`TSO500_cnv.vcf`) | `INFO/CN` integer, `FORMAT/FC` fold change | `info = {CN: 2, ...}`, `format[0] = {FC: [1.125]}` |

Both JSON blobs are stored by importer v21+ for every FORMAT/INFO key the header declares
(`upload/vcf/bulk_genotype_vcf_processor.py`), so nothing below needs a re-import.

Already in place and reused rather than rebuilt:

- The representative Variant cell is a composite (`snpdb/grids.py:get_standard_overrides`, key `'id'`) whose
  members ride along hidden on every row: contig name, position, ref, alt, `svlen`, `hgvs_c/p/g`, symbol.
  Everything the kind badge needs is already on the row.
- The analysis "Variant type" filter is `analysis/models/nodes/filters/damage_node.py:DamageNode` over VEP's
  `variant_class`, grouped by `library/genomics/vcf_enums.py:VARIANT_CLASS_GROUPS` - which already has a
  Fusion group (`GENE_FUSION`, written by `annotation/gene_level_annotation.py`).
- The All Variants page has a variant-type button group (`snpdb/variant_filters.py:VariantType`).
- Gene pages find fusions from either partner through `VariantGeneOverlap`
  (`snpdb/variant_queries.py:get_variant_queryset_for_gene_symbol`); the gene-level annotation run writes
  overlap rows for both genes.
- Search resolves a gene pair (`BCR::ABL1`, `CD74-ROS1`) to the fusion, lookup only.

Not built on this box: zero `GeneFusion` rows. The loader and its test data are in the repo
(`upload/tests/test_import_dragen_tso500_all_fusions.py`).

---

## 2. Data

No new models. One nullable text column on `VCF`, following the existing "which FORMAT field means X"
columns beside it, and one choice column on `DamageNode`.

```python
class VCF(...):                                   # snpdb/models/models_vcf.py
    ...
    genotype_quality_field = models.TextField(null=True)
    phred_likelihood_field = models.TextField(null=True)
    sample_filters_field = models.TextField(null=True)
    # The FORMAT (or, for a single-sample VCF, INFO) key carrying the sample's copy number or copy
    # ratio - CN, SM, FC... The grid labels the value with this name and the header's description
    copy_number_field = models.TextField(null=True)


class VCFSourceSettings(models.Model):            # snpdb/models/models_vcf.py
    OVERRIDABLE_SAMPLE_FIELDS = frozenset({..., "copy_number_field"})   # one more entry, no column


class DamageNode(AnalysisNode):                   # analysis/models/nodes/filters/damage_node.py
    variant_class = ArrayField(...)               # unchanged
    variant_class_exclude = models.BooleanField(default=False)
    # Structural = symbolic alt with an SVLEN (SV, CNV, and gene-level rows). Orthogonal to
    # variant_class because VEP calls a 1 Mb <DEL> and a 1 bp deletion the same class
    structural = models.CharField(max_length=1, choices=StructuralFilter.choices, default=StructuralFilter.ANY)


class StructuralFilter(models.TextChoices):       # analysis/models/nodes/filters/damage_node.py
    ANY = 'A', "Any"
    ONLY = 'O', "Structural only"
    EXCLUDE = 'E', "Exclude structural"
```

Everything else this plan shows is read at query time from what is already stored:

| Shown | Read from | Rows it applies to |
|---|---|---|
| Kind badge | `alt__seq`, `svlen` (members already on the row) | all grids |
| Copy number | `cohortgenotype_<pk>__format__<i>__<field>__0`, falling back to `__info__<field>` | analysis grid, per sample |
| Fusion calls | `cohortgenotype_<pk>__info__FUSION_OBS` (percent-encoded JSON list, one per caller row) | analysis grid, per cohort |
| Fusion partners / direction | `GeneFusion` via the variant | row expansion, all grids |

---

## 3. Phase 1 — the kind badge, and fusions reaching every grid

### 3.1 Badge in the representative Variant cell

One client-side function, `_variantKind(rowData)` in
`variantgrid/static_files/default_static/js/variantgrid_formats.js`, next to `_representativeVariantLabel`,
returns `{code, label, title}` or null:

| Row | Badge | Title |
|---|---|---|
| alt matches `<FUSION:…>` | `FUSION` | `Gene fusion` |
| alt matches `<FUSION_UNORDERED:…>` | `FUSION` with a `⇄` | `Gene fusion - direction not asserted by the caller` |
| alt `<AMP:…>` / `<LOSS:…>` (designed, not yet loaded) | `AMP` / `LOSS` | `Gene-level copy number` |
| alt `<DEL>` / `<DUP>` / `<INV>` / `<CNV>` / `<INS>` | `DEL 1.1 kb` etc, size from `svlen` | `chr:start-end <ALT>` |
| explicit ref/alt | none | |

Small variants are the overwhelming majority of rows, so they carry no badge - the badge is the signal that a
row is not one. The label the cell already draws stays as it is (`BCR::ABL1`, `chr1 g.9770512_9787106dup`,
`chr8:1000-2000 DEL`); the badge sits after it as `<span class='rv-kind rv-kind-del'>`. In two-line rows it
stays on line one - the kind is the point. Branch 3 of `_representativeVariantLabel` (coordinate for a
symbolic alt) stops repeating the type and size, since the badge now carries them for every branch.

The alt-kind vocabulary is the one `library/genomics/vcf_enums.py:GeneLevelSymbolicAlt` and
`VCFSymbolicAllele` define; the JS keeps a small table of the same strings (a gene-level alt is
`<KIND:NAMESPACE:id>` or `<KIND:UNKNOWN>`, so the kind is everything before the first `:`).

CSS: `.rv-kind` (muted pill, monospace, one colour per kind class) in
`variantgrid/static_files/default_static/css/global.scss` beside `.rv-hgvs`, hand-applied to `global.css`.

CSV: unchanged. The hidden members (`alt__seq`, `svlen`, `variantannotation__variant_class`) are already in
the export, so the CSV reader has what the badge is drawn from. The badge is presentation, not a conversion.

### 3.2 Row expansion for a fusion

`variantopedia/views.py:variant_grid_row_detail` gains, for a gene-level variant, the `GeneFusion`
(`select_related("anchor", "partner")`): the expanded row shows *5′ BCR → 3′ ABL1* or *BCR / ABL1 (direction
not asserted)* and *partner unspecified* where `partner` is null, above the existing "Genes" line
(`variantopedia/templates/variantopedia/variant_grid_row_detail.html`). This is the one place every grid
shares, so the partner/direction detail lives here rather than in five column sets.

### 3.3 Fusions reach the places they currently fall out of

- **Node CSV/VCF export** - `analysis/grids.py:ExportVariantGrid` walks `standard_contigs`, which excludes the
  gene-level contig, so fusions silently vanish from every export. Append the gene-level contig
  (`SequenceRole.VG_GENE_LEVEL_FAKE_CONTIG`) to the loop. The VCF writer emits them as the loader wrote them
  (contig `GENE_LEVEL`, ref `N`, symbolic alt), with a `##contig=<ID=GENE_LEVEL>` header line so the file
  re-imports; that is the representation `snpdb/gene_level_variants.py` chose, and the round trip is the
  test.
- **All Variants page** - add `VariantType.FUSION = "fusion"` to `snpdb/variant_filters.py:VariantType`,
  mapped to `snpdb/models/models_variant.py:Variant.get_gene_level_q`, offered as a "Fusion" button in
  `variantopedia/templates/variantopedia/variants.html`, and let `snpdb/variant_filters.py:get_contigs_q` pass
  the gene-level contig through when it is ticked (the way it already lets a chosen gene's contig through).
  Without the contig, the page's mandatory contig filter hides fusions whatever else is selected.
  Gene-level rows have `svlen = 0`, so "Structural" keeps matching them too - that is correct.
- **Gene page** - nothing to build; a test proves a fusion on either partner appears in
  `genes/grids.py:GeneSymbolVariantsGrid` with the badge.
- **Search by a single partner** - nothing to build; a bare symbol already lands on the gene page, and the
  badge now makes the fusions there recognisable. Pair search stays lookup-only.

---

## 4. Phase 2 — copy number per sample

### 4.1 Decision: read the JSON at query time, no packed column

The TSO 500 plan asked for one decision covering `CN`, `SM` and `SEGID`. It is: **read from the stored
`format` / `info` JSON at query time, keyed by a per-VCF field name**. A packed `samples_copy_number` array
would need a migration on the partitioned `CohortGenotype`, a new importer version and a whole-database
backfill, and would still have to choose between an integer (CN) and a ratio (SM, FC). A grid page is at most a
few hundred rows and sorting is already off above `ANALYSIS_GRID_SORT_MAX_ROWS`, so a JSON path per row costs
nothing that matters. `SEGID` is not surfaced: it is the caller's own gene name (`MYCL1` where we say
`MYCL`), and the gene / overlapping-symbol columns already say which gene a segment hits, resolved properly.

### 4.2 Binding the field

`upload/vcf/vcf_import.py:configure_vcf_from_header` binds `copy_number_field` by name, first match of
`CN`, `SM`, `FC` in the header's FORMAT ids (`get_format_field` in order), else the same names in INFO for a
single-sample VCF (the Pisces shape). A source whose meaning differs overrides it through
`snpdb/models/models_vcf.py:VCFSourceSettings` `sample_field_overrides`, which is why it joins
`OVERRIDABLE_SAMPLE_FIELDS`. It is editable on the VCF page like the other field names. A data migration
sets it on existing VCFs whose `VCFFormat` / `VCFInfo` rows declare one of the three ids (importer v25+
stores those; older VCFs have nothing to read anyway).

### 4.3 Showing it

The per-sample value joins the zygosity cell, which already composes AF, depths, quality marks and filters
(`VariantGridFormat.sampleZygosity`): a `CN 3` / `SM 1.49` chip after the depths, the header description of
that field in the hover. Mechanics, all in `analysis/grids.py:VariantGrid` beside the existing sample columns:

- `snpdb/grid_columns/grid_sample_columns.py:get_available_format_columns` reports `samples_copy_number`
  when any cohort VCF has `copy_number_field`.
- `get_variantgrid_zygosity_annotation_kwargs` annotates, per cohort, `<alias>_packed_samples_copy_number`
  as a `Coalesce` of the JSON paths `format__<i>__<field>__0` (the sample's dict in the per-sample list, then
  the field's one-element array) and `info__<field>` - Django's JSON key transforms handle the list index
  and the key; the sample index in `format` is the VCF column order, which is also the cohort's packed
  index for a VCF's own cohort.
- One hidden `sample_<pk>_samples_copy_number` column per sample carries it to the cell and the CSV; the
  sort menu offers "Copy number" through `_genotype_sort_func` (a `KeyTextTransform` cast to float).
- `_format_sample_value` in `analysis/grid_export.py` writes it into the VCF export's FORMAT under the
  field's own id, so a CNV VCF round-trips its copy number.

---

## 5. Phase 3 — fusion calls, and a structural filter

### 5.1 Fusion calls column (analysis grid)

The loader stores every caller row for a gene pair under `INFO/FUSION_OBS`
(`upload/tasks/import_dragen_tso500_all_fusions_task.py`) - breakpoints, caller, read counts, filters - and
merges several rows onto one `Variant`. A "Fusion calls" column, added by
`analysis/models/nodes/cohort_mixin.py:CohortMixin` next to its Filters column for each cohort whose VCF
declares `FUSION_OBS` in `snpdb/models/models_vcf.py:VCFInfo`, reads `<alias>__info__FUSION_OBS`:

- server renderer decodes the percent-encoded JSON once and returns the text the CSV gets:
  `DRAGEN chr2:29446394→chr2:42492091 (12 reads); …`
- client renderer draws `×3` with that text as the hover, so a merged pair reads as more than one call.

Blank on a row that is not a fusion, which is why the column only appears when a fusion VCF is in the node.
Grids without a cohort (gene page, All Variants) do not carry it; the row expansion (§3.2) is their fusion
detail.

### 5.2 Structural filter on the Effect node

`DamageNode.structural` (§2): *Structural only* ANDs `snpdb/models/models_variant.py:Variant.get_symbolic_q`
into `and_filters` beside the class restriction; *Exclude structural* negates it. Existing nodes default to
*Any*, so nothing changes for them. It answers the question the class filter cannot - "just the SVs/CNVs" -
in the node where variant type already lives (`analysis/forms/forms_nodes.py:DamageNodeForm`, radio beside
the grouped class checkboxes), so counts and downstream nodes see it. A grid-level toggle was considered
and rejected: the node grid deliberately has no filter toolbar (`VariantGrid.filter_builder_toolbar`), because
a filter that changes the rows without changing the node's count misleads.

---

## 6. Tests

Kept if they cover logic we wrote; the framework's JSON path lookups and choice rendering are not ours.

- `analysis/tests/test_representative_variant_column.py` - the export includes a gene-level variant (§3.3)
  and the VCF export writes its contig line; a fusion VCF cohort adds the Fusion calls column and a small
  variant row renders it blank (§5.1).
- `variantopedia/tests/test_all_variants_grid.py` - the Fusion variant type
  selects gene-level rows and passes their contig through the contig filter (§3.3).
- Copy number (§4): binding order CN → SM → FC from a header; the INFO fallback for a single-sample VCF;
  the annotation kwargs produce one alias per cohort and the sort function casts it. Fixture:
  `genes/tests/gene_fusion_test_utils.py:create_gene_fusion` for fusions,
  `snpdb/tests/utils/vcf_testing_utils.py:slowly_create_test_variant` for symbolic alts, and the TSO 500
  CNV file in `upload/test_data/tso500/` through the VCF import tests for the FORMAT shapes.
- `DamageNode.structural` - only / exclude produce the expected Q; *Any* adds nothing (§5.2).
- Gene page shows a fusion from either partner (§3.3).

The JS has no test runner (`package.json` declares none); `_variantKind` is checked by `vg page` on an
analysis holding the loaded TSO 500 files, and by eye.

---

## 7. Order, and what each step touches

| Step | Files | Migration |
|---|---|---|
| §3.1 badge | `variantgrid_formats.js`, `global.scss` + `global.css` | - |
| §3.2 row expansion | `variantopedia/views.py`, `variant_grid_row_detail.html` | - |
| §3.3 export / All Variants | `analysis/grids.py`, `analysis/grid_export.py`, `snpdb/variant_filters.py`, `variants.html` | - |
| §4 copy number | `snpdb/models/models_vcf.py`, `upload/vcf/vcf_import.py`, `grid_sample_columns.py`, `analysis/grids.py`, `variantgrid_formats.js`, `analysis/grid_export.py` | `snpdb`: add `VCF.copy_number_field`; data migration binding existing VCFs from their header rows |
| §5.1 fusion calls | `analysis/models/nodes/cohort_mixin.py`, `variantgrid_formats.js` | - |
| §5.2 structural filter | `damage_node.py`, `forms_nodes.py`, the Effect node editor template | `analysis`: add `DamageNode.structural` |

Phase 1 is independent and small; Phases 2 and 3 are independent of each other. After §4 and §5.2 run
`scripts/vg map` (a model changed). Docs in the same change: one line each in `uicore/CLAUDE.md` (grids:
the kind badge reads members riding along, and where per-sample JSON values are annotated),
`claude/domain.md` (**Variant kind**: badge vocabulary vs VEP `variant_class`), and the Phase 6 entry of
`claude/plans/tso500_overall_plan.md` pointing here.

## 8. Deferred, deliberately

- Gene-level `<AMP>` / `<LOSS>` rows - the badge and the fusion row expansion already handle the alt, and
  nothing produces them yet.
- Breakend coordinates as real `Variant`s, and any fusion equivalence / discordance - #1506 phases 2 and 3.
- A "fusions involving this gene" summary on the gene page. The grid with the badge does this; a summary is a
  follow-up if curators ask for it.
- #1706's grid half (specimen / extraction grids under the patients menu) - stays in the TSO 500 plan.
