# TSO 500: import, filter and report the way Molecular Oncology's SOP does

Written by Claude Opus 5.5 (claude-opus-5-5), 2026-10-06
Status: in progress (§1 + tumour fraction landed #2097; §2, §3 in PR 2)

[sapath#457](https://github.com/SACGF/variantgrid_sapath/issues/457). SA Pathology's Molecular Oncology (MO) analyse
TSO 500 (DRAGEN TSO 500 v2.6.2, GRCh37) by a controlled SOP, summarised in the sapath repo's
*docs/tso500_mo_sop.md*. This plan makes VG reproduce each rule in that SOP:

- the fusion set MO review
- the copy number changes they report, and the copy count they print
- BRCA1/2 Large Rearrangements
- the named splice variants
- MSI, GIS and HRD status
- the HRD tumour fraction comments

Small variant curation is unchanged.

| Phase | Repo | What |
|---|---|---|
| 1 | variantgrid | Data: SampleNode copy-ratio thresholds, EGFRvII splice event, new settings |
| 2 | variantgrid | AllFusions FILTER: kept, rescued, not kept |
| 3 | variantgrid | Tumour fraction as a fraction, GIS call, MSI-High tumour fraction rule, estimated copies on the report |
| 4 | variantgrid_sapath | Settings, HRD status / BRCA status / tumour fraction comments in the report and JSON |
| 5 | NGS-pipelines | Send *_DragenExonCNV.vcf* |
| 6 | (data) | TSO 500 analysis template nodes |
| 7 | variantgrid | `upload/test_data/tso500/README.md` |

## 1. Data

### `SampleNode` - copy ratio thresholds

`analysis/models/nodes/sources/sample_node.py:SampleNode`. A copy number call has a ratio against the normal when its VCF's
copy number field is one (`VCFConstant.COPY_NUMBER_FIELD_IS_RATIO`: `SM`, `FC`). MO report a gene amplification at a
fold change of 2.5 or more, and investigate a BRCA1/2 Large Rearrangement at a fold change below 0.5.

```python
class SampleNode(SampleMixin, GeneCoverageMixin, AnalysisNode):
    ...
    min_copy_gain_ratio = models.FloatField(null=True, blank=True)  # a gain (ratio > 1) passes at or above this
    max_copy_loss_ratio = models.FloatField(null=True, blank=True)  # a loss (ratio < 1) passes at or below this
```

- A row passes when either:
  - it has no copy ratio: every small variant, and every record from a VCF whose copy number field is absolute
    (`CN`) or absent; or
  - its ratio is > 1 and ≥ `min_copy_gain_ratio`; or
  - its ratio is < 1 and ≤ `max_copy_loss_ratio`.
- A null threshold lets every gain or loss through.
- The q is built on `snpdb/grid_columns/grid_sample_columns.py:get_copy_number_annotation`, cast to float, and appears
  in the node's description ("gain FC>=2.5", "loss FC<=0.5").
- The two fields are node-only: not in `THRESHOLD_FIELDS` and not on `SampleNodeSampleFilter`. At a group level they
  apply to each sample whose VCF has a ratio; an override would only matter for two ratio callers needing different
  thresholds in one node, and the template gives each copy number source its own node.
- Node form: two fields under the existing thresholds. Shown only when the sample's VCF has a ratio field.

### `SpliceEvent` - EGFRvII

A new data migration (`genes/migrations/0101_seed_egfr_vii_splice_event.py`) adds the second EGFR junction MO
report by name, in the same shape as `genes/migrations/0093_seed_splice_events.py`:

```python
# EGFRvII: exon 13 spliced to exon 16, ie exons 14-15 skipped (NM_005228.5), aa 521-603 in frame
{"gene_symbol": "EGFR", "label": "v_ii", "display": "EGFRvII splice variant", "contig": "7",
 "GRCh37": (55229324, 55238867), "GRCh38": (55161631, 55171174)},
```

- Coordinates come from NM_005228.5's exon table in cdot, in each build's own alignment. GRCh37 has no MANE; .5 is
  the transcript MO report, and EGFR's MANE Select in GRCh38. They are 0-based exactly as the existing rows are, and
  are checked the same way.
- Matching is on the genomic breakpoints DRAGEN writes, so the analysis's canonical transcript collection plays no
  part here.
- The label is canonicalised by `genes/gene_splice.py:canonical_splice_label`, as the migration 0095 rows are.

### Settings - `variantgrid/settings/components/tso500_settings.py`, beside `TSO500_MSI_CALL_BANDS`

```python
TSO500_GIS_CALL_BANDS = None                # over 'Genomic Instability Score', eg [(42, "POSITIVE"), (0, "NEGATIVE")]
TSO500_GIS_MIN_TUMOR_FRACTION = None        # below this caller tumour fraction a GIS under the top band gets no call, eg 0.23
TSO500_MSI_HIGH_MIN_TUMOR_FRACTION = None   # below this caller tumour fraction the top MSI band gets no call, eg 0.20
TSO500_FUSION_RESCUE_MIN_SCORE = None       # AllFusions rows not kept by DRAGEN, Score strictly above this, eg 0.5
TSO500_FUSION_RESCUE_FILTER = None          # ... and Filter exactly this, eg "FAIL;LOW_MAPQ"
```

Each one is lab policy: unset means no call / no rescue, as with the MSI and TMB bands.

## 2. AllFusions FILTER

`upload/tasks/import_dragen_tso500_all_fusions_task.py:_write_gene_level_vcf` writes one record per resolved gene pair
from one or more AllFusions rows (observations). It now passes a `vcf_filter` decided from those observations:

| Any observation with | FILTER |
|---|---|
| `KeepFusion` True | `PASS` - the CVO `[Fusions]` candidates |
| `Score` > `TSO500_FUSION_RESCUE_MIN_SCORE` and `Filter` == `TSO500_FUSION_RESCUE_FILTER` | `LowMapQRescue` |
| neither | `NotKept` |

- The first matching row of the table wins.
- Both filters are declared in the written header:
  - `##FILTER=<ID=LowMapQRescue,Description="Not kept by DRAGEN; Score and Filter match the lab's rescue rule">`
  - `##FILTER=<ID=NotKept,Description="KeepFusion is false">`
- With the rescue settings unset, `LowMapQRescue` is never written.
- DRAGEN's own `Filter` string stays per observation inside `FUSION_OBS`, as it is now.
- `Score` reads `N/A` as None, which never exceeds the threshold.
- Tests (`upload/tests/`, beside the existing AllFusions import tests):
  - one pair per row of the table;
  - a pair with a kept and a not-kept observation is `PASS`;
  - `Score` exactly 0.5 is not rescued.

## 3. Calls and copies

### Tumour fraction is a fraction

DRAGEN 2.6 writes the CVO's `[GIS]` `Tumor Fraction` as a percent (`55`), where 2.1 wrote a fraction (`0.62`).
`DragenTSO500CombinedVariantOutput.tumor_fraction` is a fraction, so
`upload/tso500/dragen_combined_variant_output_records.py:combined_variant_output_values` divides a value above 1 by
100. A tumour fraction of 1% written as `1` is indistinguishable from 1.0 and stays 1.0; DRAGEN does not call that
low.

Test: `55` → 0.55, `0.62` → 0.62.

### GIS call

`seqauto/models/models_seqauto.py:DragenTSO500CombinedVariantOutput` gets `gis_call`, shaped like `msi_call`, and
`measure_call("gis")` returns it:

- Returns None when `TSO500_GIS_CALL_BANDS` is unset or `genomic_instability_score` is None.
- Otherwise the call is `band_call(genomic_instability_score, bands)`.
- If the call is not the top band, `TSO500_GIS_MIN_TUMOR_FRACTION` is set, and `tumor_fraction` is below it, the
  result is `MeasureCall(None, threshold, source)`. A low score is not evaluable at low purity; a high one stands.
- The threshold text names both settings, as `msi_call`'s does.

### MSI-High tumour fraction

In `msi_call`, a top-band call with `TSO500_MSI_HIGH_MIN_TUMOR_FRACTION` set and `tumor_fraction` below it becomes
`MeasureCall(None, threshold, source)`. A None `tumor_fraction` leaves the call as it is: a run without the HRD arm has
no tumour fraction to hold it to.

### Estimated copies

`classification/report/case_report_context.py:ReportVariant` is built from the classification's evidence:
- `copy_number` is filled when the caller wrote an absolute count (`CN`);
- `fold_change` is filled when it wrote a ratio (`SM` / `FC`, `classification/autopopulate_evidence_keys/evidence_from_sample_and_patient.py:get_copy_number_evidence`).

Where a `COPY_NUMBER` kind row has no `copy_number` but has a `fold_change`, `copy_number` is `round(2 × fold_change)`.
The ratio is against a diploid normal, so 2.5 is 5 copies. Losses keep `copy_number` None. The sapath template and
JSON already print `copy_number` ("AR amplification 5 copies"). Test: 2.5 → 5, a loss → None, an explicit count wins.

Tests for this phase:
- `gis_call`:
  - 45 at TF 0.1 → POSITIVE
  - 30 at TF 0.5 → NEGATIVE
  - 30 at TF 0.1 → None with threshold
  - 30 at TF None → NEGATIVE
- `msi_call`:
  - 35% at TF 0.1 → None
  - 35% at TF None → MSI-High
  - 15% at TF 0.1 → MSI-Low

## 4. SA Path report (variantgrid_sapath)

### Settings

In the sapath repo's (`../variantgrid_sapath`) *variantgrid/settings/env/vg4test.py*, beside the MSI/TMB bands. It is
the one active VG4 SA Path environment, and the production settings will be copied from it once testing is done:

```python
TSO500_GIS_CALL_BANDS = [(42, "POSITIVE"), (0, "NEGATIVE")]
TSO500_GIS_MIN_TUMOR_FRACTION = 0.23
TSO500_MSI_HIGH_MIN_TUMOR_FRACTION = 0.20
TSO500_FUSION_RESCUE_MIN_SCORE = 0.5
TSO500_FUSION_RESCUE_FILTER = "FAIL;LOW_MAPQ"
```

The existing MSI bands (High ≥ 30, Low ≥ 10, MSS) and TMB bands (High ≥ 10) stay; they match the report's methods
paragraph.

### BRCA variant status

In `sapath/tso500_case_template.py`, `brca_status` becomes `{"type": "choice", "options": ["POSITIVE", "NEGATIVE"]}`.
POSITIVE means a P/LP variant in BRCA1 or BRCA2; the scientist sets it.

### HRD status and comment

A new *sapath/templatetags/sapath_report_tags.py* holds the one rule, used by the template and by
`sapath/tso500_report.py:_hrd`. `gis` is `measures.gis`, and `tumour_fraction` is `measures.tumour_fraction_sequencing`
(the caller's, as a fraction).

| GIS call | BRCA status | HRD status | Comment |
|---|---|---|---|
| POSITIVE | any | POSITIVE | TF ≥ 0.23: "The tumour purity of this specimen met the minimum 23% requirement (bioinformatic assessment)." TF < 0.23: "The low tumour purity of this specimen may result in less accuracy of the Genomic Instability Score (GIS)." |
| NEGATIVE | POSITIVE | POSITIVE | met-the-minimum sentence |
| NEGATIVE | NEGATIVE / blank | NEGATIVE | met-the-minimum sentence |
| None (not evaluable) | POSITIVE | POSITIVE | "The tumour purity of this specimen is insufficient for accurate determination of Genomic Instability Score (GIS)." |
| None (not evaluable) | NEGATIVE / blank | INCONCLUSIVE | the insufficient sentence + "Testing of different tissue block with higher tumour percentage may be considered." |
| no GIS value | any | - (not an HRD case) | none |

- `hrd_status(measures, case_values)` returns the status column.
- `hrd_comment(measures, case_values)` returns the comment column; the 23% in the comment comes from
  `TSO500_GIS_MIN_TUMOR_FRACTION`.

Changes in `sapath/data/reports/tso500_case_report.html`:
- The `HRD Status` line prints `hrd_status`.
- The `Genomic Instability Score` line keeps `(non-HRD range: 0-41)`. When the GIS call is None (not evaluable), it
  prints `Not evaluated, insufficient tumour purity.` instead of the score.
- An INCONCLUSIVE case adds `No actionable variants detected in BRCA1/2 genes.`
- `BRCA Variant Status` prints `case_values.brca_status`.
- The fixed "Tumour percentage met minimum 23% requirement for HRD assessment" sentence goes; the Comment section
  prints `hrd_comment`.

The template is loaded into the database by a migration (`sapath/migrations/0014_tso500_report_template.py`), so a new
sapath migration reloads it, and the case fields, from the file. `_hrd` sends `{"HRD": {"Category": hrd_status,
"Score": ...}}`. Tests in `sapath/tests/test_tso500_report_template.py` cover each row of the table, using the
`sapath/tests/tso500_case.py` context.

## 5. NGS-pipelines

The TSO 500 upload sends *<sample>_DragenExonCNV.vcf* from the DNA arm, with genome build GRCh37 as upload metadata.
VG imports it as coordinate SVs with `FC` as the copy-ratio field (already handled). This goes on the davmlaw fork's
TSO 500 upload branch, beside the other files it sends.

## 6. Analysis template

The SA Path TSO 500 analysis template (`TSO500_combo`) gets these nodes, edited in the UI on each deployment, and the
node table in the sapath repo's *docs/tso500_mo_sop.md* gains them:

| Node | Source | Settings |
|---|---|---|
| Gene amplifications | cnv.vcf sample, excluding DNAJB1, FANCF, FOXL2, HIST1H3A, HIST1H3C-J, HIST2H3D, TERC, TERT | `min_copy_gain_ratio` 2.5. Losses pass, for verification before any report |
| Large Rearrangements | *_DragenExonCNV.vcf* sample, gene list BRCA1, BRCA2 | `max_copy_loss_ratio` 0.5 |
| Fusions | AllFusions sample | FILTER `PASS` or `LowMapQRescue` |
| Fusions - all | AllFusions sample | none: the *_AllFusions_filtered.csv* review set, sorted by score |
| Splice variants | SpliceVariants.vcf sample | FILTER `PASS` |

## 7. Test data README

Changes to `upload/test_data/tso500/README.md`:
- The CVO copy number rule is every *cnv.vcf* record with `FILTER=PASS` and ALT `<DUP>` or `<DEL>`.
- DRAGEN 2.6 section names: `[Copy Number Variants]`, `[Large Rearrangements]`, `[Loss of Heterozygosity]`.
- Sample IDs: `<index>_MO_TSO_<DNA|RNA>_<C-number>_<initials>_<accession><container>`.
- Pair ID: `<index>_<C-number>_<initials>_<accession>`. The examples are invented ones of that shape.

## Waiting on MO

These change only numbers or wording above, never the shape:
- The 0.23 HRD tumour fraction, which matches the report's existing sentence.
- RNA on-target reads: the report's methods paragraph says > 9M, the SOP and the DRAGEN 2.6.2 file guideline say
  2.5M.
- Exon numbers in fusion descriptions. They follow MO's transcript choices, ie the analysis's canonical transcript
  collection (`genes/models/models_gene_coverage.py:CanonicalTranscriptCollection`), not MANE: TSO 500 is GRCh37,
  and MO pick per gene.
- Whether a BRCA1/2 Large Rearrangement sets BRCA status POSITIVE.
- The CNVkit output (a later file type).
