# TSO 500 test data

Illumina TruSight Oncology 500 output — DRAGEN TSO 500 v2.6.2, pipeline 2.6.2.4. One specimen,
sequenced as a DNA arm and an RNA arm.

```
MetricsOutput_orig.tsv                run-level: one column per pair on the run, both arms
ExampleSample_2600000001/
├── ..._CombinedVariantOutput.tsv     the pair's reportable calls + TMB/MSI/GIS
├── ..._MetricsOutput.tsv             run QC, analysis status and library QC per extraction
├── ExampleSample_DNA_2600000001C/
│   ├── ....hard-filtered.vcf         small variants     93 records
│   ├── ....cnv.vcf                   gene-level CNV     25
│   └── ..._DragenExonCNV.vcf         BRCA1/2 exon CNV    2
└── ExampleSample_RNA_2600000001B/
    ├── ..._SpliceVariants.vcf        splice calls       18
    └── ..._AllFusions.csv            fusions            33
```

The VCFs go through the normal VCF import, the CNV and splice ones rewritten as gene-level records
first. `AllFusions.csv` has a parser of its own (the format is DRAGEN TSO 500's, not a standard) which
writes the rows as a VCF of gene-level variants. That VCF then
goes through the normal import too - only the bcftools stages are skipped, as a gene-level locus has no
reference base to check against.

Everything is GRCh37 — b37 contigs renamed with a `chr` prefix, including `hs37d5`, `NC_007605`
and the `GL000*` decoys, with `chrM` at 16569 (rCRS, not hg19's 16571).

The names are synthetic but keep the shape a loader has to parse. `2600000001` is a ten-digit lab
accession identifying the specimen, and the trailing `C` and `B` are container suffixes naming the
two nucleic-acid extractions taken from it — so the two arm directories are the specimen's DNA and
RNA extractions. The pair is **the patient**, not the specimen: the pair ID carries the Omico C-number
(`Patient.patient_code`), written either as the whole pair sample name (`5_C0000001_FCUP_2600000001`) or
as the C-number on its own (`C0000001`) - the lab's pipeline chooses inconsistently, so both forms turn up
in the same feed and `settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX` reads either. The sample IDs are
`SA-<C-number>-<accession><container>-<D|R>`, eg `SA-C23755-2535115161C-D`. A patient re-analysed later
comes back under the same C-number with a new accession, so one pair ID spans specimens. The test files keep the older synthetic sample names
because a CVO's DNA/RNA Sample IDs have to equal the VCF sample names they link to; the pair ID
(`C0000001`) has the real shape.

## Properties to know before using this

**Coordinates in `hard-filtered.vcf` are shifted.** This is a tumour-only panel, so its small
variants were the patient's germline. `scripts/vcf_fuzz.py` moved them, re-reading REF from the
reference FASTA at each new position and rebuilding ALT so records stay internally consistent and
keep their variant type. Records describing one event move together, so the MNV group below is
intact. Positions are valid and self-consistent; they are not where the caller found them.

Positions in the other files are real. `cnv.vcf` segments are the panel's fixed segmentation
intervals — identical for every sample — and splice and fusion calls are somatic events.

**`CombinedVariantOutput.tsv` is assembled, not a caller's.** It is the pair-level summary DRAGEN writes
next to the two arm directories: `[Analysis Details]` names the patient (the code in the second field of
`Pair ID`) and the two extractions (`DNA Sample ID`, `RNA Sample ID`), then one `[Section]` per call type, each a header row then one row per
call, every line padded with tabs to the widest section (11 columns). The layout is a real 2.1.1 RUO
file's; the rows are this sample's. The real file's identifiers, dates and per-sample metrics were
replaced, so:

- `[Splice Variants]` has three rows, and they are the three `PASS` calls of `SpliceVariants.vcf`
  written the way the CVO writes them (`Breakpoint 1` = `POS`, `Breakpoint 2` = `END`,
  `Splice Supporting Reads` = `ALTDEDUP`, `Reference Reads Transcript` = `REFDEDUP`, `Affected Exon` =
  the exons skipped). The `AR` row and its VCF record are a real pair's, with the duplicate-inclusive
  counts changed (SACGF/variantgrid_sapath#457), and are what the other two are shaped on:
  `chrX:66905968` is the last base of AR exon 3 in GRCh37 (NM_000044.6) and `chrX:66914514` is in
  intron 3 at cryptic exon 3 — exon 3 spliced to CE3 is AR-V7, and `Affected Exon` is blank because
  CE3 is not an annotated exon. SpliceGirl's `REF` is not the base at `POS`: in every real record it
  is the base at `POS+2`, the second base of the intron's donor dinucleotide - `T` for a GT intron,
  `C` for the GC-AG intron at `chr2:42485683` (EML4) - and the three `PASS` records follow that.
  `bcftools norm --check-ref=s` does not touch a symbolic `<DEL>`, so a plain VCF import makes its
  Locus with the caller's base - one reason the file is imported as gene-level splice junctions
  instead (`upload/tasks/import_splicegirl_vcf_task.py`, #1903), and the section here is not a
  variant source.
  Illumina's rule for the section is passing calls on EGFR, MET and
  AR only - a gene filter over `FILTER=PASS`, not a junction whitelist - so a `PASS` call in any other
  gene would be in the VCF and not here; this file has none, so it cannot tell the two rules apart
  (source cited in `upload/tasks/import_dragen_tso500_combined_variant_output_task.py`).
- `[Fusions]` is every `KeepFusion = True` row of `AllFusions.csv`: gene pair joined with `-` where
  the caller knew the direction, supporting reads as `upload/tso500/dragen_all_fusions_parser.py`
  sums them, and `Ref A Dedup` / `Ref B Dedup` as the two reference counts.
- `[Gene Amplifications]` is the `cnv.vcf` `<DUP>` segments with `SM` ≥ 1.5, to three decimals. The
  cutoff is a guess at the caller's; Illumina does not publish it.
- `[Small Variants]` is one row, the first `PASS` record of the fuzzed `hard-filtered.vcf` with its
  HGVS recomputed at the shifted position. The section is the reportable subset of that VCF and is
  **not** what small variants are loaded from - see the upload table below.
- `[Exon-Level CNVs]` is `BRCA1 NA` / `BRCA2 NA` as every published file has it, so it disagrees
  with the constructed `_DragenExonCNV.vcf` records on purpose.
- TMB, MSI and GIS values are plausible and made up. `Coding Region Size in Megabases` is the
  panel constant.
- Module 2.1.1 writes `[Exon-Level CNVs]`; 2.6 documents that section as `Large Rearrangements` and
  adds `Gene-level Loss of Heterozygosity`. A loader keys on section names, and tolerates a section
  it does not know.

**Some rows are reconstructions, not caller output:**

- 16 fusions in `AllFusions.csv` and 2 splice calls in `SpliceVariants.vcf` are transcribed from a
  published Local App v2.2 CombinedVariantOutput
  ([AWGL/TSO500_post_processing](https://github.com/AWGL/TSO500_post_processing)) — the canonical
  oncology set (EML4-ALK, KIF5B-RET, CD74-ROS1, ETV6-NTRK3, TMPRSS2-ERG, FGFR3-TACC3 …) plus EGFR
  exons 2-7 (EGFRvIII) and MET exon 14 skipping. Gene pairs, breakpoints and read counts are as
  published, and each Gene A/Gene B split was checked against the gene's locus. Columns the source
  format does not carry (`Score`, contig/alignment fields) are `N/A`.
- Both `DragenExonCNV.vcf` records are constructed. This run had no large rearrangement, and none
  of the published TSO500 outputs contain one either — they all report `BRCA1 NA` / `BRCA2 NA`.
  The records use only the fields the file's own header declares (`END`, `GENE`, `SVTYPE`,
  `GT:FC`) over real hg19 BRCA1 and BRCA2 spans, one `PASS` `<DEL>` and one `Undetermined` to
  exercise that filter. `REF=N` and `SVTYPE=CNV` follow the sibling `cnv.vcf` convention rather
  than a documented example — Illumina publishes no VCF field specification for this file, and
  describes the equivalent data as landing in `_DragenExonCNV.json`. Treat the exact spelling as
  provisional until a real one is seen. Note also that the CombinedVariantOutput reports these as
  `<LOSS>` where the VCF header declares `<DEL>`.

**`MetricsOutput.tsv` is a lab file with its identifiers and values replaced.** DRAGEN writes one per run
beside the CombinedVariantOutput; the lab's pipeline copies it under the pair's name. The layout is a real
2.1.1 file's: `[Header]`, `[Run QC Metrics]` (one `Value` column, run-level), `[Analysis Status]`, then a
`[... QC Metrics]` section per category with `Metric (UOM)`, `LSL Guideline`, `USL Guideline` and one column
per sample, and two `[... Expanded Metrics]` sections with no guidelines. Every sample is a column in every
section - a DNA library has `NA` down the RNA sections and vice versa, so the RNA extraction's rows in
`[DNA Library QC Metrics]` are `NA` here and the DNA extraction's in `[RNA Library QC Metrics]` are too.
A metric passes when `LSL <= value <= USL`, an `NA` guideline being no bound; a sample whose column is all
`NA` in a section simply has no library of that kind. Section names and metric lists differ between
versions (2.6 adds `PCT_CHIMERIC_READS` to the small-variant section and `EXCESSIVE_TF` to GIS, and drops
`PCT_PF_UQ_READS`), so a loader keys on section and metric names and tolerates ones it does not know.
Note `MEDIAN_INSERT_SIZE` appears in both the DNA small-variant and the RNA sections, so a metric is only
unique within its section. The lab's file pads every line with tabs to the widest row; a 2.6.2 file
seen in the wild does not, so the padding is not to be relied on either way.

The identifiers, dates and every value were replaced. Its two columns are the CVO's `DNA Sample ID` and
`RNA Sample ID` - a shape DRAGEN 2.1.1 wrote and 2.6.2 does not, kept here as a parser case, since a
section's sample columns are read the same way whichever the file names. The values are plausible and made
up, kept within the guidelines so the file agrees with the CVO (a completed run with calls in every
category); `USABLE_MSI_SITES` equals the CVO's `Usable MSI Sites`. The run QC metrics are populated, as the
lab starts analysis from BCLs, not FASTQs (see `[Notes]`).

**`MetricsOutput_orig.tsv` is the run-level 2.6.2 file, and is what gets imported.** The lab's run wrapper
rewrites every `MetricsOutput.tsv` in place - inserting a sex metrics block, padding `[Analysis Status]`
rows and stripping the research-use text - before the results reach us; the one untouched DRAGEN copy it
keeps is `Results/MetricsOutput_orig.tsv`, so that is what the pipeline sends, once per sequencing run. A
column here is a **pair**, named by the CVO's `Pair ID` and carrying *both* of its arms - the DNA sections
and the RNA sections have values in the same column - so a column names a specimen (its trailing ten-digit
accession) rather than an extraction. The four columns are `5_C0000001_FCUP_2600000001` (the `ExampleSample`
pair, every metric within guideline and `USABLE_MSI_SITES` agreeing with its CVO),
`7_C0000002_ABCD_2600000002` (`MEDIAN_EXON_COVERAGE` 92 against a guideline of 150, so its Small Variants /
TMB category fails, and `COMPLETED_ALL_STEPS` `FALSE` with a `FAILED_STEPS` entry), `9_0PRI_2600000003` (a
pair with no patient C-number, as controls and research samples are named, all within guideline),
`C0000004` (the bare form the lab also writes, so the accession has to come off the run's sample sheet -
`1_TSO_DNAHRD_C0000004_2600000004C_B4`) and `11_NTC`, which nothing on a run is named for, so its rows
are kept with the claim parked. 2.6.2 adds `PCT_CHIMERIC_READS` to the
small-variant section and `EXCESSIVE_TF` to GIS, and its `TOTAL_ON_TARGET_READS` guideline is 2,500,000
where the lab's methods paragraph says 9M - which is what `settings.TSO500_LIBRARY_QC_GUIDELINES` is for.
The file names the run nowhere, and a pair column alone does not identify a pair across runs, so the run
comes in as upload metadata and is part of the `LibraryQC` key. The file names no sample either: each
row's arm is linked to its `SequencingSample` through the sheet's `Pair_ID` and `Sample_Type` columns,
which the pipeline posts per sample.

**Run dates and site paths in the caller command lines are neutralised.** Software versions,
vendor `resource_bundle/…` paths and caller arguments are real — a loader may want the pipeline
version out of them.

## Upload metadata for these files

Two facts these files don't reliably carry are supplied at upload instead, as `genome_build` and
`source` (`upload/upload_metadata.py`; API query params, or `import_vcf --genome-build/--source`).

| File | `genome_build` | `source` | `sequencing_run` |
|---|---|---|---|
| `hard-filtered.vcf` | from header contigs | `DRAGEN TSO500 SmallVariant` | — |
| `cnv.vcf` | from header contigs | `DRAGEN TSO500 CNV` | — |
| `_DragenExonCNV.vcf` | **`GRCh37` — required**, the header has no contigs and an unresolvable `##reference` | from header (`LrCalculator 1.0.0.11`) | — |
| `SpliceVariants.vcf` | from header contigs | from header (`SpliceGirl 1.0.0.614`) | — |
| `AllFusions.csv` | **`GRCh37` — required** on a multi-build deployment, the file carries no build at all | from its own `# Source =` line (`FusionProcessor 1.0.0.614`) | — |
| `CombinedVariantOutput.tsv` | **`GRCh37` — required**, no build in the file | from `Module Version` (`DRAGEN TSO500 CombinedVariantOutput 2.1.1`) | — |
| `MetricsOutput_orig.tsv` | **none accepted** - the file has no coordinates | from `Workflow Version` (`2.1.1.4`) | **required** - the only key it takes |

Send a build's **own name** (`GRCh37`), not an alias (`hg19`). These files are GRCh37 with a `chr`
prefix and `chrM` at 16569, and their `##reference` says `hg19_decoy` — which is exactly the confusion
the declared build exists to settle. Aliases do resolve, but `hg19` is both GRCh37's alias and a build
in its own right, so it only reads unambiguously while that build stays disabled.

A client-supplied `source` becomes part of VG's configuration contract, since `VCFSourceSettings.
source_regex` has to match something the client invented — so these strings want to stay stable. The
two DRAGEN ones are reserved rather than used: nothing is keyed on them yet, because the SpliceGirl
field mapping comes off the header and the copy-neutral skip is a general rule.

## Cases this data covers

- `hard-filtered.vcf` — one MNV represented three ways at overlapping positions: an SNV and a 7bp
  delins at the same POS, plus the decomposed second SNV 6bp later, tied together by `MNVTAG` and
  a shared `PS` phase set. Also `PS`-phased rows, `multiallelic`, `excluded_regions`, `hotspot`,
  all three `GermlineStatus` values, and a spread of compound `FILTER` strings.
- `cnv.vcf` — copy-neutral records with `ALT=.` alongside `<DUP>`/`<DEL>`. `FORMAT` is
  `GT:SM:BC:PE` with no AD/DP/AF, `SM` being the linear copy ratio. `SEGID=MYCL1` uses an older
  symbol than the rest of the pipeline (`MYCL`).
- `SpliceVariants.vcf` — no `GT`; `AD`/`DP` carry splice-specific meanings (see the header) and
  the sample column is literally `SAMPLE`. `chr2:47637511` appears twice with different `END`.
- `CombinedVariantOutput.tsv` — a splice call with no `Affected Exon` (AR-V7, a cryptic exon), two
  with one (`14`) and a range (`2-7`); a gene pair whose 5' side is written with a slash
  (`PPARG/AC016683.6-PAX8`) and one whose 3' side is a semicolon list (`CD74-ROS1;GOPC`); `NA`
  sections (`[Sequencing Run Details]`, both exon-level CNV genes); tab padding on every line.
- `MetricsOutput.tsv` — the `[Analysis Status]` header row starts with an empty cell; `[Run QC Metrics]`
  has a `Value` column where the others have samples; `NA` as a guideline (no bound), as a value (arm not
  sequenced) and as a `(UOM)`; a metric name repeated across sections; tab padding on every line.
- `MetricsOutput_orig.tsv` — the same, run-level and 2.6.2: a column per pair carrying both arms, a pair
  with no C-number, a pair ID in the bare C-number form, a column naming nothing at all, a metric outside
  its guideline, a pair whose run did not complete, and a banner line carrying the module version.
- `AllFusions.csv` — multi-gene partners (`RP11-458D21.5;NOTCH2NL`, `ROS1;GOPC`), a gene pair
  written with a slash (`PPARG/AC016683.6`), `SEPT14` (renamed `SEPTIN14` by HGNC, and
  date-mangled by spreadsheets), two callers, long semicolon-joined filter strings.

## Missing

Per-gene absolute copy number and minor copy number — the latter being what gene-level LOH is
derived from — are not in any file here. Illumina puts both in
`Logs_Intermediates/Gis/<sample>/<sample>.abcn_annotated.vcf`, with a companion
`<sample>.abcn_genes.tsv`, written by the PhenoHRD step. Being a VCF it is in scope for ingestion
and belongs here once obtained. The published documentation does not specify its INFO or FORMAT
fields, so it cannot usefully be mocked up in the meantime.

Sources: [DRAGEN TSO 500 v2.6 Combined Variant
Output](https://help.tso500software.illumina.com/dragen-tso-500-guides/dragen-tso-500-v2.6/analysis-output/combined-variant-output),
[DNA Analysis
Methods](https://help.tso500software.illumina.com/dragen-tso-500-guides/dragen-tso-500-v2.6/overview-1/dna-analysis-methods),
[BRCA Within Gene Large Rearrangement](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/brca-lr.htm)
