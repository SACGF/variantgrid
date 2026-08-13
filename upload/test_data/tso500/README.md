# TSO 500 test data

Illumina TruSight Oncology 500 output — DRAGEN TSO 500 v2.6.2, pipeline 2.6.2.4. One specimen,
sequenced as a DNA arm and an RNA arm.

```
ExampleSample_2600000001/
├── ExampleSample_DNA_2600000001C/
│   ├── ....hard-filtered.vcf         small variants     93 records
│   ├── ....cnv.vcf                   gene-level CNV     25
│   └── ..._DragenExonCNV.vcf         BRCA1/2 exon CNV    2
└── ExampleSample_RNA_2600000001B/
    ├── ..._SpliceVariants.vcf        splice calls       17
    └── ..._AllFusions.csv            fusions            33
```

The VCFs go through the normal VCF import. `AllFusions.csv` has a parser of its own (the format is
DRAGEN TSO 500's, not a standard) which writes the rows as a VCF of gene-level variants. That VCF then
goes through the normal import too - only the bcftools stages are skipped, as a gene-level locus has no
reference base to check against.

Everything is GRCh37 — b37 contigs renamed with a `chr` prefix, including `hs37d5`, `NC_007605`
and the `GL000*` decoys, with `chrM` at 16569 (rCRS, not hg19's 16571).

The names are synthetic but keep the shape a loader has to parse. `2600000001` is a ten-digit lab
accession identifying the specimen, and the trailing `C` and `B` are container suffixes naming the
two nucleic-acid extractions taken from it — so the pair directory is the specimen and the two arm
directories beneath it are its DNA and RNA extractions.

## Properties to know before using this

**Coordinates in `hard-filtered.vcf` are shifted.** This is a tumour-only panel, so its small
variants were the patient's germline. `scripts/vcf_fuzz.py` moved them, re-reading REF from the
reference FASTA at each new position and rebuilding ALT so records stay internally consistent and
keep their variant type. Records describing one event move together, so the MNV group below is
intact. Positions are valid and self-consistent; they are not where the caller found them.

Positions in the other files are real. `cnv.vcf` segments are the panel's fixed segmentation
intervals — identical for every sample — and splice and fusion calls are somatic events.

**Some rows are reconstructions, not caller output:**

- 16 fusions in `AllFusions.csv` and 2 splice calls in `SpliceVariants.vcf` are transcribed from a
  published Local App v2.2 CombinedVariantOutput
  ([AWGL/TSO500_post_processing](https://github.com/AWGL/TSO500_post_processing)) — the canonical
  oncology set (EML4-ALK, KIF5B-RET, CD74-ROS1, ETV6-NTRK3, TMPRSS2-ERG, FGFR3-TACC3 …) plus EGFR
  exons 2-7 (EGFRvIII) and MET exon 14 skipping. Gene pairs, breakpoints and read counts are as
  published, and each Gene A/Gene B split was checked against the gene's locus. Columns the source
  format does not carry (`Score`, contig/alignment fields) are `N/A`; the grafted splice REF bases
  were read from the FASTA.
- Both `DragenExonCNV.vcf` records are constructed. This run had no large rearrangement, and none
  of the published TSO500 outputs contain one either — they all report `BRCA1 NA` / `BRCA2 NA`.
  The records use only the fields the file's own header declares (`END`, `GENE`, `SVTYPE`,
  `GT:FC`) over real hg19 BRCA1 and BRCA2 spans, one `PASS` `<DEL>` and one `Undetermined` to
  exercise that filter. `REF=N` and `SVTYPE=CNV` follow the sibling `cnv.vcf` convention rather
  than a documented example — Illumina publishes no VCF field specification for this file, and
  describes the equivalent data as landing in `_DragenExonCNV.json`. Treat the exact spelling as
  provisional until a real one is seen. Note also that the CombinedVariantOutput reports these as
  `<LOSS>` where the VCF header declares `<DEL>`.

**Run dates and site paths in the caller command lines are neutralised.** Software versions,
vendor `resource_bundle/…` paths and caller arguments are real — a loader may want the pipeline
version out of them.

## Upload metadata for these files

Two facts these files don't reliably carry are supplied at upload instead, as `genome_build` and
`source` (`upload/upload_metadata.py`; API query params, or `import_vcf --genome-build/--source`).

| File | `genome_build` | `source` |
|---|---|---|
| `hard-filtered.vcf` | from header contigs | `DRAGEN TSO500 SmallVariant` |
| `cnv.vcf` | from header contigs | `DRAGEN TSO500 CNV` |
| `_DragenExonCNV.vcf` | **`GRCh37` — required**, the header has no contigs and an unresolvable `##reference` | from header (`LrCalculator 1.0.0.11`) |
| `SpliceVariants.vcf` | from header contigs | from header (`SpliceGirl 1.0.0.614`) |
| `AllFusions.csv` | **`GRCh37` — required** on a multi-build deployment, the file carries no build at all | from its own `# Source =` line (`FusionProcessor 1.0.0.614`) |

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
