# TSO 500 ingestion — what loads today, and what to change

Scope: getting the five DRAGEN TSO 500 output files in `upload/test_data/tso500/` through the existing
upload pipeline. Sits under the data plan in [SACGF/variantgrid_sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431);
the file set itself is [#1506](https://github.com/SACGF/variantgrid/issues/1506).

Findings below come from running each file through the real pipe stages
(`write_cleaned_vcf_header` → `vcf_clean_and_filter` → `bcftools norm` → `vcf_clean_alts`) against a
GRCh37 dev database, then feeding the output to `BulkGenotypeVCFProcessor.process_entry`.

---

## Current state

| File | Preprocess | Insert | Verdict |
|---|---|---|---|
| `hard-filtered.vcf` | 93 → 93 | 93 records | Loads |
| `cnv.vcf` | 25 → 25 | 25 records | Loads |
| `SpliceVariants.vcf` | 17 → 17 | 17 records | Loads, sample fields need mapping |
| `_DragenExonCNV.vcf` | dies at first stage | — | Needs a fix |
| `AllFusions.csv` | n/a | n/a | Needs a parser |

Genome build resolves to GRCh37 from the header contigs for the three that carry them: `chrM=16569`
matches GRCh37 (NC_012920) and separates it from hg19 (16571), leaving exactly one build. The
`GL000*` / `hs37d5` / `NC_007605` decoy contigs are stripped from the cleaned header by
`write_cleaned_vcf_header(standard_contigs_only=True)`, and no records sit on them — the skipped
contig, record and filter stats all came back empty for `hard-filtered.vcf` and `cnv.vcf`.

`<DUP>` and `<DEL>` pass `vcf_clean_alts` because both are in `settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES`,
and every symbolic record carries `SVLEN` or `END` so `vcf_get_ref_alt_svlen_and_modification` can size
it. `REF=N` is replaced from the FASTA by `bcftools norm --check-ref=s`.

---

## 1. `_DragenExonCNV.vcf` — read the declared filters with cyvcf2

`upload/management/commands/vcf_clean_and_filter.py:170` `_get_defined_vcf_filters` builds the set of
header-declared FILTER IDs with PyVCF's `Reader`. Illumina's LrCalculator writes:

```
##FILTER=<ID=Undetermined,Number=1,Type=String,Description="Both segments for this gene ...">
```

PyVCF's FILTER regex accepts `ID` and `Description` only, so the extra keys raise
`SyntaxError: One of the FILTER lines is malformed` and the whole import dies at the first pipe stage.
bcftools and cyvcf2 both read the line without complaint, and with it patched out both records go
through preprocessing and insertion.

**Change:** source the declared filter IDs from cyvcf2 — `cyvcf2_header_types(reader)["FILTER"]`, the
same helper `configure_vcf_from_header` already uses in `upload/vcf/vcf_import.py:240`. This is the
last PyVCF use in the pipe, so the `from vcf import Reader` import goes with it.

**Also needed:** this file has no `##contig` lines and `##reference=file:resource_bundle/hg19_decoy/`,
so `vcf_detect_genome_build_from_header` finds nothing to match. That path is already handled —
`upload/tasks/vcf/genotype_vcf_tasks.py:63` sets `ImportStatus.REQUIRES_USER_INPUT` and terminates the
pipeline early — so each of these needs its build set by hand until the loader supplies one. Worth
deciding whether the TSO 500 loader passes the build in explicitly alongside the file.

Two things that turned out to be fine and need no work: the file is contig-unsorted (chr17 before
chr13), which `VCFSortChecker` allows because it only requires each contig's records be grouped; and
`bcftools norm` shifts the `<DEL>` back one base (POS and END both −1) because it assumes a padding
base. The sibling `cnv.vcf` records do carry that padding base, so the shift is an artefact of the
constructed test record rather than something real output will hit — confirm against a real
`_DragenExonCNV.vcf` when one arrives.

## 2. `AllFusions.csv` — needs a parser, and needs detection to survive it

Today an upload of this file errors during file-type detection, before any parser is chosen.
`get_import_task_factory_from_extension` asks all five `csv` factories for their processing ability,
and `PatientRecordsImportTaskFactory.get_processing_ability`
(`upload/import_task_factories/import_task_factories.py:127`) calls `pd.read_csv` bare. The 31-line
`#`-prefixed comment preamble gives `pandas.errors.ParserError: Expected 1 fields in line 28, saw 2`.
If detection got past that, `GeneListImportTaskFactory` returns ability 1 and would win, importing the
fusions as a gene list.

**Changes:**
- Make `get_processing_ability` return 0 for input it cannot read, so one unreadable CSV shape stops
  affecting detection for every other CSV type.
- Add a fusion import task factory that claims this file — recognisable from the `# Source =
  FusionProcessor` first line and the `Caller,Gene A,Gene B,...` header row — with a processing
  ability above `GeneListImportTaskFactory`'s 1.

The format the parser has to handle, all present in the test file: a `#` comment preamble of column
descriptions; `N/A` as the empty value in every column; multi-gene partners semicolon-joined
(`RP11-458D21.5;NOTCH2NL`, `ROS1;GOPC`); a gene pair written with a slash (`PPARG/AC016683.6`);
`SEPT14`, which HGNC renamed to `SEPTIN14` and spreadsheets turn into a date; breakpoints as
`chr8:128806980`; long semicolon-joined filter strings (`FAIL;LOW_SCORE;MIN_SUPPORT;…`); two callers
(`DRAGEN`, `SpliceGirl`) whose scores and filters are not comparable; and per-gene semicolon lists in
`Gene A Location` / `Gene A Sense` that line up with the multi-gene partners.

Where the parsed rows land is [#1558](https://github.com/SACGF/variantgrid/issues/1558) (displaying
non-variants with variants) — this plan covers reading the file, not modelling it.

## 3. `SpliceVariants.vcf` — declare the filter, map the sample fields

Loads all 17 records, with two things that come through silently wrong.

**`LowUniqueAlignments` is dropped.** The header declares `LowQ` only, so `vcf_clean_and_filter`
strips the undeclared code from 15 of the 17 records and counts it (`Skipped 15 'LowUniqueAlignments'
FILTER`). Note the genotype processor itself already copes — `_add_undeclared_filter_code` creates a
`VCFFilter` on the fly — so the header stage is the only thing discarding it. Options are to let
undeclared FILTER values through `vcf_clean_and_filter` and rely on that, or to have the TSO 500
loader add the missing `##FILTER` line. Preserving the value matters here: it is the filter that
separates the two `PASS` oncology calls from the 15 background ones.

**`AD` and `DP` mean something else.** `configure_vcf_from_header` binds them as allele depth and read
depth by name. In this file `AD` is `Number=1` splice-supporting reads and `DP` is reference reads.
cyvcf2's `gt_ref_depths`/`gt_alt_depths` return `-1` because `AD` is not `Number=R`, so allele depth
lands empty and VAF ends up missing, while the read depth column shows the reference read count. With
no `GT` in FORMAT, every record is unknown zygosity.

**Change:** give the splice VCF its own field mapping so `AD` → alt depth and `DP` → ref depth, and
derive VAF from `ALTDEDUP`/`REFDEDUP` in INFO, which carry the deduplicated counts the caller
intends. `VCFSourceSettings` already keys off the `##source` line (`SpliceGirl 1.0.0.614`) and is the
natural place to hang this.

## 4. `cnv.vcf` — drop the copy-neutral rows, surface `SM`

Loads as-is; two refinements.

The 9 records with `ALT=.` (406 of 500 in a real run) import as reference variants — single base at
POS, segment span discarded. Filtering them at load keeps the noise out; they carry no call.

`SM` (linear copy ratio), `BC` and `PE` survive in the format JSON blob that importer v21+ stores, so
nothing is lost, but the copy ratio is not a first-class field. Absent `AD`/`DP`/`AF` the import logs
`WARNING Couldn't determine allele depth format field`, which is accurate and harmless. Promoting `SM`
to something queryable belongs with [#1558](https://github.com/SACGF/variantgrid/issues/1558) and
[sapath#304](https://github.com/SACGF/variantgrid_sapath/issues/304) — and note #304's open question
about fold-change conversion is already answered for v2.6.2, which emits `<DUP>`/`<DEL>` directly with
`SM` as the linear copy ratio.

---

## Order of work

1. cyvcf2 for declared filters in `vcf_clean_and_filter` — unblocks `_DragenExonCNV.vcf`, smallest change.
2. `get_processing_ability` returns 0 on unreadable input — unblocks uploading `AllFusions.csv` at all.
3. Fusion CSV parser + factory.
4. `VCFSourceSettings` field mapping for SpliceGirl; `LowUniqueAlignments` preserved.
5. Copy-neutral `cnv.vcf` records filtered at load.

Steps 1-2 are self-contained pipeline fixes worth doing regardless of the TSO 500 work. Steps 3-5 are
where the TSO 500 loader takes shape, and they feed the modelling issues:
[#1704](https://github.com/SACGF/variantgrid/issues/1704) (specimen → extraction → sample),
[#1558](https://github.com/SACGF/variantgrid/issues/1558) (non-variants shown with variants),
[#1559](https://github.com/SACGF/variantgrid/issues/1559) (specimen-level TMB/MSI),
[variantgrid_private#223](https://github.com/SACGF/variantgrid_private/issues/223) (PatientNode, to
analyse the DNA and RNA arms together).

## Still to obtain

`Logs_Intermediates/Gis/<sample>/<sample>.abcn_annotated.vcf` — per-gene absolute and minor copy
number, the latter being what gene-level LOH derives from. It is a VCF, so it is in ingestion scope,
but Illumina publishes no INFO/FORMAT specification for it and it cannot be usefully mocked up. See
the Missing section of `upload/test_data/tso500/README.md`.
