# Whole-chromosome copy number calls as their own kind of Variant (#2908)

Written by Claude Fable 5.2 (claude-fable-5-2), 2026-11-25

Status: in progress - implemented on branch 2908_whole_chromosome_variants, in review; importing the example VCF on a deployment is not yet verified

Cause: nothing in the coordinate canonicaliser knows a whole chromosome from a focal SV, so a DRAGEN ploidy
record (`chr21 1 N <DUP> END=46709983`) is stored as an ordinary 46.7 Mb `<DUP>`, and another caller's
whole-chromosome segment (different start and end) is a different Variant again. Fixed at the cause: the
canonicaliser every insert and lookup goes through gives a whole-chromosome gain or loss one representation.

## What happens today

The example (`upload/test_data/vcf/dragen_ploidy/trisomy21.ploidy.vcf`) run through the preprocess pipe by hand
(`vcf_clean_and_filter | bcftools norm --check-ref=s | vcf_clean_alts`, GRCh38):

- 23 of 24 records are skipped as "reference call over a span (ALT='.' with END)" (`upload/management/commands/vcf_clean_and_filter.py`).
- 1 record survives unchanged: `NC_000021.9 1 . N <DUP> 150 PASS END=46709983;SVTYPE=DUP DC:NDC 52.7611:1.49987`.
  `REF=N` at `POS=1` passes `--check-ref=s` (the fasta base there is N).
- No `GT`, so the VCF has sample columns and no genotype: every zygosity is unknown, as for SpliceGirl.
- It routes to the STRUCTURAL_VARIANT pipeline, is over `ANNOTATION_VEP_SV_MAX_SIZE`, so VEP never sees it and it
  gets a `VEPSkippedReason.TOO_LONG` VariantAnnotation row plus a VariantGeneOverlap for every gene it spans
  (`annotation/vcf_files/import_vcf_annotations.py:insert_vep_too_long_skipped`, `genes/gene_overlaps.py:SVGeneOverlapResolver`).
  So gene lists and compound het already see every gene on the chromosome. AnnotSV has no size cap and is given it.
- The variant page's nearby tab on a 223 Mb `<DUP>` (vg-test2) holds a web worker for 44 s and then reports a timeout.

## The data

A whole-chromosome event is a coordinate Variant on its real contig, with the kind in the alt:

```python
class VCFSymbolicAllele:
    ...
    # A whole chromosome gained or lost (aneuploidy) - our own subtypes of the spec's DUP/DEL
    DUP_CHROM = "<DUP:CHROM>"
    DEL_CHROM = "<DEL:CHROM>"

WHOLE_CHROMOSOME_ALTS = (VCFSymbolicAllele.DUP_CHROM, VCFSymbolicAllele.DEL_CHROM)
```

| Field | Value |
|---|---|
| `Locus.contig` | the chromosome, per build |
| `Locus.position` | 1 |
| `Locus.ref` | the fasta base at position 1 (N on every standard human contig) |
| `Variant.alt` | `<DUP:CHROM>` or `<DEL:CHROM>` |
| `Variant.svlen` | `contig.length - 1`, positive for both |
| `Variant.end` | `contig.length` |

So there are exactly two possible whole-chromosome Variants per chromosome per build, and "everyone with a gain of
chr21" is that Variant's samples. GAIN / LOSS rather than trisomy / monosomy for the reason
`library/genomics/vcf_enums.py:GeneLevelSymbolicAlt` gives: the copy number is the sample's, not the variant's.
The sample's values (`DC`, `NDC`) stay in the VCF's own cohort FORMAT JSON as for any other FORMAT field.

No new model, no schema migration. A real contig (not the gene-level contig) because the event has a chromosome:
analysis nodes narrow by contig, region filters and IGV read position and end, and the existing too-long path
resolves its genes from the interval.

`settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES` gains the two alts, so they pass `vcf_clean_alts` when our own VCFs
(liftover, analysis export re-import) carry them. That also offers them as variant types on the All Variants page and
node ("Chrom gain" / "Chrom loss", `snpdb/variant_filters.py:VARIANT_TYPE_LABELS`), and they are "Structural".

### What makes a record whole-chromosome

`<DUP>` or `<DEL>` reaching to within `settings.VARIANT_WHOLE_CHROMOSOME_END_SLACK` (1 Mb) of each end of
`snpdb/models/models_genome.py:Cytoband.get_whole_chromosome_span`: the start of the first euchromatic band (`gneg` or
`gpos*`) to the end of the last, with chrY trimmed to the region between its PARs
(`CHRY_PSEUDOAUTOSOMAL_REGIONS`, masked with N in the analysis-set references callers use). `<DUP>` becomes
`<DUP:CHROM>`, `<DEL>` becomes `<DEL:CHROM>`. `<CNV>` has no direction and stays as it is; a contig with no cytobands
(MT, scaffolds, the gene-level contig) has no whole-chromosome form.

That span is what callers cover. DRAGEN's ploidy caller writes `1..length`. A segmenting caller starts at its first
usable bin: past the 10-60 kb telomere N on most chromosomes; on an acrocentric anywhere from the first non-N base to
the q arm (GRCh38 chr21 has sequence from 5.01 Mb, euchromatin from 13.0 Mb); on chrY between the PARs, ending at
the q12 heterochromatin (26.6 Mb). The slack is smaller than any terminal band but Yp11.32 (inside PAR1), so a call
that misses one - a chr21 `<DUP>` ending before q22.3 - stays a partial gain, as does any arm.

A `<DUP>`/`<DEL>` that misses the slack but covers `VARIANT_WHOLE_CHROMOSOME_NEAR_MISS_FRACTION` (90%) of the span
(`is_whole_chromosome_near_miss`, eg the 223 Mb chr2 `<DUP>` on vg-test2) is stored as written, and the VCF import
adds a `SimpleVCFImportInfo` (`WHOLE_CHROMOSOME_NEAR_MISS`) telling the uploader to ask an administrator if it is the
whole chromosome, and reports it to Rollbar. Raising the slack in the deployment's settings and running
`one_off_whole_chromosome_variants` moves such a variant.

## Changes

### 1. Canonical form (snpdb)

- `library/genomics/vcf_enums.py`: the two alts and `WHOLE_CHROMOSOME_ALTS`.
- `snpdb/whole_chromosome_variants.py` (new file): owns the rule and the reasoning. It imports nothing from
  `snpdb.models` (models_variant imports it, and anything under `snpdb.models` would be an import cycle), so the rule is
  pure: `whole_chromosome_alt(alt, position, end, span)` returns the whole-chromosome alt or None, and
  `is_whole_chromosome_near_miss` the borderline test. The span comes from
  `snpdb/models/models_genome.py:Cytoband.get_whole_chromosome_span`, the band query cached per contig.
  `snpdb/models/models_variant.py:VariantCoordinate.whole_chromosome(genome_build, chrom, alt, ref=None)` builds the
  canonical coordinate (ref read from the fasta unless the record was already at position 1). The module also holds
  the search pattern and `parse_whole_chromosome_search`.
- `snpdb/models/models_variant.py:VariantCoordinate.as_internal_canonical_form`: a symbolic `<DUP>`/`<DEL>` that the
  rule matches returns the canonical whole-chromosome coordinate; a whole-chromosome alt already there is returned
  re-anchored to position 1 / contig length (idempotent). This is the one place: `snpdb/variant_pk_lookup.py:VariantPKLookup`
  (every VCF import), manual entry and classification matching all call it.
- `VariantCoordinate.is_whole_chromosome` and `Variant.is_whole_chromosome` properties, and
  `Variant.get_whole_chromosome_q()`, an `alt__in` Sequence subquery (the #1887 pattern, as
  `snpdb/variant_filters.py:get_structural_variant_q`).
- `VARIANT_SYMBOLIC_PATTERN` / `VariantCoordinate.from_string` accept the two alts so the stored string round-trips.
- `can_be_made_explicit`, `symbolic_hgvs_interval`, `as_vcf_coordinate` and ClinGen registration already whitelist
  DEL/DUP/INV, so a whole-chromosome variant has no g.HGVS, no c.HGVS and no ClinGen allele, like `<CNV>`. Confirm
  each with a test rather than adding a branch.

### 2. Annotation

- `annotation/annotation_pipeline_routing.py:pipeline_type_variant_q` and `pipeline_type_for_alt`: whole-chromosome
  stays STRUCTURAL_VARIANT and is subtracted from ANNOTSV (AnnotSV cannot parse the alt, and gene-level annotation of
  a whole chromosome says nothing).
- `annotation/annotation_version_querysets.py:vep_too_long_q` / `is_vep_too_long` (queryset and per-variant twins,
  used by `filter_vep_sv_max_size`, `import_vcf_annotations.py` and the Variant Types page): a whole-chromosome alt is
  always on the "left out of the VEP dump, gene overlaps resolved locally" side, whatever `ANNOTATION_VEP_SV_MAX_SIZE`
  is (VEP cannot parse the alt). So `VEPRunner.get_variants_qs` filters the SV dump even when the setting is unset.
- `annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:BulkVEPVCFAnnotationInserter.add_sv_gene_overlaps`: for
  a whole-chromosome variant write the VariantGeneOverlap rows and leave `overlapping_symbols` empty - thousands of
  symbols in one grid cell is noise, and the kind badge already says what the row is.
- `annotation/transcripts_annotation_selections.py`: the variant page's message for a whole-chromosome variant says
  "Whole chromosome - not annotated by VEP. Genes assigned by overlap (N genes)" instead of the length message.
- `variant_class`: `VariantClass.COPY_NUMBER_GAIN` / `COPY_NUMBER_LOSS` on the skipped row, as gene-level copy
  number does (`annotation/gene_level_annotation.py`). The skipped-row insert writes `variant_class` beside
  `vep_skipped_reason` for that.

- Counted as SVs elsewhere: `annotation/tasks/calculate_sample_stats.py` counts `<DUP:CHROM>` with insertions and
  `<DEL:CHROM>` with deletions (variant class copy_number_gain / copy_number_loss), and `beacon/query_targets.py`
  queries them as EFO copy number gain / loss.
- `library/genomics/vcf_enums.py:VCFConstant.COPY_NUMBER_FIELDS` gains `NDC` (a ratio), so the ploidy VCF's depth
  against the expected ploidy is the sample's copy ratio in grids, exports and SampleNode's gain/loss thresholds.

### 3. Liftover

A whole chromosome maps to the same chromosome in the other build with no tool:
`snpdb/liftover.py:_liftover_using_dest_variant_coordinate` yields the destination build's canonical coordinate
first for an allele whose variant is whole-chromosome, under a new `AlleleConversionTool.WHOLE_CHROMOSOME`
("Same chromosome", choices-only migration), and yields an error for that tool when the destination build has no
contig of that name. No other tool is tried for one (ClinGen, dbSNP and the BCFtools source lane all return early). It
rides the existing destination-coordinate lane (the one ClinGen uses), so the Variant is created by the liftover VCF
import like any other.

### 4. Display, search, nearby

- Kind badge (`VariantGridFormat.variantKind` in `variantgrid/static_files/default_static/js/variantgrid_formats.js`):
  `<DUP:CHROM>` / `<DEL:CHROM>` show `GAIN` / `LOSS` with the title "Whole chromosome gain/loss", in their own solid
  css class beside the dup/del ones in the scss. `representativeVariantLabel` names the row by its chromosome, so it
  reads "chr21 [GAIN]" as a gene-level call reads "RAF1 [GAIN]", and the row's text is a search that finds it.
- Variant Types page (`variantopedia/variant_types.py`): a "Whole chromosome gain / loss" kind with its outcome
  (symbolic, SV pipeline, genes by overlap, liftover by same chromosome).
- Search (`variantopedia/signals/variant_search.py:search_variant_whole_chromosome`, beside the variant coordinate
  searches, registered in `variantopedia/apps.py:VariantopediaConfig.ready`): "trisomy 21", "monosomy X",
  "chr21 gain", "gain of chr21", "chrX loss" resolve to the whole-chromosome Variant in each of the search's builds
  (trisomy -> gain, monosomy -> loss). Lookup only. The result's samples grid is the "everyone with trisomy 21" list.
- Nearby tab (`variantopedia/views.py:nearby_variants_tab`): new `settings.VARIANT_DETAILS_NEARBY_MAX_VARIANT_LENGTH`
  (10 Mb). A longer variant renders the tab with "Nearby variants are not shown for a variant longer than 10 Mb"
  and runs no query. This covers whole-chromosome variants and every other SV that today waits out
  `VARIANT_DETAILS_NEARBY_STATEMENT_TIMEOUT_SECONDS`.

### 5. Stored variants

manage.py one_off_whole_chromosome_variants (new command, `category = "one-off"`, shipped as a `ManualOperation`):
finds stored `<DUP>`/`<DEL>` Variants the rule matches (candidates narrowed in SQL by `abs(svlen)` over half the
contig length) and, per variant, moves it to the canonical locus / alt / svlen / end in place, keeping its pk so
genotypes, tags and classifications follow. Where the canonical Variant already exists it reports the pair and
leaves both. `--dry-run` lists what it would change. The pipeline type is unchanged (SV, too long), so annotation
needs no re-route; the moved variant's existing VariantAnnotation row keeps its `overlapping_symbols` and null
`variant_class` (only new annotation writes the whole-chromosome form). The ManualOperation is in
`snpdb/migrations/0019_whole_chromosome_variants.py`, registered only where such a candidate exists. vg-test2 has none.

### 6. Tests kept

- snpdb: the rule and near miss (DRAGEN full span; chr21 from its first non-N base and from 21q; a chr21 `<DUP>`
  stopping short of q22.3; a long arm; chrY between its PARs; chr2 ending short, a near miss; `<CNV>`; contig without bands),
  canonical form idempotent, string round trip, no HGVS / ClinGen.
- annotation: routing (SV yes, AnnotSV no), always-too-long predicate, overlaps written with empty `overlapping_symbols`.
- snpdb liftover: destination coordinate from the same chromosome; missing contig is an error.
- upload: `trisomy21.ploidy.vcf` through `vcf_clean_and_filter` keeps the one record, and the processor test
  (`upload/tests/vcf/test_vcf_processors.py:TestVCFProcessors` pattern) resolves it to the `<DUP:CHROM>` coordinate
  with NDC as the copy number field; `whole_chromosome_near_miss.vcf` gets the import message and Rollbar report.
- variantopedia search (`test_search.py`): "trisomy 21" finds it. variantopedia: nearby tab over the length cap runs no
  query. Variant page: the whole-chromosome warning (`annotation/tests/test_vep_too_long.py`).
- The one-off command's test lives with it.

### 7. Docs in the same PR

`snpdb/AGENTS.md` (pattern line beside the symbolic alt one), `upload/AGENTS.md` (a DRAGEN ploidy VCF is an ordinary
genotype VCF: euploid records are reference spans, an aneuploid one becomes a whole-chromosome Variant, a euploid
sample's file imports with no variants), `annotation/AGENTS.md` (routing sentence), `claude/domain.md` (glossary:
whole-chromosome variant; variant kind). The changelog is not edited here: its cards are written per deploy from the
git log, and a deployed card is never edited.

## Decisions on the issue's open questions

- **Linking**: nothing new. The ploidy VCF is another VCF from the backend upload, its own Sample, tied to the
  extraction / specimen / patient by the same upload metadata as the sample's other VCFs.
- **Findable and classifiable**: yes, as the Variant above.
- **Other callers**: normalised by the band rule in the canonicaliser, so the stored form is caller-neutral and no
  DRAGEN-specific import path exists.
- **General guard**: the canonical form plus the AnnotSV subtraction and the nearby length cap bound the cost for any
  source. A chromosome arm is still a focal SV.
- **Mosaicism and the euploid chromosomes' depths**: a mosaic and a full trisomy are the same Variant; the called
  record's NDC is the sample's copy ratio, so SampleNode's ratio thresholds tell them apart (about 1.15 vs 1.5). The
  euploid records are reference spans and stay dropped, so per-chromosome NDC below DRAGEN's calling threshold is not
  stored. That, the header's
  `estimatedSexKaryotype` / `referenceSexKaryotype` / `autosomeDepthOfCoverage`, and a sex-mismatch flag are a
  per-sample QC record (the #1904 shape), a separate piece of work from the variant representation.
- **XXY**: a chrX gain finds XXY and XXX alike; telling them apart needs the sex karyotype above.

## Verification left for a person

Importing the example on a deployment (it writes `snpdb_variant`) and `vg page --js-writes-to-database` on a grid
showing the badge are ask-first on vg-test2, so they are testing notes on the issue.
