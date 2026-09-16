# Issue #1875 — TSO 500 splice calls from the CombinedVariantOutput

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-16

Status: landed 404d06228 (phase 1 9fc889512, genome build fix 40af1b76e)

## Decision

Splice calls are loaded from the pair's CombinedVariantOutput TSV (the CVO) and from nothing else. The
SpliceGirl VCF (`upload/test_data/tso500/ExampleSample_2600000001/ExampleSample_RNA_2600000001B/ExampleSample_RNA_2600000001B_SpliceVariants.vcf`) is left out of the TSO 500 file set: its only extra content is the LowQ
background (15 of 17 records in the test file) and a QUAL score, and loading both would put the same call
in one sample twice, once as a genomic `<DEL>` and once as a splice event. The CVO's `[Splice Variants]`
section is already what a scientist reports from - Illumina writes only passing calls on AR, EGFR and MET
into it - and it names the gene and the two breakpoints instead of pretending to be a deletion.

The other CVO sections are not sources. Small variants and copy number stay on their VCFs, and fusions on
`upload/test_data/tso500/ExampleSample_2600000001/ExampleSample_RNA_2600000001B/ExampleSample_RNA_2600000001B_AllFusions.csv` (the CVO fusion rows are the `KeepFusion` subset with the caller, score, filters and
split/pair breakdown dropped, all of which `upload/tso500/dragen_all_fusions_parser.py` keeps).

A splice call becomes a gene-level `Variant` (`snpdb/gene_level_variants.py`), the same treatment as a
fusion or a whole-gene copy number call: no coordinate on a real contig, so no FASTA expansion, no VEP,
no gnomAD-SV and no AnnotSV class, and one `Allele` that means "AR-V7" rather than "a deletion of
chrX:66905968-66914514". The report's `^SpliceGirl` sniff on `sample.vcf.source`
(`classification/report/case_report_context.py`) goes, because the alt says what the variant is.

What this loses against a coordinate variant: `c.HGVS`, search by position, and liftover - all of which
a junction does have in principle. The report-side extract on the issue shows 0 of 40 reported splice rows
carry a transcript or `c.`, and every one carries a label (`V7`, `vIII`, `Exon 14 skipping`, in nine
spellings), so the label is what classifications will meet on, and the coordinates ride along in INFO.

Test data: `upload/test_data/tso500/ExampleSample_2600000001/ExampleSample_2600000001_CombinedVariantOutput.tsv`
(what is real and what is reconstructed is in `upload/test_data/tso500/README.md`).

## Data

### The alt

`library/genomics/vcf_enums.py:GeneLevelSymbolicAlt` gains a kind whose alt carries the event label as a
third segment, so two events in one gene are two variants:

```
<SPLICE:HGNC:644:V7>          AR-V7
<SPLICE:HGNC:3236:vIII>       EGFRvIII
<SPLICE:HGNC:7029:ex14skip>   MET exon 14 skipping
```

`GENE_LEVEL_ALT_PATTERN` accepts the extra segment for `SPLICE` only; a label is `[A-Za-z0-9._-]+`.
`parse` returns it as a fourth element (`None` for the other kinds) so `format`/`parse` stay one pair.

### Naming a junction

```python
class SpliceEvent(models.Model):
    """ A recurrent splice junction and the name a report gives it. Seeded; a lab adds rows for junctions
        its panel reports that we have not named """
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    label = models.TextField()                    # what the alt carries: V7, vIII, ex14skip
    display = models.TextField()                  # what the report writes: "AR-V7 splice variant"
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)
    donor = models.IntegerField()                 # Breakpoint 1: last base of the 5' exon
    acceptor = models.IntegerField()              # Breakpoint 2, as the caller writes it

    class Meta:
        unique_together = (("genome_build", "contig", "donor", "acceptor"),
                           ("gene_symbol", "label", "genome_build"))
```

Lives in `genes/` beside `genes/gene_fusions.py`. Seeded by data migration for GRCh37 from the test file's
three rows, and for GRCh38 by lifting the two points (AR exon 3 end / CE3, EGFR exon 1 / exon 8, MET exon 13 /
exon 15 - checked against the transcript's exon table in the migration, since these are exon boundaries).

A junction with no `SpliceEvent` row still imports, with a label made from the coordinates as written,
`X_66905968_66914514`. That label is build-specific by construction and reads as raw coordinates on the
report, which is the prompt to add a row. The `splice_label` evidence key
(`classification/migrations/0182_splice_label_ekey.py`) is autopopulated from `SpliceEvent.display` the way
the gene symbol is for gene-level events (15dd91a46), and left for the scientist when the label is a
coordinate one.

### What the sample carries

The rows become a VCF on the gene-level contig, as `upload/tasks/import_dragen_tso500_all_fusions_task.py`
does for fusions:

| | |
|---|---|
| sample | `RNA Sample ID` from `[Analysis Details]` (the splice caller runs on the RNA arm) |
| `##source` | `DRAGEN TSO500 CombinedVariantOutput <Module Version>` |
| INFO `SPLICE` | `AR V7` - gene and label, for the grid and search |
| INFO `SPLICE_OBS` | the CVO row as JSON: breakpoints, affected exon, both read counts |
| FORMAT `ALT_READS` / `REF_READS` | `Splice Supporting Reads` / `Reference Reads Transcript` |

A `VCFSourceSettings` row for `^DRAGEN TSO500 CombinedVariantOutput` binds the two FORMAT fields as alt and
ref depth, as `snpdb/migrations/0255_vcf_source_settings_fusion_processor_reads.py` does for fusions, so the
sample node's minimum-reads threshold and allele frequency work. The frequency it yields is the junction
ratio; the grid column keeps saying VAF, as it does for fusions.

The file carries no genome build. Declared at upload like the AllFusions CSV
(`upload/upload_metadata.py`); the create-VCF step needs it to look junctions up in `SpliceEvent`.

## Loader

- `upload/tso500/dragen_combined_variant_output_parser.py` - the file only, no database. A `[Section]`
  line starts a section; the next line is its header; rows follow to the blank line; every line is
  right-padded with tabs to the widest section and `NA` alone on a row means the section is empty.
  Returns `{section name: (header, rows)}` and tolerates a section it has no name for, since 2.6 renames
  `[Exon-Level CNVs]` to `Large Rearrangements` and adds `Gene-level Loss of Heterozygosity`. Recognised by
  its first line, `DRAGEN TruSight Oncology 500 Analysis Software - Combined Variant Output`.
- `upload/uploaded_file_type.py`: `DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT`, a VCF-loading type.
- `upload/import_task_factories/import_task_factories.py`: a factory beside
  `DragenTSO500AllFusionsImportTaskFactory`, same skipped bcftools stages.
- `upload/tasks/import_dragen_tso500_combined_variant_output_task.py`: create-VCF step reads
  `[Splice Variants]`, resolves `Gene` to a `GeneLevelId`, looks each (contig, donor, acceptor) up in
  `SpliceEvent` for the declared build, writes the VCF above. No post-insert step: a splice event has no
  `GeneFusion`-style row of its own, the alt and INFO carry everything.
- `annotation/`: `VariantAnnotationPipelineType.GENE_LEVEL` claims every gene-level variant, but it walks the
  `GeneFusion` and `GeneCopyNumberEvent` rows, and a splice call has neither - so
  `annotation/gene_level_annotation.py:_gene_level_events_for_run` also yields the splice events
  `genes/gene_splice.py:splice_event_variants` reads off the alts, with `splicing_variant` (SO:0001568) as their
  consequence and variant class. That gives them the `VariantAnnotation` and `VariantGeneOverlap` rows gene lists
  and the case report's gene lookup need.

## Display and report

- `classification/report/case_report_context.py:_kind_and_alteration` returns `SPLICE` when the parsed alt
  kind is `SPLICE`; `_is_splice_call` and `SPLICE_CALLER_SOURCE_PATTERN` are deleted.
- The grid's kind badge shows `SPLICE` and the event label, and a "Splice calls" column reads `SPLICE_OBS`
  off `CohortGenotype.info` the way the fusion column does (`analysis/models/nodes/cohort_mixin.py`).
- Variant page title for a gene-level splice variant is `<gene> <display>` (`AR-V7 splice variant`).

## Phase 2 - the rest of the file

Same parser, same upload step, after the splice calls land:

- `[Analysis Details]` names the whole chain: `Pair ID` is the patient (the Omico C-number,
  `Patient.patient_code`; the same one comes back on re-analysis), the ten-digit accession in the sample IDs
  is the specimen, and the container suffixes are the DNA and RNA extractions - what
  `seqauto/models/models_seqauto.py` matches from filenames today. Each level is resolved with
  `patients/external_references.py:resolve_reference` and created if absent, so a CVO arriving first
  leaves a stub Patient holding only its code, a Specimen and two named Extractions; the API keeps what
  the file lacks (name, DOB, sex, tissue, dates) and fills the stub in. In practice Mocha has pushed the
  patient before sequencing starts, so the stub is the exception the loader tolerates rather than the
  path it is built for. Both arms' samples are linked to their extraction by exact sample ID, which
  replaces the filename regex for this file set.
- The same two sample IDs are the seqauto join. DRAGEN writes them from the SampleSheet's `Sample_ID`, so
  each is an exact `SequencingSample.sample_name`, and the sample sheet gives the `SequencingRun`. The
  post-insert step writes `VCFFromSequencingRun` and `SampleFromSequencingSample` for the splice VCF
  against the RNA arm's sequencing sample, the rows `upload/vcf/vcf_import.py:link_samples_and_vcfs_to_sequencing`
  writes for a path-matched VCF. A sequencing sample not yet registered means no link rows, and the
  patient chain and measures still land. The upload's `path` is kept as every other TSO 500 file's is,
  for dedup and provenance; `create_backend_vcf_links` is not consulted for this file type, since the CVO
  is one file for two sequencing samples and nothing seqauto registers has that shape.
- `[TMB]`, `[MSI]` and `[GIS]` (tumour fraction, ploidy) are the five `SpecimenMeasureType` values, written
  through `patients/serializers.py:upsert_specimen_measure` against the specimen and DNA extraction.

## Tests

- parser: sections and padding, the `NA` forms, blank `Affected Exon`, an unknown section name.
- create-VCF task on the test file: three variants with the alts above, `ALT_READS:REF_READS` of
  `27:573`, `64:1`, `91:1`, sample named from `RNA Sample ID`, source from `Module Version`; a junction not in
  `SpliceEvent` gets the coordinate label.
- `GeneLevelSymbolicAlt` round trip with and without the label segment.
- report kind for a `<SPLICE:...>` variant, and `splice_label` autopopulation from `SpliceEvent.display`.
- The SpliceGirl tests in `upload/tests/vcf/test_vcf_processors.py:TestVCFProcessors` stay: the VCF is still importable on its own, it is just not part of the file set.
