# Gene amplifications and losses as gene-level variants

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-09
Status: in progress

Issue: [#1836](https://github.com/SACGF/variantgrid/issues/1836). Builds on the gene-level variant design in
`snpdb/gene_level_variants.py` (#1506) and is the "coordinate-free gene-level CNV" follow-on that
`claude/plans/tso500_overall_plan.md` deferred.

## Decision

A whole-gene copy number call is a gene-level Variant, exactly as a fusion is: identity is gene plus direction, the
copy ratio is per observation, and there is no coordinate. The segment coordinates the caller writes are the panel's
target window, not the event, and are not stored as a Variant. What decides that a `<DUP>`/`<DEL>` record is a
gene-level event is a segment field in the VCF that names the gene (`SEGID` in DRAGEN TSO500 CNV output), which the
import factory reads from the header.

Classifications are made against the gene-level variant, so the analysis grid, tags, discordance, export and
classification import all work with no special-case code, which is the reason for choosing this over the two
alternatives considered in the issue (a grouping record linking coordinate dups, or overlap queries at every lookup).

## Why the coordinates are not the identity

The DRAGEN TSO500 CNV caller reports one record per gene, on a segment fixed by the panel manifest. Every sample gets
the same segment for the same gene, so within one pipeline version the coordinate dup already matches across samples.
It fails as soon as anything else is involved: a manifest or DRAGEN version change moves the segment, another assay
reports a different span for the same event, and MO's existing classifications name the gene with no coordinates at
all. The segment is also not the gene. From the Illumina example fixture in
`upload/test_data/tso500/ExampleSample_2600000001/ExampleSample_DNA_2600000001C/ExampleSample_DNA_2600000001C.cnv.vcf`:

```
chr1  9770511  <DUP>  SVLEN=16595;SVTYPE=CNV;END=9787106;REFLEN=16595;SEGID=PIK3CD  GT:SM:BC:PE  ./1:1.49428:22:15,10
```

PIK3CD on GRCh37 spans 9711789-9789172 in Ensembl. The segment starts 59 kb inside the gene, at the first targeted
exon. A rule of "SV fully contains the gene" would produce nothing for this file, so the caller's gene claim in
`SEGID` is the signal, and geometric containment for callers that don't name genes stays a follow-on.

Treating the segment as a structural variant would also mislead downstream: `annotation/models/models.py:VariantGeneOverlap`
rows for every neighbouring gene inside a large target window, gnomAD SV and ClinVar CNV matches against manifest
coordinates, and a GRCh37/GRCh38 pair of Variants for a build-free event.

The DragenExonCNV file also names a gene (`GENE=BRCA1`) but on partial-gene, exon-level calls. It keeps importing as
coordinate SVs. Only the segment field, not any gene-naming field, triggers the gene-level path.

## Data

### Rename `FusionGeneId` to `GeneLevelId`

`genes/models/models_gene_level.py:GeneLevelId` (`FusionGeneId` as it was) is already the right thing: pk is the HGNC id where the gene has
one, a local number above one million otherwise, with `symbol_str` for anything that leaves the deployment. Only its
name and docstring say fusion. It moves to its own module with the fields unchanged:

```python
# genes/models/models_gene_level.py
class GeneLevelId(models.Model):
    """ The number a gene-level variant carries as Locus.position and inside its alt """
    CUSTOM_ID_START = 1_000_000

    symbol_str = TextField(unique=True, db_collation='case_insensitive')
    gene_symbol = ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    hgnc = ForeignKey(HGNC, null=True, on_delete=SET_NULL)
    genes = ManyToManyField(Gene, blank=True)
```

A `RenameModel` migration in `genes/`; the fusion migrations up to `0090` are on `origin/master` so none are edited.
`genes/models/models_gene_fusion.py:GeneFusion` keeps its `anchor` / `partner` FKs and their `related_name`s.

### The event

```python
# genes/models/models_gene_level.py
class GeneCopyNumberEventKind(models.TextChoices):
    GAIN = "G", "Gain"
    LOSS = "L", "Loss"

class GeneCopyNumberEvent(models.Model):
    """ Twin of GeneFusion: what a gene-level copy number Variant is """
    variant = OneToOneField('snpdb.Variant', on_delete=CASCADE)
    gene = ForeignKey(GeneLevelId, on_delete=PROTECT)
    kind = CharField(max_length=1, choices=GeneCopyNumberEventKind.choices)
```

The kind is GAIN rather than amplification because the threshold that makes a gain an amplification is the lab's, not
ours. Identity never encodes a threshold; the copy ratio stays on each sample's `CohortGenotype` where the copy number
column already reads it.

### The Variant

Same contig, ref and svlen as a fusion (`snpdb/gene_level_variants.py`). `Locus.position` is the `GeneLevelId` pk and
the alt is `<GAIN:HGNC:3236>` / `<LOSS:HGNC:3236>` (namespace `GENE` for a custom id). The alt repeats the gene so the
alt alone says what the variant is, and `clean()` checks position and alt agree, as `GeneFusion.clean` does.
`library/genomics/vcf_enums.py:GeneLevelSymbolicAlt` renames its unused `AMP` value to `GAIN`; nothing has written it.

### The VCF

```python
# snpdb/models/models_vcf.py:VCF
gene_level_segment_field = TextField(null=True)   # the INFO key gene-level events were read from, e.g. SEGID
```

Recorded for the VCF page, next to `copy_number_field`. The known spellings are a setting,
`VCF_GENE_LEVEL_SEGMENT_FIELDS = ("SEGID",)`, beside `library/genomics/vcf_enums.py:VCFConstant` copy number fields.
The decision has to be made when the file is claimed, before any `VCF` row exists, so it is not a per-VCF override.

## Design

### 1. Import

A new factory beside `upload/import_task_factories/import_task_factories.py:DragenTSO500AllFusionsImportTaskFactory`,
claiming a `.vcf` whose header declares one of the segment fields, with a processing ability above
`upload/import_task_factories/import_task_factories.py:GenotypeVCFImportFactory`. The header sniff is the same shape
as the AllFusions `can_process_file`.

Its pipeline is the fusion one:

- **Pre-VCF task** rewrites the file, the peer of the AllFusions create-VCF step in
  `upload/tasks/import_dragen_tso500_all_fusions_task.py`. Every called `<DUP>` / `<DEL>` record whose segment field
  resolves to a gene becomes one record on the gene-level contig: position and alt from the `GeneLevelId`, END = POS
  so svlen is 0, the sample's FORMAT columns copied through unchanged so `SM` is still there for
  `upload/vcf/vcf_import.py:get_copy_number_field` to bind. The caller's INFO is dropped. No-call rows (alt `.`) are
  not written. A called record with a blank segment field is not written either and is counted on the upload step,
  so a file this rule does not fit says so rather than half-importing. A name HGNC does not know is not a reason to
  drop a call: it gets a custom `GeneLevelId`, as an unknown fusion partner does.
- **Preprocess** is `upload/tasks/vcf/import_vcf_tasks.py:GeneLevelPreprocessVCFTask`, split only.
- **Post-insert** creates the `GeneCopyNumberEvent` rows for the VCF's variants, then `VCFCheckAnnotationTask`, as the
  fusion factory does. Zygosity counts are left off as for fusions.

Gene resolution reuses the single-gene half of `genes/gene_fusions.py:GeneFusionResolver` (matcher, then HGNC
previous and alias symbols, then a custom id), extracted so the two callers share it. `MYCL1` in the caller's output
resolves to `MYCL` this way.

A mixed file is not supported: the standard preprocess runs `bcftools norm` against a reference and cannot carry a
gene-level record, and the gene-level preprocess skips the reference stages. A VCF with a segment field is gene-level
throughout.

### 2. Event rows

`genes/gene_fusions.py:create_gene_fusions_for_variants` today takes every variant on the gene-level contig with no
`GeneFusion` and parses its alt as a fusion. It filters on the fusion alt kinds from now on. A sibling
`create_gene_copy_number_events_for_variants` does the same for GAIN / LOSS, and
`upload/tasks/vcf/import_vcf_tasks.py:GeneLevelInsertEventsTask` (`GeneLevelInsertGeneFusionsTask` as it
was) becomes the task that runs both, so the
classification-driven `GENE_LEVEL_INSERT_VARIANTS_ONLY` pipeline mints events as well as fusions.

### 3. Annotation

`annotation/gene_level_annotation.py:annotate_gene_level_run` already claims every variant on the gene-level contig.
It gains a copy number branch beside the fusion one, and `annotation/gene_level_annotation.py:GeneLevelIdResolver`
is renamed with the model. For an event: one representative `VariantAnnotation` for the gene, per-transcript rows for
the release's transcripts, a `VariantGeneOverlap` for the one gene, `symbol` and `overlapping_symbols` the gene,
`hgvs_c` / `hgvs_g` the canonical string. Consequence is the SO term `transcript_amplification` for GAIN and
`transcript_ablation` for LOSS, impact HIGH, variant class `library/genomics/vcf_enums.py:VariantClass`
COPY_NUMBER_GAIN / COPY_NUMBER_LOSS, which already exist and are what the grid's kind display keys on.

### 4. Canonical string, classification import and search

Fusions have `BCR::ABL1` as canonical output and accept `BCR-ABL1`. Copy number events get:

- Canonical: `EGFR amplification`, `EGFR loss`. These are the words on the reports and in Illumina's combined output.
- Accepted on input, case-insensitive, aliases applied: `amplification`, `amp`, `gain`; `loss`, `deletion`, `del`.
  "Deletion" stays input-only because as output it reads as a coordinate event.

The classification import's `resolve_gene_fusion` in `classification/models/classification_variant_info_models.py`
becomes `resolve_gene_level`, trying the fusion string then the copy number string, and sets the variant coordinate
the same way, so a historical classification with no coordinates goes through the existing gene-level insert pipeline
and mints the variant if no sample has brought it in yet. The `gene_fusion` property that feeds `gene_symbol` on the
resolved variant info gets a gene-level sibling that answers for either kind. The Allele side needs nothing: one
Allele serves every build and ClinGen is skipped, as `genes/gene_fusions.py:get_gene_fusion_allele` already does.

Search: a receiver beside `snpdb/signals/variant_search.py:search_variant_gene_fusion`, lookup only, so
`EGFR amplification` finds the variant and typing it never mints one.

### 5. Display

- `variantopedia/templates/variantopedia/variant_details.html` and the grid row detail show a "Copy number" row
  through a tag beside `genes/templatetags/gene_fusion_tags.py`, reading the event the way `_get_gene_fusion` reads
  the fusion.
- The grid's kind badge and the JS alt parser in `variantgrid/static_files/default_static/js/variantgrid_formats.js`
  take the kind from everything before the first colon, so `<GAIN:...>` needs no change there.
- `snpdb/management/commands/delete_unused_variants.py` adds `GeneCopyNumberEvent` to the derived rows it deletes with
  a variant; `GeneLevelId` rows stay, as they do for fusions.

### 6. Existing data on vg-test2

Three TSO500 CNV uploads (VCF pks 28, 32, 37) hold 16 coordinate `<DUP>`/`<DEL>` variants between them. Re-import
the three files once the new factory claims them; the orphaned coordinate variants fall to `delete_unused_variants`.
No conversion command for 16 variants.

### 7. Housekeeping

- `CACHE_VERSION` in `variantgrid/settings/components/default_settings.py` is bumped for the model rename.
- `snpdb/gene_level_variants.py` docstring: the "designed for but not yet built" sentence and the fusion-only wording.
- `claude/plans/tso500_overall_plan.md`: the follow-on bullet and the line saying cnv.vcf dups stay structural
  variants.
- `snpdb/CLAUDE.md` gene-level line, `upload/CLAUDE.md` for the segment field rule, `claude/domain.md` entries for
  gene-level id and gene copy number event.

## Out of scope

- **Geometric containment** for callers that write no gene name. Needs its own definition of whole gene (the example
  above shows gene span is wrong even for a panel that means the whole gene; "every coding exon of the canonical
  transcript" is the candidate), must be opt-in per source and gated by a gene list because an arm-level WGS gain
  contains hundreds of genes.
- **Per-event zygosity counts** ("how many samples carry EGFR gain"). Fusions don't have them either; one change for
  both when wanted.
- **Amplification threshold filtering** in analysis. The copy ratio is on the genotype and the copy number column shows
  it; a node filter on it is a separate piece of work.

## Tests

Mirror the fusion suites, keeping only tests of our logic:

- Factory claims the example cnv.vcf fixture over the generic VCF factory, and does not claim the example
  DragenExonCNV fixture (`GENE=`, no `SEGID`).
- Rewrite: called records land on the gene-level contig with the right alt and svlen 0; `SM` survives in FORMAT and
  binds as the copy number field; no-call rows are absent; a record with a blank segment field is absent and counted;
  an alias symbol resolves to the current one; a name HGNC lacks gets a custom id and is still written.
- Post-insert creates one `GeneCopyNumberEvent` per variant and is idempotent; `create_gene_fusions_for_variants`
  leaves a GAIN variant alone.
- Annotation writes the overlap for the gene, the consequence and variant class per kind, and a gene list on the gene
  finds the event.
- Classification import resolves `EGFR amplification` and `egfr amp` to the same variant coordinate and leaves
  `BCR::ABL1` resolving as a fusion; search finds the variant and mints nothing.

## For the next agent

The `1506` breakpoint-resolution plan is in progress in the same tables; land or rebase on it before the rename.
The alt enum value and the event kind are GAIN, and the human string is "amplification": that is deliberate, not an
inconsistency to tidy.
