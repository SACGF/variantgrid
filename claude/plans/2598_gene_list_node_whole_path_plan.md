# GeneListNode tested through the real query path (#2598)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-11-03
Status: in progress

Issue: https://github.com/SACGF/variantgrid/issues/2598

## What is known

- The existing `analysis/tests/test_gene_list_node.py` applies `_get_node_q()` to `Variant.objects`. The node's real
  queryset starts from `annotation/annotation_version_querysets.py:get_variant_queryset_for_annotation_version`,
  whose SQL transformer rewrites every annotation table to the version's inheritance partition
  (`annotation_variantgeneoverlap` becomes `annotation_variantgeneoverlap_version_<n>`).
- Rows a test creates with the ORM (`VariantGeneOverlap.objects.create`) land in the base table, which the rewritten
  query never reads. A whole-path test written that way fails, with the gene list node returning nothing (exclude
  returning everything), for a fixture reason. Written into the partition the way production does
  (`library/django_utils/django_partition.py:bulk_create_in_partition`), the same test passes. The debugging version of
  this test, with print statements, is the uncommitted change to `analysis/tests/test_gene_list_node.py` in the main
  checkout at `/home/dlawrence/localwork/variantgrid`; it is the starting point.
- On this box (fromerd, dev database) every READY custom-text gene list node's live `get_queryset().count()` equals its
  stored count, so the live symptom ("doesn't seem to work") has not been reproduced here. The deployment it was seen
  on is not named in the issue.

## Design

### 1. A partition-writing fixture helper

`annotation/gene_level_annotation.py`'s private partition writer moves to `library/django_utils/django_partition.py` as `bulk_create_in_partition(version,
klass, base_table_name, records)` next to `temporary_db_table` (it uses only that and the version's
`get_partition_table`); `annotation/gene_level_annotation.py` imports it from there.

`annotation/fake_data.py` gains `create_fake_variant_gene_overlaps(variant_annotation_version, gene, variants,
annotation_run=None) -> list[VariantGeneOverlap]`, writing through `bulk_create_in_partition`. It is listed in
`claude/guides/testing.md`'s fixture index.

### 2. The whole-path tests

`analysis/tests/test_gene_list_node.py` keeps its editor and warning tests, and its filter tests become whole-path
tests: an `AllVariantsNode` parent set READY with a count, a `GeneListNode` child joined with `add_parent`, and the
assertion on `node.get_queryset()` pks (restricted to the fixture's variants). Fixture: one fake transcript version
(RUNX1, chr21, from `genes/fake_data.py:create_fake_transcript_version`), two variants overlapping it (written with the
helper above), one on chr21 outside it, one on chr1.

Cases, each a test method:
- named gene list: include, exclude; both with `analysis.node_queryset_filter_contigs` on and off (`_get_node_contigs`
  narrows to the gene's contig when on);
- custom text gene list: through `GeneListNodeForm` with `custom_gene_list_text`, then `node._load()` (which creates
  the hidden list), then the queryset;
- sample QC gene list (`accordion_panel = SAMPLE_GENE_LIST`, sample mode) and patient mode (`_get_sample_qc_gene_lists`);
- pathology test version;
- PanelApp panel with a `PanelAppPanelLocalCache` already holding the gene list;
- the full load: `node.load()` (or `update_analysis` with celery eager, as `analysis/tests/test_urls.py` runs it)
  stores `count` equal to the queryset size, and `NodeVersion.variant_ids` holds the pks.

### 3. Other tests writing partitioned annotation rows with the ORM

`analysis/tests/test_clone_nodes.py`, `analysis/tests/test_sv_comp_het.py`, `annotation/tests/test_vep_too_long.py`
and `genes/tests/test_gene_symbol_variants_grid.py` create `VariantGeneOverlap` (and some `VariantAnnotation` /
`VariantTranscriptAnnotation`) rows with the ORM. Each is checked: where the test reads through a partition-rewritten
queryset it switches to the helper (and its assertion is confirmed to exercise the rows); where it only reads the model
directly it stays as it is.

### 4. The live symptom

If a whole-path test fails for a code reason, the fix goes in this PR. If every path passes, the PR is the tests and the
helper, and the issue gets a comment saying so and asking which deployment showed the failure, so it can be looked at
there (its node's stored count against `get_queryset().count()`, and the node's errors, via
`manage.py shell` as in the probe above).

## Docs

One line in `annotation/AGENTS.md` under gotchas: rows created with the ORM in an inheritance-partitioned annotation
table land in the base table and are invisible to every annotation queryset; tests write them with
`create_fake_variant_gene_overlaps` / `bulk_create_in_partition`. `claude/guides/testing.md` lists the new fixture.
