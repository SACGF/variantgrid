# #1862 — HPO terms not pulling associated genes

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-11
Status: draft

[#1862](https://github.com/SACGF/variantgrid/issues/1862): every HPO term on a phenotype node warns "have no
associated genes" and the node passes all variants through.

## Cause

`ontology/management/commands/ontology_import.py:load_phenotype_to_genes` was written for HPO's old
phenotype_to_genes.txt - a `#Format:` comment line then seven columns (`hpo_id, hpo_name, entrez_gene_id,
entrez_gene_symbol, status, source, omim_id`) - and keeps only rows whose `source` is `mim2gene`. HPO changed the
file between mid 2022 and late 2023 to five columns with a plain header row:

```
hpo_id	hpo_name	ncbi_gene_id	gene_symbol	disease_id
HP:0004808	Acute myeloid leukemia	1050	CEBPA	OMIM:601626
HP:0004808	Acute myeloid leukemia	2672	GFI1	ORPHA:486
```

Reading it with the seven old names leaves `source` NaN on every row (the header row is read as data), the
`mim2gene` filter empties the frame, and the import completes "successfully" with zero relations. Reproduced on
this box: the July 2024 import owns 0 `OntologyTermRelation` rows and took under a second, so no HPO term reaches
a gene. Re-running with the same file is a no-op because `OntologyBuilder.ensure_hash_changed` sees a completed
import with the same hash.

The phenotype node is behaving as designed for a term with no genes (`analysis/models/nodes/filters/phenotype_node.py:PhenotypeNode.modifies_parents`
returns False, so it is a pass-through). Nothing changes there.

## Data

No model changes. Re-importing creates a new `OntologyImport` row (its own `OntologyTermRelation` partition),
and `ontology/models/models_ontology.py:OntologyVersion.latest` then creates a new `OntologyVersion` and
annotation sub-version. Relations written are the same two kinds as before: HPO -> OMIM `associated` and
OMIM -> HGNC `mim2gene`; ORPHA rows are skipped since `OntologyService.ORPHANET` terms are not stored locally.

## Steps

1. **Read either format** in `load_phenotype_to_genes`. Peek at the first line: a `#` comment is the old
   seven-column file (keep today's `names=` read and the `source == "mim2gene"` filter); otherwise read with
   `header=0` and rename `ncbi_gene_id` / `gene_symbol` / `disease_id` to the internal `entrez_gene_id` /
   `entrez_gene_symbol` / `omim_id` so the rest of the loader is shared. Filter to `omim_id.str.startswith("OMIM:")`,
   which is what `mim2gene` meant. Raise a clear error if the frame is empty after filtering - an import that
   writes nothing is the bug being fixed.
2. **Bump `processor_version` to 2** on that builder, so the previous completed import with the same hash no
   longer counts as "already done" and the file re-imports without `--force`.
3. **Test**: a five-column fixture in `ontology/tests/test_data/` (a dozen rows, both OMIM and ORPHA) loaded by
   `ontology/tests/test_data_ontology.py:create_ontology_test_data`, asserting one HPO term resolves to a gene
   symbol through `OntologyVersion.gene_symbols_for_terms`. That covers the format detection, the OMIM filter and
   the HPO -> OMIM -> HGNC path in one go.
4. **Deploy step** as a `ManualOperation.operation_other` in a new `ontology` migration (see
   `ontology/migrations/0008_load_ontology.py`), with a `test` that fires when the phenotype_to_genes import
   owns no relations: download the current phenotype_to_genes.txt, run
   `ontology_import --phenotype_to_genes`, then `gene_annotation --missing` for the new ontology version
   (`GeneAnnotation.hpo_terms` / `omim_terms`, which the node's free-text branch searches, come from the
   same relations), then move analyses to the new annotation version.
5. Note the format history in `claude/research/ontology.md` next to the existing `load_phenotype_to_genes` sentence.

## Verification

- `ontology/tests/test_ontology.py` and the new assertion pass with `--keepdb`.
- Locally: run the import against `/data/annotation/variantgrid_setup_data/update_dec_2023/phenotype_to_genes.txt`,
  then `OntologyVersion.latest().gene_symbols_for_terms(("HP:0004808",))` lists CEBPA, GFI1, BRCA2 among others,
  and a phenotype node with that term filters instead of warning.
