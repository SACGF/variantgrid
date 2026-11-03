# ACMG evidence the annotation captures, and the analysis node that filters on it (#2600)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-11-03
Status: in progress

Issue: https://github.com/SACGF/variantgrid/issues/2600

## What already exists

The dedicated filters the issue asks for landed in #2561 (PR #2568): the Effect node
(`analysis/models/nodes/filters/damage_node.py:DamageNode`) filters on AVI (`avi_phred_min`, via
`annotation/pathogenicity_predictions.py:TOOLS`), AlphaGenome splicing (`alphagenome_splicing_min`), SpliceVault
(`splicevault_out_of_frame`) and UTRAnnotator (`utr_high_impact`), each gated off the annotation version's
`visible_columns` and skipped with a warning where the version's backfill has not run. #2596 added the PTC / NMD choice
(`ptc_nmd`). The first step is confirming these on a GRCh38 analysis on this box: each filter present in the editor,
applied in the SQL (`manage.py profile_analysis_nodes --analysis <id> --rerun --explain`), and shown as a chip.

## The audit

| ACMG criterion | Evidence we annotate | Node today | Status |
|---|---|---|---|
| PVS1 null variant | `impact` HIGH, `nmd_escape_status` / PTC (`annotation/ptc.py`) | Effect: `impact_min`, `ptc_nmd` | covered |
| PVS1 LoF is a known mechanism | gnomAD LoF constraint (`genes/models/models_gnomad_gene_constraint.py:GnomADGeneConstraint`: `lof_oe`, `lof_oe_ci_upper` LOEUF, `lof_pli`); `GeneAnnotation.gnomad_oe_lof` is the point estimate only | none | **gap: this PR** |
| PS1 / PM5 same residue known pathogenic | ClinVar records are matched by variant, not by protein position | none | out of scope (needs a protein-position index over ClinVar; raise as a follow-up issue) |
| PM1 hotspot / functional domain | `interpro_domain`, `domains`, `cosmic_count` | Effect: `protein_domain`, `cosmic_count_min` | covered |
| PM2 / BA1 / BS1 population frequency | gnomAD AF, popmax, FAF95/99, 1kg, UK10K, TOPMed, internal counts | Population node | covered |
| PM4 in-frame indel / stop loss outside a repeat | `variant_class`, consequence; `repeat_masker` | Effect: variant class and impact; repeat region not filterable | **gap: this PR** |
| PP2 / BP1 missense constraint | `GnomADGeneConstraint.mis_oe`, `mis_z_score` | none | **gap: this PR** (same join as PVS1 constraint) |
| PP3 / BP4 computational | REVEL, AlphaMissense, CADD, BayesDel, ClinPred, MetaLR, VEST4, EVE, MutPred etc., SpliceAI, dbscSNV, AlphaGenome splicing, SpliceVault, UTRAnnotator, AVI, conservation | Effect node, Conservation node | covered |
| PP5 / BP6 reputable source | ClinVar significance, stars, conflicts; internal classifications | ClinVar node, Classifications node | covered |
| BP7 synonymous without splice effect | consequence, splice scores | Effect node filters for a splice score, not for its absence | covered by exclusion (an Effect node with `impact_min` keeps synonymous out); no change |
| BS2 observed in healthy adult | `gnomad_hom_alt`, popmax hom alt | Population: `gnomad_hom_alt_max` | covered |
| PS2 / PM6 de novo | trio genotypes | Trio node (de novo modes) | covered |
| PM3 / BP2 in trans / in cis | trio genotypes (compound het, unphased) | Trio / Quad / Duo comp het | covered (phase not annotated) |
| PP1 / BS4 segregation | pedigree genotypes | Pedigree node | covered |
| PP4 phenotype | HPO / OMIM / MONDO on patient and gene | Phenotype node | covered |
| Somatic (AMP tiers, oncogenicity) | ClinVar somatic tiers and oncogenicity, COSMIC | ClinVar node, Effect `cosmic_count_min` | covered |

The table, with the gaps resolved as "covered", goes into `analysis/__analysis_readme.md` as the section
"ACMG evidence and the node that filters on it" (that readme is three lines today), and is posted on the issue as the
🤖 comment that reports the work.

## The gaps this PR fills

All three are Effect node fields: a migration on `DamageNode`, editor form (`analysis/forms/forms_nodes.py`) and
template (`analysis/templates/analysis/node_editors/damagenode_editor.html`, static JS in
`variantgrid/static_files/default_static/js/analysis/damagenode_editor.js`), `_get_node_q`, `get_node_chips` and the method summary, in the
same shape as the existing `*_min` / `*_required` / `*_allow_null` triples.

### Gene constraint (PVS1 mechanism, PP2 / BP1)

```python
class DamageNode(AnalysisNode):
    ...
    gene_loeuf_max = models.FloatField(null=True, blank=True)          # LOEUF <= this keeps the variant
    gene_loeuf_required = models.BooleanField(default=False)
    gene_loeuf_allow_null = models.BooleanField(default=True)
    gene_mis_z_min = models.FloatField(null=True, blank=True)          # missense Z >= this
    gene_mis_z_required = models.BooleanField(default=False)
    gene_mis_z_allow_null = models.BooleanField(default=True)
```

The Q joins the representative annotation's symbol to the constraint table's MANE Select rows
(`annotation/effect_filters.py:gene_constraint_q`): gnomAD v4 has a row per transcript (~210k rows), and a MANE gene
has an Ensembl and a RefSeq MANE row with the same values, so a filter over every row would pass a tolerant gene for
any one constrained transcript. `variantannotation__symbol` is the text symbol, `GnomADGeneConstraint.gene_symbol`
the FK, both the HGNC symbol. Allow-null is a `NOT EXISTS` against the constraint table, not a negated `__in` across
the VariantAnnotation relation (that is a correlated scan of the partition, 35x slower). Null handling follows `*_allow_null`: with allow-null, genes with no constraint row
pass. Defaults on the slider: LOEUF 0.6 (gnomAD v4's constrained-gene cutoff; v2 used 0.35), missense Z 3.09; both are
recorded on the Pathogenicity Thresholds page the way `TOOLS` entries are, with the gnomAD flagship paper
(Karczewski 2020) as the LOEUF source with the gnomAD v4 cutoff note, and Lek 2016 (ExAC) for missense Z 3.09. The filter is offered only when the box has a
`GnomADGeneConstraint` import (`GnomADGeneConstraint.objects.exists()`, a `cached_property` on the node); otherwise
the editor hides it and the node warns.

### Repeat region (PM4)

```python
    outside_repeat_region = models.BooleanField(default=False)   # keep only variants with no RepeatMasker hit
```

Q: `Q(variantannotation__repeat_masker__isnull=True)`, ANDed (it is a restriction, not an OR'd evidence). Chip
"no repeat" with title "Outside RepeatMasker repeat regions". Offered when `repeat_masker` is in the version's
`visible_columns`, like the plugin columns.

## Tests

`analysis/tests/test_damage_node_v5.py` is where the Effect node's column filters are tested; new tests there: the
constraint Q keeps a variant in a constrained gene and drops one in an unconstrained gene, allow-null keeps a gene with
no row and required drops it; the repeat filter drops a variant with `repeat_masker` set. Fixtures: `GnomADGeneConstraint`
rows created directly (the model has no fake-data helper), `VariantAnnotation` rows written into the version's
partition with `library/django_utils/django_partition.py:temporary_db_table` as `analysis/tests/test_grid_genotype_join.py`
does. The migration's test is `manage.py makemigrations --check`.

## Docs

`analysis/__analysis_readme.md` gets the table. `annotation/AGENTS.md` gets nothing new unless the constraint join
turns out to need a line. The Pathogenicity Thresholds page carries the two constraint cutoffs and their source.
