# Analysis

Analysis allows users to take thousands upon thousands of variants and filter them down into meaningful ones, with the use
of an interactive UI.
## ACMG evidence and the node that filters on it

Which node filters on the evidence VariantGrid annotates for each ACMG/AMP criterion (#2600). The Effect node is
`analysis/models/nodes/filters/damage_node.py:DamageNode`; its cutoffs and their sources are on the Pathogenicity
Thresholds page (`annotation/pathogenicity_predictions.py:EFFECT_NODE_FILTERS`).

| ACMG criterion | Evidence we annotate | Node | Status |
|---|---|---|---|
| PVS1 null variant | `impact` HIGH, `nmd_escape_status` / PTC (`annotation/ptc.py`) | Effect: `impact_min`, `ptc_nmd` | covered |
| PVS1 LoF is a known mechanism | gnomAD v4 LOEUF (`genes/models/models_gnomad_gene_constraint.py:GnomADGeneConstraint` `lof_oe_ci_upper`, the gene's MANE Select row) | Effect: `gene_loeuf_max` (starts at 0.6) | covered |
| PS1 / PM5 same residue known pathogenic | ClinVar records are matched by variant, not by protein position | none | out of scope (needs a protein-position index over ClinVar) |
| PM1 hotspot / functional domain | `interpro_domain`, `domains`, `cosmic_count` | Effect: `protein_domain`, `cosmic_count_min` | covered |
| PM2 / BA1 / BS1 population frequency | gnomAD AF, popmax, FAF95/99, 1kg, UK10K, TOPMed, internal counts | Population node | covered |
| PM4 in-frame indel / stop loss outside a repeat | `variant_class`, consequence, `repeat_masker` | Effect: variant class, impact, `outside_repeat_region` | covered |
| PP2 / BP1 missense constraint | gnomAD v4 missense Z (`GnomADGeneConstraint.mis_z_score`, MANE Select row) | Effect: `gene_mis_z_min` (starts at 3.09) | covered |
| PP3 / BP4 computational | REVEL, AlphaMissense, CADD, BayesDel, ClinPred, MetaLR, VEST4, EVE, MutPred etc., SpliceAI, dbscSNV, AlphaGenome splicing, SpliceVault, UTRAnnotator, AVI, conservation | Effect node, Conservation node | covered |
| PP5 / BP6 reputable source | ClinVar significance, stars, conflicts; internal classifications | ClinVar node, Classifications node | covered |
| BP7 synonymous without splice effect | consequence, splice scores | Effect node keeps synonymous out with `impact_min`; it filters for a splice score, not its absence | covered by exclusion |
| BS2 observed in healthy adult | `gnomad_hom_alt`, popmax hom alt | Population: `gnomad_hom_alt_max` | covered |
| PS2 / PM6 de novo | trio genotypes | Trio node (de novo modes) | covered |
| PM3 / BP2 in trans / in cis | trio genotypes (compound het, unphased) | Trio / Quad / Duo comp het | covered (phase not annotated) |
| PP1 / BS4 segregation | pedigree genotypes | Pedigree node | covered |
| PP4 phenotype | HPO / OMIM / MONDO on patient and gene | Phenotype node | covered |
| Somatic (AMP tiers, oncogenicity) | ClinVar somatic tiers and oncogenicity, COSMIC | ClinVar node, Effect `cosmic_count_min` | covered |

The gene constraint filters read the representative transcript's `symbol` against gnomAD's MANE Select rows (gnomAD
v4 has a row per transcript, and a MANE gene has an Ensembl and a RefSeq row with the same values). They are offered
when the gnomAD gene constraint has been imported; Allow Null keeps a variant with no symbol or whose gene has no
MANE Select value. Outside repeat region is a restriction (ANDed), offered where the version populates `repeat_masker`.

The Effect, Population and ClinVar editors offer an ACMG preset dropdown (`analysis/node_presets.py`, #2605) that
fills the form with one criterion's cutoffs for the user to adjust and save: PVS1 (impact HIGH, LOEUF), PP3 (any
ClinGen-calibrated predictor at its supporting cutoff), PP2 (missense Z), PM2 and PP5 (ClinVar P/LP at two or more
stars). The PM2 lines depend on mode of inheritance, penetrance and severity and are lab policy, so they are
admin-edited `PopulationNodePreset` rows - seeded as dominant severe early onset (gnomAD allele count ≤ 1), dominant
adult onset (popmax / FAF95 ≤ 0.01%), recessive (≤ 1% and hom alt ≤ 2) and X-linked (≤ 0.01% and hemizygous ≤ 1) -
and the editor links superusers to them. Each preset is one criterion
because the criteria on a node pull in different directions (PM2 vs BA1, PVS1 vs PP3), and they look for pathogenic
evidence only. Nothing is stored on the node, so its name and chips describe what it filters on.
