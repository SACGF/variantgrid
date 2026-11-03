# Analysis node labels say what the node reads (#2599)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-11-03
Status: in progress

Issue: https://github.com/SACGF/variantgrid/issues/2599

## The rule

A node card has three places to say what it is: the name (`get_node_name`, up to three short lines on the card), the
badge icon (`get_node_icon`) and the chips under the name (`get_node_chips`, `analysis/models/nodes/node_display.py:NodeChip`,
each with hover `title` for the full text). The card is drawn from `analysis/models/nodes/node_utils.py:get_rendering_dict`.

- Line 1 of the name is the thing the node reads: the source object (sample, cohort, family, pedigree), the gene lists,
  the tag, the phenotype source. It is the most important fact and is never dropped for a count ("1 x PanelApp panels")
  or a mode.
- Line 2 is the mode or the one headline filter (inheritance mode, zygosity, "≥ Moderate impact").
- Everything else is a chip: one per detail that differs between two nodes of the same class, with the full value in
  `title`. A chip is omitted when the setting is at its default (not filtering).
- A name line longer than 30 characters is cut to 27 and "…" (`library/utils` has a truncating helper if one exists;
  otherwise a module-level function in `node_display.py` used by every node). When several things share line 1
  (several gene lists), line 1 is the first and a "+N more" suffix says the rest, which the chips list.

This rule goes in `analysis/AGENTS.md` (one bullet under Patterns, next to the node editor bullet).

## Changes the issue names

### Family nodes (Duo / Trio / Quad)

`analysis/models/nodes/family_inheritance.py:FamilyInheritanceNodeMixin.get_node_name` currently returns only the
inheritance label. New name: line 1 the family (`str(self.family)`: its name, or "Quad <pk>" when unnamed; truncated as above; the family model is the
one `FAMILY_SHAPE.source_field` names), line 2 the inheritance label with "?" when a required zygosity is unset, line 3
the filter description in parentheses as now. No family selected: the name stays the inheritance label.

### Effect node

`analysis/models/nodes/filters/damage_node.py:DamageNode.get_node_name` currently returns "1 of 15". New name, lines in
this order, only for settings that are set: impact ("≥ Moderate" from `PathogenicityImpact` label), predictions
("1 of 15 predictions" where `num_prediction_fields` is known, else "1 damage predictions"), PTC/NMD label when
`ptc_nmd` is set. Capped at three lines.

Chips (`get_node_chips`, which already covers variant class, splice, SpliceVault, 5' UTR): add one per raw-score tool
with a minimum set (text the tool's `node_label`, title "≥ value, required / allows null"), CADD phred, COSMIC count,
protein domain, published, and the "required" marker where `*_required` is on. The chips use the
`annotation/pathogenicity_predictions.py:TOOLS` metadata rather than a hand-kept list, so a new tool gets its chip.
`_get_method_summary` keeps returning the name.

### Gene list node

`analysis/models/nodes/filters/gene_list_node.py:GeneListNode.get_node_name`: line 1 the first list's name, truncated,
plus "+N more" when there are several; the "N x gene lists" / "N x PanelApp panels" fallback goes. For PanelApp,
`_get_gene_list_names` already resolves the cached name (with version) or the panel. Custom text keeps
"Custom: GENE1, GENE2…", sample gene list keeps "Sample Gene List: …", pathology test keeps "PathologyTest: …".
"Exclude: " prefix stays. A second line gives the PanelApp confidence where it is above LOW.

Chips: one per gene list or panel (text the name truncated to ~20 characters, title the full name, version and, for
PanelApp, the confidence), so the panel names are on the card even when there are several.

### Sample node: restricted to QC gene list

`restrict_to_qc_gene_list` on `analysis/models/nodes/sources/sample_node.py:SampleNode` has no mark on the card today.
The lock icon the issue remembers is found with `git log -S` over `analysis/` (search for `fa-lock`, `lock` in the
sample node's old css classes and the sample editor template `analysis/templates/analysis/node_editors/samplenode_editor.html`;
`#locked-warning-icon` in `analysis_nodes.css` and `LOCK_INPUT_SOURCES` in that template are the analysis-level lock,
a different thing). Whatever its old form, it comes back as a chip in `SampleNode.get_node_chips` when
`restrict_to_qc_gene_list` is on: icon `fa-solid fa-lock`, text "QC gene list", title the gene list's name, and
`css_class="node-chip-warning"` with title "Restricted to Sample Gene List, but none specified!" when
`sample_gene_list` is None (the same text `_get_configuration_errors` raises). If the old icon lived in the editor
rather than the card, it is restored there as well.

## The audit of every node class

Every concrete node (`analysis/models/nodes/filters/*.py`, `analysis/models/nodes/sources/*.py`; the list is in
`analysis/AGENTS.md`) is checked against the rule: for each, what `get_node_name` and `get_node_chips` return for a
configured node, read off real nodes on this box (`manage.py shell`: `AnalysisNode.objects.select_subclasses()` per
`node_type`, print `get_node_name()` and `get_node_chips()`) and from the code where the box has none. The outcome is a
table in the PR description: node class, line 1 today, line 1 after, chips added. A node whose line 1 is already the
thing it reads (SampleNode, CohortNode, PedigreeNode, TagNode, PhenotypeNode) is left alone; one that leads with a
count, a mode or nothing is changed to the rule, in the same style as the three above. Likely candidates: Zygosity,
AlleleFrequency, MOI (sample on line 1), Population (the threshold and which frequencies), ClinVar / Classifications
(significances as chips via `significance_chips`), Conservation, Tissue, BuiltInFilter, Filter, Intersection,
SelectedInParent, Merge / Venn.

## Tests

analysis/tests/test_node_names.py (new file): one test per changed `get_node_name` / `get_node_chips`, built on
`analysis/tests/utils.py:AnalysisSetupMixin` and the fake family / gene list fixtures in `claude/guides/testing.md`
(`snpdb/fake_data.py:create_fake_cohort` and friends for a Trio / Quad). The truncation helper has its own test.
Existing tests asserting old names (grep `get_node_name` and the old strings such as "x PanelApp") are updated.

## Verification

The canvas is `vg page /analysis/<id>` (names and chips come through `get_rendering_dict`), and
`scripts/browser/README.md` screenshots the cards against a working-tree `runserver`; a Quad analysis exists on this
box's dev database if the audit needs one to look at.
