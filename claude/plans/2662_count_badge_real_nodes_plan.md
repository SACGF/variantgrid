# Count badge click creates the real node; built-in filter nodes stop being made (#2662)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-11-05
Status: in progress

Clicking a node's count badge drops the real node configured to match the label - an Effect, ClinVar or
Classifications node - the way a tag badge already drops a `TagNode`. The parent-count reuse that was
`BuiltInFilterNode`'s one trick moves up to `AnalysisNode`, so the dropped node still shows the badge's number
without a count query. Each label's filter has one definition: the node configured for it. `BuiltInFilterNode`
stays only for the rows that exist (`disabled = True`, so the palette and the badge click never make another);
templates still carrying one lose it as people edit them.

## Data

### `DamageNode` (`analysis/models/nodes/filters/damage_node.py`) - two presence flags, migration 0013

```python
# "Has a COSMIC id", separate from cosmic_count_min: the id comes from VEP's Existing_variation on every
# deployment, the count only from the licensed COSMIC custom VCF (locally 5.1M annotation rows have an id,
# 0.79M a count)
in_cosmic = models.BooleanField(default=False)
in_cosmic_required = models.BooleanField(default=False)
# The variant's gene has OMIM phenotype terms (GeneAnnotation.omim_terms for the analysis's gene annotation version)
omim_phenotype = models.BooleanField(default=False)
omim_phenotype_required = models.BooleanField(default=False)
```

They behave as `protein_domain` / `published` do: an OR term in the scoring pool, `_required` moves it to the AND
set; listed in `modifies_parents`, `has_required`, `get_node_chips` (`COSMIC` / `OMIM` chips) and the editor
(`analysis/templates/analysis/node_editors/damagenode_editor.html`, rows beside protein domain). `get_omim_q`
moves from `analysis/models/nodes/node_counts.py` into `damage_node.py`, its only caller.

### `ClinVarNode` (`analysis/models/nodes/filters/clinvar_node.py`) - no fields change

`_get_node_q` ORs the narrowed axes instead of ANDing them: a record matches when any axis the node narrows says
so. An axis with every pill on stays out of the OR (it filters nothing, as `_filtering_values` already says), so
"germline P/LP, somatic untouched" still means germline P/LP; no narrowed axis means any ClinVar record. This is
what `ClassificationsNode._classification_match_q` does for germline/somatic, what `_review_status_q` already
does for `stars_min`, and what the C badge counts. The node is new (#2605) so no stored node changes meaning.

### `ClassificationsNode` (`analysis/models/nodes/filters/classifications_node.py`) - nothing changes

Each pill row is already gated on its own allele origin bucket (`_germline_q`, `_somatic_q`), so a node with
origin "All" and a narrowed germline row means germline P/LP or any somatic record, and a somatic row with every
pill off excludes somatic. The P and S badge nodes use origin Germline / Somatic: the same SQL as "All" with the
other row off, and the chips read that way.

### Count labels declared on the nodes

Each node class that can stand in for a badge declares the field values that make it that badge:

```python
class DamageNode(AnalysisNode):
    COUNT_LABEL_FIELDS = {
        BuiltInFilters.IMPACT_HIGH_OR_MODERATE: {"impact_min": PathogenicityImpact.MODERATE},
        BuiltInFilters.COSMIC: {"in_cosmic": True},
        BuiltInFilters.OMIM: {"omim_phenotype": True},
    }

class ClassificationsNode(SignificanceFilterNodeMixin, AnalysisNode):
    COUNT_LABEL_FIELDS = {
        BuiltInFilters.CLASSIFIED: {},
        BuiltInFilters.CLASSIFIED_PATHOGENIC: {"allele_origin": AlleleOriginFilterDefault.GERMLINE,
                                               "other": False, "benign": False, "likely_benign": False, "vus": False},
        BuiltInFilters.CLASSIFIED_TIER_1_2: {"allele_origin": AlleleOriginFilterDefault.SOMATIC,
                                             "tier_3": False, "tier_4": False},
    }

class ClinVarNode(SignificanceFilterNodeMixin, AnalysisNode):
    COUNT_LABEL_FIELDS = {
        BuiltInFilters.CLINVAR: {"germline_uncertain": False, "germline_likely_benign": False,
                                 "germline_benign": False, "germline_other": False,
                                 "somatic_tier_3": False, "somatic_tier_4": False,
                                 "oncogenicity_uncertain": False, "oncogenicity_likely_benign": False,
                                 "oncogenicity_benign": False},
    }
```

`AnalysisNode.COUNT_LABEL_FIELDS = {}`. Fields not named take the model default, so the I node is impact-only
and the G node is every pill on, any origin, parent matching, no labs.

## Code

### `AnalysisNode` (`analysis/models/nodes/analysis_node.py`)

- `count_label_node(analysis, label)` classmethod: the subclass whose `COUNT_LABEL_FIELDS` names the label,
  found with `get_subclasses(AnalysisNode)` at call time as `analysis/models/nodes/node_types.py` does (the
  concrete nodes import this module, so it can't import them), instantiated unsaved with
  `analysis=analysis, **fields`. `count_label_q(analysis, label)` is that node's `_get_node_q()`; a tag label
  and `TOTAL` are not its business.
- `matches_count_label() -> Optional[str]`: the label whose configured node this one equals - every field the
  subclass declares (`_meta.get_fields()` below `AnalysisNode`) equal, `node_input` parent-matching where the
  class has it. `None` for everything else.
- `_get_cached_label_count`: after the existing cloned/parent checks, if `matches_count_label()` is `label` or
  the asked label is `TOTAL`, return the single parent's stored count for the matched label. Edit the node so it
  no longer matches and it drops to a real count, as `BuiltInFilterNode` does today.

### `analysis/models/nodes/node_counts.py`

`get_extra_filters_q(analysis, label)` and its per-label branches go. Every caller has a node, so
`get_node_extra_filters_q`, `get_extra_filters_count` and `_aggregate_node_counts` call
`node.count_label_q(node.analysis, label)` (tags keep `get_analysis_tagged_variant_ids`). The imports of
`Classification`, `Pathogenicity`, `ClinVarOncogenicity`, `PathogenicityImpact`, `GeneAnnotation` and
`AnnotationVersion` leave with the branches.

`ClassificationsNode.get_labs` reads a reverse relation, which an unsaved instance can't: it returns
`Lab.objects.none()` when `self.pk is None`.

### Badge click (`analysis/views/views_json.py:create_extra_filter_child`)

The non-tag branch becomes `AnalysisNode.count_label_node(node.analysis, extra_filters)` with `x`, `y`,
`ready=False` set, saved, then `add_parent` as now. The node's name is its auto name; its chips show the
configuration.

### `BuiltInFilterNode` (`analysis/models/nodes/filters/built_in_filter_node.py`)

- `disabled = True`: out of the palette and `get_node_types_hash`, as `SelectedInParentNode` is. Loading,
  editing, cloning, templates and the serializer keep working for the rows that exist.
- `_get_node_q` is `self.count_label_q(self.analysis, self.built_in_filter)`, ANDed with the stars / COSMIC
  count extras the fields still hold.
- `matches_count_label()` returns `built_in_filter` when neither extra is set; the `_get_cached_label_count`
  override goes, the base does it.
- `get_help_text` says it is kept for existing analyses; the Effect, ClinVar and Classifications nodes are the
  replacements.
- `get_node_class_choices` (`analysis/models/nodes/node_types.py`) includes disabled classes so the node grid
  still labels old rows; the palette keeps using `get_nodes_by_classification`, which doesn't.

### Presets

The Effect editor's generic preset fill resets model fields to their defaults, so the two new booleans need no
entry in `PRESET_FIELD_DEFAULTS`.

## Counts users will see change

| Label | Before | After |
|---|---|---|
| P | any bucket with germline significance LP/P | germline or unknown bucket, LP/P |
| S | any bucket with somatic significance Tier I/II | somatic or unknown bucket, Tier I/II |
| C, I, O, M, G | | the same |

The P/S shift is the bucket now agreeing with the significance axis, which is the more correct count.
`annotation/models/models_cohort_stats.py` keeps its own T/I/O/C aggregates for sample and cohort nodes; their
definitions already match and are covered by the tests below.

## Tests

- analysis/tests/test_count_label_nodes.py (new file): fake annotation (impact, cosmic_id, OMIM gene
  annotation, ClinVar across the three axes) and a few published classifications per bucket; for every label in
  `BuiltInFilters.FILTER_CHOICES` the badge count on a node equals `count_label_node(...)` dropped under it, and
  the dropped node reports `matches_count_label() == label` and reuses the parent's count; editing one field
  returns `None`.
- `analysis/tests/test_clinvar_node.py`: two narrowed axes match a record satisfying either; an all-on axis
  adds nothing.
- `analysis/tests/test_classifications_node.py`: the P/S expectations above.
- `analysis/tests/test_built_in_filter_node.py`: the historical node still filters and reuses the parent count;
  it is absent from `get_nodes_by_classification`.
- The badge click view's tests (`analysis/views/views_json.py:create_extra_filter_child`): each non-tag label
  creates the right class, and a sample node's C badge (stats cache path) equals the ClinVar node's count.

## Docs

- `analysis/AGENTS.md`: one line under patterns - a count badge is a node configuration, declared in that node's
  `COUNT_LABEL_FIELDS`; `AnalysisNode.count_label_node` builds it and `matches_count_label` lets it reuse the
  parent's count.
- `variantgrid/templates/default_templates/changelog.html`: the badge now drops a real node; built-in filter
  nodes are retired.
