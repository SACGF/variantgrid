# Split ClassificationsNode into Classifications + ClinVar nodes (#1789)

## Goal

`ClassificationsNode` currently answers two questions in one editor — "what have our labs said about
this variant" and "what has ClinVar said" — and ORs them together internally. Split it into two nodes
that each answer one question, and let the analysis graph compose them (Merge for OR, chaining for
AND, "Parent variants that DO NOT match" for NOT).

Both nodes ship as **filter** nodes: they take a parent by default, and the user picks
`Matching variants (no parent)` on the Input dropdown to make one a source — the pattern `TagNode`
already uses.

The ClinVar fields landed today (migrations `0121`, `0122`) and are not on prod, so this is a
straight move rather than a data migration.

## Conventions

These two rules apply identically in both nodes, so the pill rows behave the same wherever the user
sees them.

**Every significance row defaults to all-on, and all-on means the row applies no filter.**
`_selected_values()` returns the ticked fields; when that count equals the row's field count the row
contributes nothing, otherwise it contributes `Q(field__in=selected)`. This is what
`ClassificationsNode._germline_q()` already does — ClinVar's rows adopt it.

Each ClinVar axis gets a catch-all pill at the end of the row for records ClinVar has given no call
on, so "all pills on" is literally "everything in ClinVar":

| Axis | Column | Catch-all pill covers |
|---|---|---|
| Germline classification | `clinvar.highest_pathogenicity` | `0` — drug response, risk factor, etc |
| Somatic clinical impact | `clinvar.somatic_tier` | `NULL` — no somatic call |
| Oncogenicity | `clinvar.highest_oncogenicity` | `0` — no oncogenicity call |

**A bare node means "this source has something to say".** `ClassificationsNode` with every pill on is
"classified in this database"; `ClinVarNode` with every pill on is "has a ClinVar record"
(`Q(clinvar__isnull=False)`, which the annotation-version SQL transformer scopes to the analysis's
`ClinVarVersion`). Setting Input to `Parent variants that DO NOT match` on a bare `ClinVarNode` gives
"not in ClinVar", which is how the old `ClinVarRecordFilter.ABSENT` reading is expressed.

## §1 Shared input enum

Rename `analysis.models.enums.ClassificationsNodeInput` to `NodeMatchInput` and use it for both
nodes — its labels are already source/filter generic:

```python
class NodeMatchInput(models.TextChoices):
    """ Whether the node is a source (matching variants) or filters its parent """
    MATCHING_VARIANTS = 'C', 'Matching variants (no parent)'
    PARENT_MATCHING = 'I', 'Parent variants that MATCH'
    PARENT_NOT_MATCHING = 'E', 'Parent variants that DO NOT match'
```

Update the references in `classifications_node.py`, `forms_nodes.py`, and the two test modules.

## §2 ClassificationsNode after the split

Move `analysis/models/nodes/sources/classifications_node.py` →
`analysis/models/nodes/filters/classifications_node.py`, updating
`analysis/models/nodes/{sources,filters}/__init__.py` and the imports in
`analysis/views/nodes/node_views.py`, `analysis/forms/forms_nodes.py` and the tests.

Fields it keeps: `node_input` (**default `PARENT_MATCHING`**), `allele_origin`, the six germline pills,
the four somatic tier pills, and `ClassificationsNodeLab`.

`_get_node_q()` becomes:

```python
def _get_node_q(self) -> Q:
    q = self._classifications_q()
    if self.node_input == NodeMatchInput.PARENT_NOT_MATCHING:
        q = ~q
    return q
```

`has_classification_filters()` keeps its job of telling `ClassificationsNodeView` whether to compute
the out-of-date message, and drives the "changes over time" warning in the template.

`get_node_chips()`, `_get_method_summary()` and `_clinvar_summary()` shed their ClinVar branches;
`_clinvar_summary`, `_clinvar_q`, `_clinvar_somatic_tiers`, `_clinvar_review_status_q`,
`has_clinvar_filters` and the three `FIELD_CLINVAR_*` maps move to §3.

Help text: `"Variants classified in this database. Filters its parent, or set Input to 'Matching
variants' to read the database directly."`

## §3 ClinVarNode (new)

`analysis/models/nodes/filters/clinvar_node.py`.

```python
class ClinVarNode(AnalysisNode):
    node_input = models.CharField(max_length=1, choices=NodeMatchInput.choices,
                                  default=NodeMatchInput.PARENT_MATCHING)
    allele_origin = models.CharField(max_length=1, choices=AlleleOriginFilterDefault.choices,
                                     default=AlleleOriginFilterDefault.SHOW_ALL)

    germline_pathogenic = models.BooleanField(default=True, blank=True)
    germline_likely_pathogenic = models.BooleanField(default=True, blank=True)
    germline_uncertain = models.BooleanField(default=True, blank=True)
    germline_likely_benign = models.BooleanField(default=True, blank=True)
    germline_benign = models.BooleanField(default=True, blank=True)
    germline_other = models.BooleanField(default=True, blank=True)

    somatic_tier_1 = models.BooleanField(default=True, blank=True)
    somatic_tier_2 = models.BooleanField(default=True, blank=True)
    somatic_tier_3 = models.BooleanField(default=True, blank=True)
    somatic_tier_4 = models.BooleanField(default=True, blank=True)
    somatic_tier_none = models.BooleanField(default=True, blank=True)

    oncogenicity_oncogenic = models.BooleanField(default=True, blank=True)
    oncogenicity_likely_oncogenic = models.BooleanField(default=True, blank=True)
    oncogenicity_uncertain = models.BooleanField(default=True, blank=True)
    oncogenicity_likely_benign = models.BooleanField(default=True, blank=True)
    oncogenicity_benign = models.BooleanField(default=True, blank=True)
    oncogenicity_none = models.BooleanField(default=True, blank=True)

    stars_min = models.IntegerField(default=0)
    conflicting = models.BooleanField(default=False, blank=True)
    conflicting_significance = models.TextField(null=True, blank=True)
    variation_ids = ArrayField(models.IntegerField(), default=list, blank=True)
```

`min_inputs`/`max_inputs` mirror `ClassificationsNode`: 0 for `MATCHING_VARIANTS`, else 1.

`allele_origin` decides which axes apply, matching `germline_enabled` / `somatic_enabled` on
`ClassificationsNode`: germline gates the pathogenicity row, somatic gates both the somatic tier and
oncogenicity rows. A gated-out row applies nothing regardless of its pills, so switching the toggle
takes effect immediately.

### `_get_node_q`

```python
def _get_node_q(self) -> Q:
    and_filters = [Q(clinvar__isnull=False)]

    if self.germline_enabled:
        if selected := self._filtering_values(self.FIELD_PATHOGENICITY):
            and_filters.append(Q(clinvar__highest_pathogenicity__in=selected))
    if self.somatic_enabled:
        if selected := self._somatic_tier_values():
            and_filters.append(self._somatic_tier_q(selected))
        if selected := self._filtering_values(self.FIELD_ONCOGENICITY):
            and_filters.append(Q(clinvar__highest_oncogenicity__in=selected))

    if self.stars_min:
        and_filters.append(self._review_status_q(ClinVarReviewStatus.statuses_gte_stars(self.stars_min)))
    if self.conflicting:
        and_filters.append(Q(clinvar__conflicting_clinical_significance__isnull=False))
    if self.conflicting_significance:
        and_filters.append(Q(clinvar__conflicting_clinical_significance__icontains=self.conflicting_significance))
    if self.variation_ids:
        and_filters.append(Q(clinvar__clinvar_variation_id__in=self.variation_ids))

    q = reduce(operator.and_, and_filters)
    if self.node_input == NodeMatchInput.PARENT_NOT_MATCHING:
        q = ~q
    return q
```

`_filtering_values()` is `_selected_values()` plus the all-on check — it returns `[]` when every pill
in the row is ticked, so the row drops out. The catch-all pills map to their sentinel values:
`germline_other` → `0`, `oncogenicity_none` → `0`.

`_somatic_tier_values()` returns the selected tiers, appending `TIER_1_OR_2` when Tier I or Tier II is
selected (a record recorded as "Tier I/II" might be either — the rule
`ClassificationsNode._somatic_q()` already uses). `_somatic_tier_q()` builds
`Q(clinvar__somatic_tier__in=selected)`, ORing `Q(clinvar__somatic_tier__isnull=True)` when
`somatic_tier_none` is ticked.

`_review_status_q()` is today's `_clinvar_review_status_q()`: it ORs the review-status column of each
axis the node is actually filtering on, and reads all three when no axis is filtering — matching
`ClinVar.is_expert_panel_or_greater` taking the max across them. The "is filtering" predicate becomes
the same `_filtering_values()` used above.

### Display

- `get_node_class_label()` → `"ClinVar"`; `get_node_class_icon()` → `NodeIcon(fa="fa-solid fa-book-medical")`
- `get_node_name()` → `"ClinVar"`, or `"Not ClinVar"` for `PARENT_NOT_MATCHING`
- `get_help_text()` → `"Variants ClinVar has a record for. Filters its parent, or set Input to 'Matching variants' to read ClinVar directly."`
- `get_node_chips()` — allele origin chip when set, one chip per selected value on any filtering axis
  (lift `ClassificationsNode._significance_chips` to a module-level helper in
  `analysis/models/nodes/node_display.py` so both nodes call it), plus a `★`-suffixed chip for
  `stars_min`
- `_get_method_summary()` — today's `_clinvar_summary()` text, with the exclude/record branches gone
- `auditlog.register(ClinVarNode)`

## §4 Migrations

One migration on top of `0123`:

1. `CreateModel` `ClinVarNode`
2. `RemoveField` × 21 for the `clinvar_*` fields on `ClassificationsNode`
3. `AlterField` `ClassificationsNode.node_input` for the new default

Stored `node_input` values on existing rows are untouched, so `ClassificationsNode`s already saved on
dev keep behaving as sources; the default change applies to newly added nodes.

## §5 Forms, views, serializers, templates

**`analysis/forms/forms_nodes.py`** — `ClassificationsNodeForm` keeps `lab`, drops every `clinvar_*`
field and the `clinvar_variation_ids` clean/save handling. New `ClinVarNodeForm` carries the
`variation_ids` `CharField` + `clean_variation_ids()` + `save()` (lifted as-is, minus the lab block),
with widgets `{'allele_origin': RadioSelect(), 'stars_min': StarsWidget(),
'conflicting_significance': TextInput(...)}`.

**`analysis/views/nodes/node_views.py`** — add `ClinVarNodeView(NodeView)`. `ClassificationsNodeView`
keeps its out-of-date-classifications context. Both are picked up automatically by
`get_node_views_by_class()`.

**`analysis/serializers.py`** — add `ClinVarNodeSerializer`; `AnalysisNodeSerializer.__subclasses__()`
registers it.

**Templates** — `classificationsnode_editor.html` keeps the input select, allele origin toggle, the two
significance rows, the lab select and the changes-over-time warning. New
`clinvarnode_editor.html` carries the three ClinVar rows plus stars / conflicting / variation IDs, and
a line of muted help saying a bare node keeps everything ClinVar has a record for.

The pill CSS (`.significance-pills`, `.pill-toggle`, `.row-label`) and the `updatePillToggles` /
`updateSignificanceRows` JS are now wanted by both editors — move the CSS into
`variantgrid/static_files/default_static/css/analysis_nodes.css` and the JS into a shared
`analysis/node_editors/_significance_pills.html` include that each editor pulls in.
`analysis_nodes.css` has no `.scss` source, so it is edited directly.

The `#classification-node-form-widgets` id the pill CSS is currently scoped to becomes a
`significance-pills`-scoped ruleset so both editors pick it up.

## §6 Tests

`analysis/tests/test_classifications_node_clinvar.py` becomes `test_clinvar_node.py`, rewritten
against `ClinVarNode`:

- a bare node is `Q(clinvar__isnull=False)`, and `PARENT_NOT_MATCHING` negates it
- all pills on for a row contributes nothing; deselecting one gives `Q(...__in=[remaining])`
- the catch-all pills select `0` / `isnull=True`
- Tier I selected pulls in `TIER_1_OR_2`; Tier III alone leaves it out
- `stars_min` alone reads all three review-status columns; with one axis filtering it reads that axis
- `variation_ids` form cleaning: separators, blank, non-numeric rejected

`test_classifications_node.py` and `test_classifications_node_clinvar.py`'s
`ClassificationsNodeInputTest` keep the input/parent-count coverage for the slimmed node.

`test_node_editors_render.py` — replace `test_classifications_node_clinvar_section` with the
equivalent against `ClinVarNode`, and add the editor-posts-back-valid round trip for `ClinVarNode`
alongside the existing `ClassificationsNode` one.

`test_clone_nodes.py` — confirm it covers `ClinVarNode` via its node-class sweep.

## §7 Verification

```bash
python3 manage.py test --keepdb analysis.tests.test_clinvar_node \
    analysis.tests.test_classifications_node analysis.tests.test_node_editors_render \
    analysis.tests.test_clone_nodes
```

Then in the app: add each node to an analysis from the Filter group of the add-node dropdown, confirm
each takes a parent by default and releases it when set to `Matching variants (no parent)`, and check
a Merge node over the two reproduces the old OR behaviour.
