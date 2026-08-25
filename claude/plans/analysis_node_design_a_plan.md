# Analysis node redesign — "Design A" cards, FA/SVG icons, generic chips

Visual spec: [`claude/analysis_node_icons_v2.html`](../analysis_node_icons_v2.html) (open in a browser —
the "full palette in Design A" section is the target look; card CSS values can be lifted from it).
Decisions already made there: per-node icon picks, class-strip labels ("AF%" for allele frequency),
globe for Population, bands+arrows chromosome for Intervals, 3-lane merge glyph, Sample keeps the
pedigree sex-shape as its badge with extraction/VCF chips. The Merge node's top drop-surface mechanics
are explicitly **out of scope** (being handled separately).

## Current architecture (what we're replacing)

A node on the canvas is `div.window.<ClassName>` whose look comes entirely from CSS
(`variantgrid/static_files/default_static/css/analysis_nodes.css`):

- `.window.<ClassName>` rules set per-class fixed sizes (`:258-357`), `.node-overlay.<ClassName>`
  rules set PNG `background-image`s (`:264-374`), `.ui-icon.<ClassName>` rules (`:149-244`) give the
  add-node dropdown its icons via the same PNGs.
- `get_rendering_dict()` (`analysis/models/nodes/node_utils.py:125`) serializes a node as
  `{attributes, node_class, overlay_css_classes, name, args}`; `args` comes from
  `get_rendering_args()` (`analysis_node.py:313`, default `{}`) and is passed to the client-side
  `updateState()`.
- `analysis_nodes.js` builds DOM: `createDefaultNode()` (`:14`) makes `.node-overlay` + `.node-name`;
  `NODE_FACTORIES` (`:51`) special-cases SampleNode/VennNode/TrioNode/QuadNode with bespoke renderers
  (D3 sex-shapes in `samplenode.js`, pedigree drawing in `pedigree_node.js`, live venn widget).
  `updateNodeFromData()` (`:70`) re-applies `overlay_css_classes`/name/args on appearance updates
  (`updateNodeAppearance` `:934`, versioned by `appearance_version`).
- Counts float in an absolutely-positioned `.count-overlay` (`analysis_nodes.js:550`,
  css `:517-547`).
- The dropdown is `AnalysisNodeClassesForm` (`analysis/forms/forms.py:77`, choices from
  `get_nodes_by_classification()` in `analysis/models/nodes/node_types.py:41`) rendered at
  `analysis/templates/analysis/analysis.html:290`, iconified by `setupNodeTypeSelect()`'s
  `_renderItem` (`analysis.js:215`) which emits `<span class="ui-icon <ClassName>">`.

`samplenode.js:80` (`sourceBadge`) is the precedent for the chips idea — a one-off VCF-count badge;
this plan generalises it and retires the one-off.

## Target architecture

One generic card, everything data-driven from the rendering dict. Per-class knowledge lives on the
node classes (Python), not in CSS selectors.

### 1. Server-side display metadata

New module `analysis/models/nodes/node_display.py`:

```python
@dataclass(frozen=True)
class NodeIcon:
    fa: Optional[str] = None       # e.g. "fa-solid fa-filter"
    symbol: Optional[str] = None   # SVG sprite symbol id, e.g. "node-icon-gene-list"

@dataclass(frozen=True)
class NodeChip:
    text: str
    icon: Optional[str] = None     # FA classes
    title: Optional[str] = None    # hover text
    css_class: Optional[str] = None
```

On `AnalysisNode`:

- `get_node_icon(self) -> NodeIcon` — instance method so SampleNode can vary it; base returns the
  class default. Every concrete class overrides (picks per the v2 mock: `fa-database`, `fa-users`,
  `fa-clipboard-check`, `fa-filter`, `fa-bolt`, `fa-chart-simple`, `fa-earth-oceania`,
  `fa-stethoscope`, `fa-lungs`, `fa-sliders`, `fa-tags`, `fa-square-check`; sprite symbols for
  gene-list, zygosity, moi, intervals, conservation, merge, trio, quad, pedigree, venn, sample
  shapes).
- `get_node_class_label_short() -> str` — class-strip text; default `get_node_class_label()`.
  Overrides where the long label doesn't fit: AlleleFrequencyNode → `"AF%"`,
  IntersectionNode → `"Intervals"`, BuiltInFilterNode → `"Built-in"`,
  SelectedInParentNode → `"Selected"`, MOINode → `"MOI"`, TissueNode → `"Tissue"`.
- `get_node_chips(self) -> list[NodeChip]` — default `[]`. Initial implementations:
  - **SampleNode**: extraction chip (`fa-vial`, extraction id, title = full extraction description)
    when `source_level == EXTRACTION`; one `fa-file-lines` "VCF" chip per entry of
    `get_source_vcf_names()` (title = vcf name), collapsing to a single "VCF ×N" chip above ~3.
    Replaces the ad-hoc `sourceBadge` in `samplenode.js`.
  - **DamageNode**: `fa-scissors` "splice" chip when splicing filters are enabled — replaces the
    `EffectNodeSplicing` CSS-class/alternate-PNG mechanism (`damage_node.py` sets it via
    `get_css_classes`; css `:363`).
  - **GeneListNode** (optional, same PR or follow-up): convert the `span.gene-list-security` lock
    overlay (css `:539`) into a `fa-lock` chip.

`get_rendering_dict()` gains: `"icon": asdict(node.get_node_icon())`,
`"class_label_short": ...`, `"chips": [asdict(c) for c in node.get_node_chips()]`.
Chips derive from persisted node config, so the existing `appearance_version` flow already
invalidates them (same staleness window as today's `sourceBadge` for out-of-band VCF changes).

SampleNode's `get_node_icon()` returns `node-icon-sample-male` / `-female` (+ `-deceased` variants)
from `self.sample.patient`; `get_rendering_args()` then no longer needs to ship `patient`/`source`
and can revert to the base `{}`.

### 2. SVG sprite

`analysis/templates/analysis/tags/node_icon_sprite.html`: one hidden `<svg>` of `<symbol>`s
(viewBoxes and paths copied from the v2 mock), strokes/fills in `currentColor` so the same symbol
renders white on the badge and dark in the dropdown. `{% include %}` it from `analysis.html`.
Badge markup is `<svg class="node-icon"><use href="#node-icon-x"/></svg>`; FA icons are
`<i class="...">`.

### 3. Dropdown

Extend the per-class dicts in `get_nodes_by_classification()` (`node_types.py:48`) with
`icon` (class-default `NodeIcon`, via a `get_node_class_icon()` classmethod that
`get_node_icon()` defaults to) and `class_label_short`. Emit the whole structure as
`window.NODE_TYPES` via `json_script` in `analysis.html`, and rewrite `_renderItem` in
`setupNodeTypeSelect()` (`analysis.js:215`) to build the FA `<i>`/sprite `<use>` from it instead of
`<span class="ui-icon ...">`. Delete the `.ui-icon.<ClassName>` CSS block (`analysis_nodes.css:130-244`).

### 4. Client card (Design A)

Rewrite `createDefaultNode()` to:

```html
<div class="window design-a-node">
  <div class="node-overlay">              <!-- class list still driven by overlay_css_classes -->
    <span class="node-badge"></span>      <!-- corner circle, icon injected -->
    <span class="node-klass"></span>      <!-- class strip: class_label_short -->
    <div class="node-name-holder"><span class="node-name"></span></div>
    <div class="node-chips"></div>
    <div class="node-color-overlay"></div>
  </div>
  <div class="node-counts-strip"></div>
</div>
```

`updateNodeFromData()` additionally populates badge icon, `node-klass` text, chips (a small
`renderChip(chip)` helper: icon + text + `title`), and sets the card's `title` to the full name
(body clamps to 2 lines with CSS `line-clamp`; hover gives the long text). Keep the
`node-overlay` element + `overlay_css_classes` contract untouched so `updateNodeAppearance`
(`:934`) and the output/variable colour overlays (`variable-node`/`output-node`,
css `:499-505`) keep working.

- **Uniform sizing**: one `.window.design-a-node { width: 128px; }`, height auto. Delete all
  per-class `.window.<X>` size rules. Endpoints are anchored `TopCenter`/`BottomCenter`
  (`analysis_nodes.js:621-641`), so call `jsPlumb.revalidate(node)` at the end of
  `updateNodeFromData` — chip rows change node height (the `jsPlumb.repaint(this)` at `:212` shows
  the pattern).
- **Counts**: repoint the count rendering (`:550-570`) at `.node-counts-strip` — a flex row of the
  existing `.node-count-<type>` divs, keeping `clickable-count`, `node-counts-selected` and the
  legend-colour classes so toggling/legend code is untouched. Delete the floating
  `.count-overlay`/`span.node-counts` positioning CSS (`:517-535`).
- **Factories**: `NODE_FACTORIES` shrinks to `{VennNode: createVennNode}` — the venn widget stays,
  wrapped between the name and counts strip inside a Design A frame. SampleNode/TrioNode/QuadNode
  become default cards (badge carries the pedigree glyph); delete the D3 drawing path in
  `samplenode.js` (`maleSVG`/`femaleSVG`/`addDeceasedStroke`/`sampleNodeUpdateState`/`sourceBadge`)
  and the equivalent trio/quad drawing in `pedigree_node.js`. `PedigreeNode` was already a default
  node; it just gets the new 3-generation badge symbol.
- **Category colour**: strip + badge background from `is_source` — serialize
  `node.get_node_classification()` in `attributes` and colour via
  `.design-a-node.source` / `.design-a-node.filter` CSS (green `#0e7a4a` / blue `#3455a4` as in the
  mock). Selected state stays the existing `.ui-selected` border; analysis-variable orange and
  output-node green keep working through `node-color-overlay` (a corner-chip refinement like the
  mock's is optional polish, not this pass).

### 5. Cleanup

- Remove the replaced CSS blocks from `analysis_nodes.css` (ui-icon rules, per-class
  window/overlay rules, `.grid-link-icon.GeneListNode` → FA `fa-list-check`-based rule).
- Delete images under `images/node_icons/` **only after grepping each filename repo-wide** — several
  are shared: `trio.png`/`trio_*_affected.png`/`quad*.png` are used by `trio_wizard.html:48`,
  `quad_wizard.html`, cohort wizard tags and `view_vcf.html`; `tags_colored.png`,
  `database-1-1.png`, `diagnose_icon.png` live outside `node_icons/` and are used elsewhere. Those
  all stay; converting the wizard pages to the sprite is a follow-up, not this change.
- `./scripts/linting/run_pylint.sh` on touched Python.

## Tests

- New `analysis/tests/test_node_display.py`: iterate `get_node_types_hash()` (instantiating each
  class as `get_nodes_by_classification()` already does) and assert each returns a usable
  `NodeIcon` (exactly one of `fa`/`symbol` set), a non-empty `get_node_class_label_short()`, and
  that `get_node_chips()` returns a list — this is the guard that a future node class can't ship
  without card metadata. For sprite-based icons, assert the symbol id exists in the sprite template.
- SampleNode chip test on existing fake-data fixtures: extraction + N VCFs → expected chip list.
- Existing `analysis/tests/test_urls.py` URL tests cover the serialization path still returning 200.

## Suggested phases (each independently landable)

1. **Server metadata** — `node_display.py`, the three methods on all node classes, rendering-dict
   fields, tests. No visible change.
2. **Sprite + dropdown** — sprite template, `NODE_TYPES` json_script, `_renderItem` rewrite, delete
   ui-icon CSS. Canvas still old-style.
3. **Card rebuild** — new `createDefaultNode`/`updateNodeFromData`/counts strip/CSS, factory
   removals, `revalidate`. The big visible change.
4. **Cleanup + manual QA** — dead CSS/JS/image removal; walk the checklist below on a dev analysis.

Manual QA checklist: create every node type from the dropdown; connect/disconnect edges (top/bottom
endpoints land sensibly with auto heights); counts render and the count-type toggle + legend still
work; select/multi-select highlight; drag; node with long name (clamp + hover); SampleNode with
patient (sex shapes, deceased), with extraction, with 2 VCFs; DamageNode with splicing (chip);
GeneListNode lock; Venn clicking segments; analysis template variables (orange) and output nodes
(green); dual-screen mode; an old saved analysis renders without a hard refresh
(`appearance_version` bump on deploy is *not* automatic — old browsers cache `analysis_nodes.js`,
so verify static file hashing serves the new JS/CSS pair together).
