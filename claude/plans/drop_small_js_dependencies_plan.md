# Plan: Drop lodash and d3 — small in-house replacements for whole libraries

Companion to issue #1790 (re-download cut-down jQuery UI / Plotly builds), which needs no
application code change. This plan covers the libraries where a small, self-contained fix
retires the whole dependency.

## Decision

Two vendored libraries are each doing a job that modern browser APIs — or a handful of
lines of our own code — now cover completely. Retire them:

| Library | Size | Loaded from | Call sites |
| --- | --- | --- | --- |
| `lodash/4.18.1/lodash.min.js` | 73KB | `base.html` (every page) | 13 |
| `d3.v2.js` | 240KB | 4 templates | 2 files, `select`/`append`/`attr`/`style`/`on` only |

A third, `highlight/highlight.pack.js` (13KB), stays but stops loading on every analysis
page — see §3.

Total: ~313KB of vendored JS retired, 73KB of it off the critical path of every page.

`jquery.blockUI.js` (20KB) was in this plan and has already been deleted — its only caller
was `upload/templates/upload/upload_forms.html`, an orphaned template no view rendered. Both
are gone, along with the unused script tag in `analysis_editor_and_grid.html`. If a use
resurfaces, `LoadingOverlay` is loaded globally and is the established idiom
(`global.js:369`, `676`, `1254`; `flags.js:425`).

Each section below is independently shippable and independently revertable. Do them in
order — §1 touches the most files, §2 is the most delicate.

`eslint.config.mjs` declares `_`, `d3` and friends as readonly globals for the vendored
`<script>` tags. Each section must remove its own entry from that list, so a stray
reference to a retired library becomes a lint error rather than a runtime one.

---

## §1 — lodash

73KB loaded on every page for 13 calls, every one of which has a native equivalent in the
ES2022 baseline the project already targets (`eslint.config.mjs` sets `ecmaVersion: 2022`).

### Call sites and replacements

| File:line | Call | Replacement |
| --- | --- | --- |
| `global.js:1319` | `_.isNumber(jsonObj)` | `typeof jsonObj === "number"` |
| `global.js:1321` | `_.isBoolean(jsonObj)` | `typeof jsonObj === "boolean"` |
| `global.js:1323` | `_.isString(jsonObj)` | `typeof jsonObj === "string"` |
| `global.js:1335` | `_.isArray(jsonObj)` | `Array.isArray(jsonObj)` |
| `variant_details.html:186` | `_.isNumber(value)` | `typeof value === "number"` |
| `variant_details.html:223` | `_.isString(value)` | `typeof value === "string"` |
| `variant_details.html:228` | `_.isBoolean(value)` | `typeof value === "boolean"` |
| `vc_keys.js:340` | `_.some(ekey.options, fn)` | `ekey.options.some(fn)` |
| `vc_keys.js:343` | `_.cloneDeep(ekey.options)` | `structuredClone(ekey.options)` |
| `datatable_definition.js:105` | `_.random(0, 100000)` | `Math.floor(Math.random() * 100001)` |
| `vc_form.js:1771` | `_.flatten(_.partition(src, pred))` | see below |
| `condition_matching.html:104` | `_.maxBy(messages, fn)` | see below |
| `condition_matchings.html:34` | `_.debounce(fn, 250, {maxWait: 1000})` | see §1.1 |

Both type-guard groups are `else if` chains over values that arrive from `JSON.parse`, so
they are always primitives — the boxed-object cases lodash also matches (`new Number(1)`)
cannot occur. `typeof NaN === "number"` matches `_.isNumber(NaN)`, so the numeric branch is
unchanged for NaN.

`structuredClone` matches `_.cloneDeep` for the JSON-shaped `EvidenceKey.options` arrays it
is used on.

`_.random(0, 100000)` is inclusive at both ends, hence the `* 100001`. The value only seeds
a DOM id (`"tid" + n`), so the exact distribution is immaterial — keep it inclusive anyway
so the change is a pure substitution.

**`vc_form.js:1771`** — `_.partition` returns `[matching, notMatching]` and `_.flatten`
concatenates one level, i.e. "stable sort, non-excluded first". Directly:

```js
optionSources = [
    ...optionSources.filter(o => !o.exclude_namespace),
    ...optionSources.filter(o => o.exclude_namespace),
];
```

**`condition_matching.html:104`** — `_.maxBy` returns the *first* element attaining the
maximum, and `undefined` for an empty array. `reduce` with a strict `>` preserves both:

```js
let maxMessage = suggestion.messages.reduce((best, message) => {
    const rank = severity[message.severity.toLowerCase()[0]] || 0;
    return (best === null || rank > best.rank) ? {rank, message} : best;
}, null);
return maxMessage ? maxMessage.message.severity : null;
```

### §1.1 — `_.debounce`

The one call needing real code. `condition_matchings.html` uses two lodash features beyond
plain debouncing, and both are load-bearing:

- `{maxWait: 1000}` — guarantees the datatable reloads at least once a second while the
  user keeps typing, rather than only after they stop.
- `.flush()` — the Filter button (`condition_matchings.html:110`) calls
  `debouncedApplyFilters.flush()` to apply a pending edit immediately.

Add a `debounce` helper to `global.js` alongside the other shared utilities:

```js
// Trailing-edge debounce with lodash-compatible maxWait and flush()
function debounce(func, wait, {maxWait} = {}) {
    let timer = null;
    let maxTimer = null;
    let pendingArgs = null;
    let pendingThis = null;

    function invoke() {
        clearTimeout(timer);
        clearTimeout(maxTimer);
        timer = maxTimer = null;
        if (pendingArgs) {
            const [args, self] = [pendingArgs, pendingThis];
            pendingArgs = pendingThis = null;
            func.apply(self, args);
        }
    }

    function debounced(...args) {
        pendingArgs = args;
        pendingThis = this;
        clearTimeout(timer);
        timer = setTimeout(invoke, wait);
        if (maxWait && !maxTimer) {
            maxTimer = setTimeout(invoke, maxWait);
        }
    }
    debounced.flush = invoke;
    debounced.cancel = function() {
        clearTimeout(timer);
        clearTimeout(maxTimer);
        timer = maxTimer = pendingArgs = pendingThis = null;
    };
    return debounced;
}
```

`invoke` is a no-op when nothing is pending, so a `flush()` on an idle debouncer does not
fire a spurious reload — matching lodash.

### §1.2 — Finish up

- `global.js:895` — delete the stale `// deprecated use _.pad(num + "", size, '0');` comment
  inside `zero_pad()`. It is the last mention of lodash and refers to a migration that never
  happened; `zero_pad` itself stays.
- Remove the `<script src=".../lodash.min.js">` tag and its `<!-- lodash for our own usage -->`
  comment from `uicore/templates/uicore/page/base.html` (~line 97-98).
- `git rm -r variantgrid/static_files/default_static/js/lib/lodash`
- Remove `_: "readonly"` from `eslint.config.mjs`, then run `npx eslint` to confirm nothing
  else reaches for it.

### §1.3 — Verifying

Every replacement is a pure expression swap with no I/O, so the check is visual rather than
automated — exercise each page once:

- Variant details: the transcript info table renders numbers with threshold colouring,
  strings with code-word colouring, booleans as high/low badges (`variant_details.html`).
- Any page rendering a JSON blob through `formatJson` (`global.js`) — numbers, booleans,
  strings, arrays and nulls all still get their `js-*` classes.
- Classification form: an evidence key with namespaced options still shows namespace
  filtering, and select options put excluded ones last (`vc_form.js`, `vc_keys.js`).
- Condition matchings: typing in the text filter reloads the table on a 250ms trailing edge
  and at least once per second while typing; the Filter button applies immediately.
- Condition matching: a suggestion with several messages shows the most severe one's colour.

---

## §2 — d3 v2

`d3.v2.js` is 240KB, dated 2012, and loaded by `analysis.html`, `analysis_includes.html`,
`patient_samplenode_head.html` and `upload.html`. Two of our files use it, and between them
they touch exactly five d3 APIs:

```
d3.select(nodeOrSelector)   .append("svg:tag")   .attr(name[, value])   .style(name, value)   .on("click", fn)
```

No data joins, no `selectAll`, no scales, axes, layouts or transitions. d3 is there purely
because jQuery cannot create namespaced SVG elements — `$("<circle>")` produces an HTML
element that never renders inside an `<svg>`.

### §2.1 — Delete `venn3` first

`venn_intersect.js:125-214` defines `venn3()`, which is called from nowhere in the repo — its
click handlers are literally `alert(2)` through `alert(7)`, i.e. unfinished scaffolding. It is
90 of the file's 215 lines and roughly half its d3 usage. Delete it before doing any
conversion work, so none of it needs converting.

Only `venn2`, `venn_select` and `vennAddToggleCallbacks` are live, all from
`analysis/templates/analysis/node_editors/vennnode_editor.html` and `analysis_nodes.js:105`.

### §2.2 — The replacement

Add a small chainable SVG builder to `global.js`. It needs to reproduce four d3 behaviours
the calling code depends on:

- `append` returns a selection wrapping the **new child**, so chains keep descending.
- `attr`/`style` return the **same** selection, so chains stay flat.
- `attr(name)` with one argument **reads** — `venn_intersect.js` reads back `toggled` and
  `venn_bit`, and gets a string (`"true"` / `"1"`), which its `==` comparisons rely on.
- Inside an `on` handler, `this` is the **DOM node** — `vennAddToggleCallbacks` calls
  `toggleSelect(d3.select(this))` from a click handler.

```js
// Minimal chainable SVG element builder - replaces the d3 v2 subset we used.
// jQuery can't do this: SVG elements need createElementNS to render.
const SVG_NS = "http://www.w3.org/2000/svg";

class SvgSelection {
    constructor(node) {
        this.node = node;
    }

    // Accepts "svg:rect" or "rect" - the d3 v2 namespace prefix is kept so call sites read the same
    append(tagName) {
        const child = document.createElementNS(SVG_NS, tagName.replace(/^svg:/, ""));
        this.node.appendChild(child);
        return new SvgSelection(child);
    }

    attr(name, value) {
        if (value === undefined) {
            return this.node.getAttribute(name);
        }
        this.node.setAttribute(name, value);
        return this;
    }

    style(name, value) {
        this.node.style.setProperty(name, value);
        return this;
    }

    on(eventName, handler) {
        this.node.addEventListener(eventName, handler.bind(this.node));
        return this;
    }
}

function svgSelect(target) {
    const node = typeof target === "string" ? document.querySelector(target) : (target.jquery ? target[0] : target);
    return new SvgSelection(node);
}
```

Two notes on fidelity:

- `handler.bind(this.node)` reproduces d3's `this`-is-the-node contract. Without it,
  `toggleSelect(d3.select(this))` in `vennAddToggleCallbacks` receives the `window`.
- `style.setProperty` stringifies numbers, which matters for
  `.style("stroke-width", 2)` in `venn2`.

`d3.select` also accepts a DOM node, a selector string, or (via `analysis_nodes.js:105`
passing `vennHolder[0]`) an already-unwrapped node — `svgSelect` handles all three, plus a
jQuery object for good measure.

### §2.3 — Convert the call sites

Mechanical: `d3.select(` → `svgSelect(` in

- `venn_intersect.js` lines 7, 16, 25, 31, 55 (after `venn3` is gone)
- `samplenode.js:94`

Everything downstream — `femaleSVG`, `maleSVG`, `addDeceasedStroke`, the `addCirc1`/`addCirc2`
/`setRing` helpers in `venn2` — receives a selection and needs no change at all, because the
chainable surface is identical.

Then drop the four `<script src=".../d3.v2.js">` tags:

- `analysis/templates/analysis/analysis.html:13`
- `analysis/templates/analysis/analysis_includes.html:4`
- `patients/templates/patients/patient_samplenode_head.html:2`
- `upload/templates/upload/upload.html:8`

…`git rm variantgrid/static_files/default_static/js/lib/d3.v2.js`, and remove
`d3: "readonly"` from `eslint.config.mjs`.

Check `global.js` is loaded before `samplenode.js`/`venn_intersect.js` on each of those four
pages — it is on anything extending `base.html`, but `analysis_editor_and_grid.html` builds
its own header for `stand_alone` mode and should be confirmed.

### §2.4 — Verifying

The two widgets are small and visual:

- Venn node editor (`vennnode_editor.html`): the two-circle Venn renders with its rings;
  clicking each of the three regions toggles it red/white; the region flags round-trip
  (reopen the editor and the selection is still there — that path goes through
  `venn_select` reading the `toggled` attribute back).
- Analysis canvas: a Venn node's badge renders (`analysis_nodes.js:105`).
- Patient form example node and the upload page: male (square) and female (circle) pedigree
  shapes render, with the deceased diagonal stroke where applicable.

---

## §3 — highlight.js: load it where it is used

Not a removal — a scope reduction. `highlight/highlight.pack.js` (13KB) plus its stylesheet
load from `analysis/templates/analysis/analysis_includes.html:2,5`, i.e. on **every analysis
page**, for a single `hljs.highlightBlock` call in
`analysis/templates/analysis/node_editors/grid_editor_debug_tab.html:26` — a developer-only
tab.

Move the cost to the tab. The debug tab is rendered into the page by ajax
(`analysis/views/views.py:557`), and a `<script src>` inside injected HTML loads
asynchronously — so its inline `$(document).ready(...)` would run before `hljs` exists. Load
it explicitly instead:

```js
$.getScript("{% static 'js/lib/highlight/highlight.pack.js' %}", function() {
    $('pre code.sql').each(function(i, block) {
        hljs.highlightBlock(block);
    });
});
```

`$.getScript` caches by URL within a page, so reopening the tab does not refetch. The
stylesheet link moves to the same template.

Then remove lines 2 and 5 from `analysis_includes.html`.

**Verify:** the node editor's Debug tab still syntax-highlights the SQL blocks, and the tab
works when opened twice.

While in `analysis_includes.html`, note that `hljs.highlightBlock` was deprecated in
highlight.js 10 in favour of `highlightElement` — our pinned build predates that, so leave
the call as-is unless the library is upgraded.

---

## Out of scope

- **moment (58KB) + jquery.timeago (10KB).** Genuinely replaceable by `Intl.DateTimeFormat`
  and `Intl.RelativeTimeFormat`, but it is 28 call sites across `global.js`, `flags.js`,
  `vc_diff.js`, `datatable_definition.js` and `vc_form.js`, each with its own format string
  (`'DD-MMM-YYYY HH:mm'`, `'hh:mm A'`, …) that has to be re-expressed as `Intl` options and
  eyeballed. That is a plan of its own, not a small fix.
- **js.cookie (1KB, 11 call sites)** and **jquery.form / `ajaxForm` (17KB, 38 call sites)** —
  both earn their keep at their size.
- **jQuery UI and Plotly** — issue #1790. Those are download-a-smaller-build tasks with no
  code change, deliberately kept separate from this one.
