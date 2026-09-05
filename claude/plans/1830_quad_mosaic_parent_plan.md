# #1830 - Quad node "Dominant (mosaic parent)" inheritance mode

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-05

PR #1833 added `Dominant (mosaic parent)` to the Trio and Duo nodes. This plan brings the same mode
to the Quad node, and fixes the Quad node's stats-cache handler while we're there.

## Why Quad wants this mode

The textbook presentation of parental mosaicism is two affected siblings sharing a variant that the
germline caller wrote off as `0/0` in both parents. No current Quad mode returns it:

- Denovo requires the sibling to be NO_VARIANT, so an affected sibling carrying the variant is excluded.
- Dominant errors out unless a parent is marked affected, and with unaffected parents requires them
  NO_VARIANT, so a parent called HET at a low VAF is excluded too.
- Any Affected returns it, but with no parent constraint at all.

The sibling gives Quad a lever the Trio lacks: an affected sibling must carry the variant, an
unaffected sibling must lack it. That is exactly QuadDominant's existing affected-status rule.

## Data

### `QuadNode` - two new fields, `'M'` added to the inheritance choices

```python
class QuadNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    quad = models.ForeignKey(Quad, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=QuadInheritance.choices,
                                   default=QuadInheritance.RECESSIVE)
    require_parent_zygosity = models.BooleanField(default=True)
    require_sibling_zygosity = models.BooleanField(default=True)
    # Mosaic parent mode only - the low VAF band a mosaic parent's alt reads have to fall in (#1830)
    mosaic_max_af = models.FloatField(default=0.35)
    mosaic_min_alt_reads = models.IntegerField(default=2)
```

### `QuadInheritance` enum (`analysis/models/enums.py`)

```python
class QuadInheritance(models.TextChoices):
    RECESSIVE = 'R', 'Recessive'
    ALL_RECESSIVE = 'A', 'All Recessive (AR + XLR)'
    COMPOUND_HET = 'C', 'C. Het'
    DOMINANT = 'D', 'Dominant'
    MOSAIC_PARENT = 'M', 'Dominant (mosaic parent)'
    DENOVO = 'N', 'Denovo'
    XLINKED_RECESSIVE = 'X', 'X-Linked Recessive'
    ANY_AFFECTED = 'Y', 'Any Affected (variant in ≥1 affected)'
```

`'M'` is the same letter Trio and Duo use, and sits in the same position after Dominant.

### Migration `analysis/migrations/0135_quad_mosaic_parent.py`

Depends on `0134_mosaic_parent_inheritance`. Two `AddField`s on `quadnode` (same defaults as Trio) and
an `AlterField` on `quadnode.inheritance` for the new choices. Generate with `makemigrations`, then
rename to that filename.

## The filter - `QuadMosaicParent`

Lives in `analysis/models/nodes/sources/quad_node.py`, subclassing `AbstractQuadInheritance` and
mirroring Trio's `MosaicParent` (`trio_node.py:76`). The shared helpers from
`analysis/models/nodes/family_inheritance.py` do the AD/AF work: `mosaic_evidence_q`,
`mosaic_absent_q`, `mosaic_evidence_description`, `MOSAIC_PARENT_WARNINGS`, `MOSAIC_ZYGOSITIES`.

Each side is the Trio's zygosity tuple with a sibling term appended, then run through `_get_zyg_q`
so `require_parent_zygosity` / `require_sibling_zygosity` apply as they do for every other Quad mode:

| Side           | mother             | father             | proband     | sibling                                            |
|----------------|--------------------|--------------------|-------------|----------------------------------------------------|
| mother mosaic  | MOSAIC_ZYGOSITIES  | NO_VARIANT         | HAS_VARIANT | `UNAFFECTED_AND_AFFECTED_ZYGOSITIES[sibling_affected]` |
| father mosaic  | NO_VARIANT         | MOSAIC_ZYGOSITIES  | HAS_VARIANT | same                                               |

Each side is additionally ANDed with `mosaic_evidence_q(cgc, mosaic_parent_sample, mosaic_max_af,
mosaic_min_alt_reads)` and `mosaic_absent_q(cgc, other_parent_sample, mosaic_min_alt_reads)`, then the
two sides are ORed. The sibling is constitutional (inherited it fully or not at all), so its genotype
call is the whole constraint - the AD/AF evidence applies to the parents only.

A `_sibling_zyg()` helper returning `UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.sibling_affected)]`
keeps the two sides in lockstep and is what `get_zygosity_table_data` reads.

`get_method()` follows the Trio wording with the sibling added, e.g.
`Proband: HET, HOM_ALT, and either parent (HOM_REF, HET, MISSING) with ≥2 alt reads at AF ≤ 0.35 while the other (HOM_REF, MISSING) has <2 alt reads; Sibling: ...`.
`get_other_filters_description()` returns `mosaic_evidence_description(...)` as Trio does.

## QuadNode wiring

- `INHERITANCE_CLASSES[QuadInheritance.MOSAIC_PARENT] = QuadMosaicParent`.
- `get_warnings()` - new override in the Trio's shape: `super().get_warnings()` plus
  `MOSAIC_PARENT_WARNINGS` when `self.quad and self.inheritance == QuadInheritance.MOSAIC_PARENT`.
  `FamilyInheritanceNodeMixin._load` already turns warnings into the yellow shadow.
- `get_quad_inheritance_errors` - no entry for the new mode. Mosaic parent has no affected-parent
  requirement (the mosaic parent is typically unaffected).
- `get_help_text()` - append the Trio's sentence about the mode.
- `get_zygosity_table_data()` - a `klass is QuadMosaicParent` branch (before the CompHet fallback)
  emitting `mother`/`father` as `fmt(MOSAIC_ZYGOSITIES)`, `proband` as `fmt(HAS_VARIANT)`, and
  `sibling_affected` / `sibling_unaffected` from `_sibling_zyg()` with the stub's `sibling_affected`
  set each way. The existing `description` loop then adds the AD/AF text to `other_filters_*`.

## Stats cache - `QuadNode` is reading the raw cohort row

`get_handler_for_node` in `analysis/models/nodes/stats_cache.py` has no `QuadNode` entry, so it
falls to `NoFilterHandler`, and `QuadNode` does not override `_has_filters_that_affect_label_counts`.
Every Quad inheritance mode therefore reports the whole cohort's count in the node's TOTAL/OMIM/etc
labels - the bug #1833 fixed for Duo. Same fix:

- `class QuadInheritanceHandler(FilterKeyHandler)` returning `UNCACHEABLE` from `filter_key_for_node`
  (nothing is precomputed for quads - the buckets are built off a cohort's Trio).
- Register `QuadNode: QuadInheritanceHandler()` in the handler map, importing `QuadNode` alongside the
  others there.
- A test in `analysis/tests/test_cohort_stats_cache.py` asserting `get_handler_for_node(QuadNode())`
  yields `UNCACHEABLE`, next to the Duo one.

## Form and editor

- `QuadNodeForm` (`analysis/forms/forms_nodes.py:1146`) - add `"mosaic_max_af": WIDGET_UNIT_INTERVAL`
  and `"mosaic_min_alt_reads": WIDGET_INTEGER_MIN_1` to `widgets`, as the Trio/Duo forms do.
- `analysis/templates/analysis/node_editors/quadnode_editor.html` - copy the Trio editor's pieces:
  `updateMosaicFields()` toggling `#mosaic-fields` on `$("#id_inheritance").val() === 'M'`, called on
  inheritance change and at the end of `docreadyjs`; the `genotype-row#mosaic-fields` div with the
  two inputs, placed directly above the existing AD/DP/GQ/PL row.
- The zygosity table already has a Sibling row; the `sibling_affected` / `sibling_unaffected` keys
  flow through the existing `memberZyg('sibling', siblingAffected, requireSiblingZyg)`.

## Tests - `analysis/tests/test_quad_node.py`

The shared-mode tests come from `InheritanceNodeTestsMixin`; the mosaic ones go in the Quad class,
built with `make_cohort_genotype`'s per-sample AD/AF (packed order `[proband, mother, father,
sibling]`). Both quads in `setUpTestData` are available - `cls.quad` (sibling unaffected, `cls.cgc`)
and `cls.quad_aff` (sibling affected, `cls.cgc_aff`).

Port the Trio's mosaic tests (matches HOM_REF-called mother with 3 alt reads at 6%; matches low-VAF
HET father; excludes full HET parent; excludes below-threshold reads; excludes reads in both parents;
excludes variants the proband lacks; missing AF value and NULL AF column both match; both thresholds
configurable; source node; always warns; needs no affected parent). Then the sibling-specific ones:

- Unaffected-sibling quad: proband HET + mosaic mother + sibling HET is excluded; sibling HOM_REF
  matches.
- Affected-sibling quad: proband HET + mosaic mother + sibling HET matches (the two-affected-sibs
  case that motivates the mode); sibling HOM_REF is excluded.
- `require_sibling_zygosity=False` lets a sibling UNKNOWN through; `True` excludes it.
- Zygosity table: the mode's entry has `sibling_affected`, `sibling_unaffected`, and an
  `other_filters_mother` mentioning alt reads.

Audit the tests at the end per CLAUDE.md - keep the ones covering our branches (the two sides, the
sibling rule, the thresholds, the AF-absent shapes), drop any that only restate a field declaration.

## Changelog

`variantgrid/templates/default_templates/changelog.html` - extend the two existing #1830 lines
(the prose one near line 44 and the issue list one near line 268) from "Trio/Duo nodes" to
"Trio/Duo/Quad nodes".

## Out of scope

`QuadCompHet`'s sibling TODO and any precomputed stats buckets for quads.
