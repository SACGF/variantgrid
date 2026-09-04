# Duo analysis - proband + one parent

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-04

Design for [#1829](https://github.com/SACGF/variantgrid/issues/1829). A medical scientist wants an easy
way to analyse a **duo**: a proband plus one sequenced parent (father unavailable, prenatal case entered
under the mother's record, deceased parent, cost). Today the only route is a 2-sample Cohort node with
hand-set zygosities - no comp het, no inheritance labels, no pedigree figure.

`Duo` becomes a third family model beside `Trio` and `Quad`, with its own wizard and `DuoNode`, built on
the family backend Quad already shares (`analysis/models/nodes/family_inheritance.py`,
`FamilyGroupMixin` in `snpdb/models/models_cohort.py`). Quad is the template throughout: every file
listed in §5 has a Quad counterpart to copy from.

Two-sibling duos (no parents) stay on the Cohort node - shared-variant search already works there.

---

## 1. Data

```python
class DuoRelationship(models.TextChoices):   # snpdb/models/models_enums.py
    MOTHER = 'M', 'Mother'
    FATHER = 'F', 'Father'


class Duo(FamilyGroupMixin, GuardianPermissionsAutoInitialSaveMixin, SvgSymbolPreviewIconMixin,
          PreviewModelMixin, SortByPKMixin, TimeStampedModel):
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    proband = models.ForeignKey(CohortSample, related_name='duo_proband', on_delete=CASCADE)
    parent = models.ForeignKey(CohortSample, related_name='duo_parent', on_delete=CASCADE)
    relationship = models.CharField(max_length=1, choices=DuoRelationship.choices)
    parent_affected = models.BooleanField(default=False)
    # Set in the duo wizard when the scientist resolves patient.sex vs sample.detected_sex
    proband_sex = models.CharField(max_length=1, choices=Sex.choices, null=True, blank=True)

    preview_icon_symbol = "node-icon-duo"
    pedigree_icon_members = ("parent",)
```

`parent` + `relationship` rather than nullable `mother`/`father` columns: one FK that is always set,
and the inheritance code needs the role anyway (X-linked recessive is only meaningful when the parent
is the mother; comp het is described as "one from mother, one not from mother").

```python
class DuoInheritance(models.TextChoices):   # analysis/models/enums.py
    RECESSIVE = 'R', 'Recessive'
    ALL_RECESSIVE = 'A', 'All Recessive (AR + XLR)'
    COMPOUND_HET = 'C', 'C. Het (half phased)'
    DOMINANT = 'D', 'Dominant'
    ABSENT_IN_PARENT = 'N', 'Absent in parent'
    XLINKED_RECESSIVE = 'X', 'X-Linked Recessive'
    ANY_AFFECTED = 'Y', 'Any Affected (variant in ≥1 affected)'


class DuoSample(models.TextChoices):        # analysis/models/enums.py - wizard roles
    PARENT = 'A', 'Parent'
    PROBAND = 'P', 'Proband'


class DuoNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    duo = models.ForeignKey(Duo, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=DuoInheritance.choices,
                                   default=DuoInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)   # parent only - proband always required (#947)
```

Same letter codes as `TrioInheritance` so the analysis-template machinery and anyone reading the DB sees
the same thing mean the same thing. `ABSENT_IN_PARENT` keeps Trio's `'N'` (de novo) slot because it is
the one-parent version of that filter.

Migrations: `snpdb` (Duo), `analysis` (DuoNode). New tables only - nothing existing changes.

## 2. Inheritance modes

Each mode is the Trio mode with the missing parent's constraint dropped. Zygosity sets are the ones
already on `AbstractFamilyInheritance` (`NO_VARIANT`, `HAS_VARIANT`,
`UNAFFECTED_AND_AFFECTED_ZYGOSITIES`).

| Mode | Parent | Proband | Extra |
|---|---|---|---|
| Recessive | HET | HOM_ALT | |
| X-linked recessive | HET | HOM_ALT | chrX only; error when `relationship` is FATHER (see §2.1) |
| All recessive | AR ∪ XLR | | XLR branch chrX only; collapses to AR alone when the parent is the father |
| Dominant | `UNAFFECTED_AND_AFFECTED_ZYGOSITIES[parent_affected]` | HAS_VARIANT | |
| Absent in parent | NO_VARIANT | HAS_VARIANT | |
| C. Het (half phased) | `_mum_but_not_dad` = (parent HET, proband HET); `_dad_but_not_mum` = (parent NO_VARIANT, proband HET) | | two-hit gene logic from `AbstractCompHetInheritance` unchanged |
| Any affected | present in ≥1 of the affected members (proband always counts) | | same shape as `TrioAnyAffected` |

`_get_zyg_q` builds `[(duo.parent.sample, parent_zyg, require_zygosity), (duo.proband.sample, proband_zyg, True)]`
through `_build_family_zyg_q`. `get_method()` strings name the parent by its relationship
(`Mother:{...} Proband:{...}`) so the node's method summary reads naturally.

`get_other_filters_description()` for comp het: `"≥2 hits in same gene, one from the {mother}, one not"`.

### 2.1 Errors and warnings

Through `FamilyInheritanceNodeMixin` exactly as Trio/Quad, so the ignore-as-warning button works:

- Dominant: `_dominant_requires_affected_parent_error(parent_affected, False)` - message reads
  "Dominant inheritance requires an affected parent" as today.
- X-linked recessive: `_xlinked_recessive_errors(proband_sample, effective_proband_sex, mother_affected)`
  when the parent is the mother; when the parent is the father a new error
  "X-linked recessive needs the mother - the father's chrX is not passed to a son".
- Absent in parent: a **warning** rather than an error, always present, so the scientist sees it on the
  node: "One parent only - de novo cannot be confirmed; variant may be inherited from the missing
  {father}". `FamilyInheritanceNodeMixin._load` already turns warnings into the yellow shadow.

## 3. Wizard

`duo_wizard(request, cohort_id, sample1_id, sample2_id)` in `analysis/views/views_wizard.py`,
`UserDuoWizardForm` beside `UserQuadWizardForm`, template `analysis/duo_wizard.html`. Roles are
Proband / Parent with a Mother-or-Father radio and an affected tick on the parent; proband sex resolution
(`proband_sex_mismatch.js`) and the pedigree figure (`family_wizard_roles.js`) work as for Trio - the
figure gets the mother/father symbol from the radio.

Entry points mirror `cohort_quad_wizard_tag.html`: a `cohort_duo_wizard_tag.html` shown on the cohort
details tab and VCF page when `sample_count >= 2`, and the Cohort node editor's family shortcut if Trio
has one there.

The `Duo` name defaults to `"{parent}/{proband} from {cohort}"`; `get_or_create` on
(cohort, proband, parent, relationship, parent_affected) like the Quad wizard, updating `proband_sex` on
an existing row.

## 4. UI

- **Node**: `DuoNodeForm` (`analysis/forms/forms_nodes.py`, `duo_autocomplete`), `DuoNodeSerializer`,
  `duonode_editor.html` copied from `quadnode_editor.html` with two member rows. The zygosity table is
  driven by `get_zygosity_display_data()` on the node with `members = ['parent', 'proband']`; the parent
  row is labelled by relationship. `initially_show_zygosity_table` user setting covers it (rename the
  help text to "Trio/Quad/Duo").
- **Icon**: `node-icon-duo` symbol in `uicore/templates/uicore/tags/svg_icon_sprite.html` - square or
  circle above, proband below; `parent-affected` CSS fills the parent shape (`global.scss`, and the same
  minimal edit in `global.css`). One symbol with both parent shapes and a `duo-mother`/`duo-father` class
  choosing which is drawn keeps the sprite static.
- **Pages**: `duos` listing (`DuosListColumns` on `FamilyGroupListColumns` with
  `FAMILY_MEMBERS = [("parent", "Parent", True), ("proband", "Proband", False)]` plus a relationship
  column), `view_duo` with `duo_table` tag, sharing tab, related analyses
  (`get_related_analysis_details_for_duo`), `duo_search` signal, `DuoAutocompleteView`, `DuoSerializer` +
  `api/duo/<pk>`.
- **Analysis templates**: add `"duo"` to `ANALYSIS_TEMPLATE_VCF_SOURCE_FIELDS` and the source-field
  loop in `models_analysis.py` (`'snpdb.Duo': {'snpdb.Sample'}` in the type map), so a template built on
  a Duo can be instantiated against another Duo.
- **Invalidation**: `handle_duo_pre_delete` in `analysis/signals/source_data_invalidation.py`, connected
  in `analysis/apps.py`.

## 5. Files (Quad counterpart → Duo)

| Area | Quad file | Duo |
|---|---|---|
| model | `snpdb/models/models_cohort.py` (`Quad`) | `Duo` after `Quad`; `DuoRelationship` in `models_enums.py` |
| enums | `analysis/models/enums.py` (`QuadInheritance`, `QuadSample`) | `DuoInheritance`, `DuoSample` |
| node | `analysis/models/nodes/sources/quad_node.py` | `duo_node.py`, exported from `sources/__init__.py` |
| node form/serializer | `forms_nodes.py` `QuadNodeForm`, `serializers.py` `QuadNodeSerializer` | same |
| node editor | `node_editors/quadnode_editor.html`, `node_views.py` (`quad_loadable`) | `duonode_editor.html`, `duo_loadable` |
| wizard | `views_wizard.py` `quad_wizard`, `analysis/forms` `UserQuadWizardForm`, `quad_wizard.html`, `cohort_quad_wizard_tag.html` + templatetag | `duo_*` |
| listing/view | `views_cohort.py` `quads`/`view_quad`, `quads.html`, `view_quad.html`, `quad_table.html` + `model_tags.quad_table`, `grids.py` `QuadsListColumns` | `duo*` |
| urls | `snpdb/urls.py` (`quads`, `view_quad`, `quad_datatable`, `quad_autocomplete`, `api/quad`), `analysis/urls.py` (`quad_wizard`) | `duo*` |
| search/autocomplete/rest | `snpdb/signals/quad_search.py`, `views_autocomplete.QuadAutocompleteView`, `views_rest.QuadView`, `snpdb/serializers.QuadSerializer` | `Duo*` |
| related analyses | `analysis/related_analyses.py`, `related_analyses_tags.py` | add duo |
| invalidation | `source_data_invalidation.py`, `analysis/apps.py` | `handle_duo_pre_delete` |
| templates | `models_analysis.py` source fields | add `"duo"` |
| icon/css | `svg_icon_sprite.html`, `global.scss`/`.css` | `node-icon-duo`, `parent-affected` |
| tests | `analysis/tests/test_quad_node.py`, `inheritance_node_mixin.py`, `snpdb/tests/utils/fake_cohort_data.py` `create_fake_quad`, `test_urls.py` (both apps), `variantopedia/tests/test_search.py` | `test_duo_node.py`, `create_fake_duo`, url + search rows |
| changelog | `changelog.html` | entry |

Run `scripts/vg map` after the models and urls land.

## 6. Tests

`analysis/tests/test_duo_node.py` on `inheritance_node_mixin.py`, following `test_quad_node.py`:

- one test per mode asserting the parent/proband zygosity sets, both relationships where they differ
  (XLR and All Recessive with mother vs father).
- comp het: half-phased pairing returns the gene with one hit in the parent and one not, and excludes a
  gene where both hits are in the parent.
- errors: dominant with unaffected parent; XLR with a father parent; XLR with a female proband; absent
  in parent always carries the warning.
- `require_zygosity=False` widens the parent set only.

`snpdb/tests/test_urls.py` and `analysis/tests/test_urls.py` gain the duo rows; `test_search.py` gains
a duo name search.

## 7. Follow-ups

**Mosaic parent - [#1830](https://github.com/SACGF/variantgrid/issues/1830).** A "Dominant (mosaic
parent)" mode (parent zygosity unconstrained, parent alt AD / AF in a low range, proband HAS_VARIANT)
is planned for Trio and belongs on Duo too. With one parent the "other parent clean" clause drops, so
the Duo version is the simpler one; leave `DuoInheritance` room for it and keep `require_zygosity`
per-parent so the mode can widen the parent's zygosity set without touching the proband.

**Members-driven family backend.**
Duo is the third copy of `_get_zyg_q` / `get_method` / `get_zygosity_display_data` / the editor template,
each hard-coding member names. All three are "list of (member, sample, zyg_set, require)" underneath; a
members-driven `AbstractFamilyInheritance` and editor would make Duo (and Quad) mostly configuration.
Separate issue once Duo is in.
