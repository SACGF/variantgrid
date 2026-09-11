# #1861 — Duo: proband + sibling

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-11
Status: draft

[#1861](https://github.com/SACGF/variantgrid/issues/1861) asks for the Duo to cover the two kinds of pair a lab
sees: parent + child (parent may or may not be affected) and two siblings (usually both affected). Today the
Duo's second member is always a parent: `snpdb/models/models_cohort.py:Duo` stores it in a `parent` FK with a
`relationship` of Mother or Father, and every `analysis/models/nodes/sources/duo_node.py:DuoNode` mode is a
(parent, proband) zygosity pair. The Quad already filters a sibling by its affected tick
(`analysis/models/nodes/sources/quad_node.py:QuadRecessive`), so a sibling Duo is those sibling rules with the
parents dropped.

Decision: extend the Duo with a third relationship rather than add a sibling-pair model and node. One relationship
value is not worth a second copy of the family-inheritance machinery, and the wizard, editor, listing and
pedigree figure all key off `relationship` already.

## Data

`snpdb/models/models_cohort.py:Duo` - the second member stops being a parent by definition, so the two fields
that name it as one are renamed. Existing rows are all Mother/Father and need no backfill.

```python
class Duo(FamilyGroupMixin, GuardianPermissionsAutoInitialSaveMixin, SvgSymbolPreviewIconMixin, PreviewModelMixin,
          SortByPKMixin, TimeStampedModel):
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    proband = models.ForeignKey(CohortSample, related_name='duo_proband', on_delete=CASCADE)
    relative = models.ForeignKey(CohortSample, related_name='duo_relative', on_delete=CASCADE)  # was parent
    relationship = models.CharField(max_length=1, choices=DuoRelationship.choices)
    relative_affected = models.BooleanField(default=False)  # was parent_affected
    proband_sex = models.CharField(max_length=1, choices=Sex.choices, null=True, blank=True)
```

`snpdb/models/models_enums.py:DuoRelationship` and its wizard twin `analysis/models/enums.py:DuoSample` (the
values are shared so the role picked in the wizard is stored as the relationship):

```python
class DuoRelationship(models.TextChoices):
    MOTHER = 'M', 'Mother'
    FATHER = 'F', 'Father'
    SIBLING = 'S', 'Sibling'

class DuoSample(models.TextChoices):
    MOTHER = 'M', 'Mother'
    FATHER = 'F', 'Father'
    SIBLING = 'S', 'Sibling'
    PROBAND = 'P', 'Proband'
```

Migration snpdb/migrations/0262_duo_sibling.py (new): `RenameField` parent → relative, `RenameField`
parent_affected → relative_affected, `AlterField` relationship (new choice). `DuoNode` itself does not change
shape - `require_zygosity` keeps applying to the relative, whoever they are.

## Zygosity rules per mode

`R` = relative, `P` = proband. Parent columns are today's behaviour, unchanged. "not HOM_ALT" is
`{HET, HOM_REF, MISSING}`, the Quad's unaffected-sibling set.

| Mode | Parent (as now) | Sibling, affected | Sibling, unaffected |
|---|---|---|---|
| Recessive | R HET, P HOM_ALT | R HOM_ALT, P HOM_ALT | R not HOM_ALT, P HOM_ALT |
| Dominant | R has/lacks by affected, P has | R has, P has | R lacks, P has |
| X-linked recessive | R HET, P HOM_ALT, chrX | R HOM_ALT, P HOM_ALT, chrX | R not HOM_ALT, P HOM_ALT, chrX |
| All recessive | AR, plus XLR branch when mother | AR or XLR (both branches) | AR or XLR (both branches) |
| C. Het | one hit R HET, one hit R lacks | R HET on every hit | R unconstrained on every hit |
| Any affected | OR of affected members | OR of both | P alone |
| Absent in parent | R lacks, P has | unavailable | unavailable |
| Dominant (mosaic parent) | R low-VAF, P has | unavailable | unavailable |

Notes on the choices:
- Recessive with an affected sibling asks for HOM_ALT, the same as `QuadXLinkedRecessive` does. `QuadRecessive`
  takes HAS_VARIANT for an affected sibling, which is looser than its own XLR mode; leave the Quad as it is,
  it's out of scope here.
- Dominant with an unaffected sibling is a valid discordant-pair filter (proband has it, sibling doesn't), so
  the "Dominant inheritance requires an affected parent" error applies to parent relationships only.
- X-linked recessive: the "needs the mother" error applies only when the relative is the father. The existing
  proband-sex checks in `analysis/models/nodes/family_inheritance.py` stay. No check on the sibling's sex - an
  affected HOM_ALT female sibling is rare but possible, and the Duo only records the proband's sex.
- Compound het with siblings is unphased - two shared HET hits in a gene may be in cis. Both halves of the
  two-pass query (`AbstractCompHetInheritance._mum_but_not_dad` / `_dad_but_not_mum`) return the same tuple, so
  the OR collapses to one branch. An unaffected sibling carrying *both* hits is not excluded - the same gap as
  `QuadCompHet` (variantgrid_private#1263) - say so in a warning.
- Absent in parent and mosaic parent have no meaning without a parent: they raise an inheritance error for a
  sibling Duo, and the editor disables the options the way it disables X-linked for a father Duo today
  (`setModeAvailability` in the editor template). `ignore_field_errors` still lets a user run them.

## Node (`analysis/models/nodes/sources/duo_node.py`)

- `AbstractDuoInheritance._get_zyg_q` reads `duo.relative`; `parent_label` becomes `relative_label`
  (still `duo.relationship_label`). `get_zygosities_method` keys the second member by that label.
- `SimpleDuoInheritance._get_parent_proband_zygosities` → `_get_relative_proband_zygosities`. Each subclass
  branches on `duo.relative_is_sibling` per the table. A small helper on `AbstractFamilyInheritance`,
  `sibling_zygosities(affected, affected_zyg)` returning `affected_zyg` or the not-HOM_ALT set, keeps the
  recessive/XLR branches one line each and is where the Quad can converge later.
- `DuoAllRecessive._has_xlinked_branch` returns `relationship != FATHER`.
- `DuoCompHet` returns `({HET}, {HET})` for both halves when the sibling is affected, `(set(), {HET})` when
  not; the method and other-filters strings say "shared HET hits, unphased" for a sibling.
- `DuoNode.get_duo_inheritance_errors` gains the sibling rules: ABSENT_IN_PARENT and MOSAIC_PARENT error on a
  sibling Duo ("needs a parent"); DOMINANT's unaffected check and XLINKED's father check are guarded by
  relationship. `get_warnings` adds the unphased comp-het warning for a sibling Duo.
- `get_zygosity_table_data`: the row key `parent` becomes `relative`, and modes whose sibling row differs from the
  parent row emit `relative_S`, `relative_S_affected` / `relative_S_unaffected` and matching
  `other_filters_relative_S` keys. The editor's `lookup()` tries, in order, `key_<rel>_<affected>`, `key_<rel>`,
  `key_<affected>`, `key` - one candidate added at the front of the existing list.
- `get_help_text` gains a sentence on the sibling pair.

## Duo model, wizard and pages

- `Duo`: `relative_is_sibling` property; `parent_is_mother` stays (used by the X-linked branch);
  `missing_parent_label` is only meaningful for a parent Duo - callers branch on `relative_is_sibling` first.
  `pedigree_icon_members = ("relative",)`; `get_preview_icon_css_class` adds `duo-sibling` for the third shape.
  `parent_details` → `relative_details`. `get_cohort_samples` returns `[relative, proband]`.
- Wizard (`analysis/views/views_wizard.py:DuoWizardView`, `analysis/forms/forms.py:UserDuoWizardForm`):
  `parent_role` → `relative_role`; `role_classes`/`role_shape_classes` gain the SIBLING entries
  (`relative-affected`, `duo-sibling`). `_roles_for_sex` already leaves non-parent roles open to either sex, so
  Sibling appears for every sample. Sex-based auto-fill stays as the default (most duos are parent + child); in
  `variantgrid/static_files/default_static/js/family_wizard_roles.js`, choosing Sibling ticks that sample's
  affected box, since that's the usual case in the issue. `role_help` mentions the sibling option.
- Pedigree symbol `node-icon-duo` in `uicore/templates/uicore/tags/svg_icon_sprite.html`: add a sibling shape
  (overlapped square and circle like the parent, sex unknown) beside the proband under a sibship bar, stroked
  `none` by default; `.duo-sibling` in `variantgrid/static_files/default_static/css/global.scss` hides the
  parent shapes and connector and shows the sibling shape, `.relative-affected` (renamed from
  `.parent-affected`) fills whichever relative shape is visible. `--pedigree-parent-fill` keeps its name for the
  Trio/Quad.
- `snpdb/templatetags/model_tags.py:duo_table` already labels the row with `relationship_label`; the footer in
  `snpdb/templates/snpdb/tags/duo_table.html` reads "No parent was sequenced - de novo can't be confirmed and
  compound het is unphased" for a sibling Duo.
- `snpdb/serializers.py:DuoSerializer`: `parent` → `relative` (the node editor's `displayDuoInfo` reads it).
  `snpdb/grids.py:DuosListColumns` picks up the new choice through `DuoRelationship.choices`.
- Editor template `analysis/templates/analysis/node_editors/duonode_editor.html`: `.parent`/`.zyg-parent`/
  `.other-parent` → `relative`; the "Require parent zygosity calls" label takes the relationship label the way
  `.parent-role-zyg` does; `setModeAvailability` disables `N` and `M` when `relationship === 'S'`.
- Test fixture `snpdb/tests/utils/fake_cohort_data.py:create_fake_duo` takes `relative_affected`; the second
  sample is named after the relationship.

## Tests

`analysis/tests/test_duo_node.py` (sibling Duo over the same two samples, as `duo_father` is today):
recessive affected sibling wants HOM_ALT in both and rejects a HET sibling; recessive unaffected sibling rejects
HOM_ALT; dominant unaffected sibling raises no error and drops shared variants; X-linked with a sibling raises no
"needs the mother" error; all-recessive with a sibling keeps the X branch; comp-het with an affected sibling
needs the sibling HET on both hits and warns it's unphased; absent-in-parent and mosaic error on a sibling Duo.
`analysis/tests/test_family_wizard.py`: Sibling is on offer for both sexes; posting Sibling stores relationship
`S` and the affected tick. `snpdb/tests` URL test renders the view page for a sibling Duo.

## Docs

- `analysis/CLAUDE.md` wizard bullet: a Duo's relationship is Mother, Father or Sibling; the parent-only modes
  and what "relative" means.
- `claude/domain.md` Trio / Duo / Quad entry and the `Duo` docstring: "proband + one relative (a parent or a
  sibling)".
- `claude/research/analysis.md` DuoNode paragraph.
