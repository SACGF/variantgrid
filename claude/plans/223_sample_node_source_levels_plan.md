# private#223 — SampleNode at patient / specimen / extraction level: remaining work

[SACGF/variantgrid_private#223](https://github.com/SACGF/variantgrid_private/issues/223). Phase 7 of
`tso500_overall_plan.md` landed the extraction level (PR #1718). User testing on that raised two things:

> SampleNode - I don't like the interface here - where you choose extraction/sample. Not sure how better
> to do it though. Also, we need to be able to select by patient, and specimen.

> The sample icon is in effect a patient icon (always has been). I do like it showing the patient icon …
> We now have the concept of "chips" - might be good to show eg if you have a patient what
> extractions/samples you are using, if extraction how many samples etc.

Mockups for the decisions below are in `claude/mockups/223_sample_node_canvas.html` (card, icon options,
add-node menu) and `claude/mockups/223_sample_node_editor.html` (editor, current vs proposed).

---

## Where things stand

Everything level-blind already works and is tested (`analysis/tests/test_sample_node_extraction.py`):
query-time resolution ORing `pk IN (subquery)` per sample, per-sample thresholds, node-level VCF filter
ids resolved per VCF, exclusion reporting, live counts, grid Source column, clone.

What is level-specific, and where Specimen/Patient stop today:

| Piece | File | State |
|---|---|---|
| `SampleNodeSourceLevel` enum, `extraction`/`specimen`/`patient` FKs | `analysis/models/enums.py:23`, `sample_node.py:44-49` | all four exist |
| `get_sample_group()` | `sample_node.py:94` | Sample + Extraction resolve; others return an empty group |
| `IMPLEMENTED_SOURCE_LEVELS` | `sample_node.py:73` | `{SAMPLE, EXTRACTION}` → "not supported yet" config error |
| `SampleNodeForm` | `forms_nodes.py:837` | `specimen`/`patient` in `Meta.exclude`; level is a `<select>` that resubmits the form to swap widgets |
| `samplenode_editor.html` | `analysis/templates/analysis/node_editors/` | hardcoded `#id_extraction`, `Urls.extraction_samples`, `Urls.view_extraction` |
| Preview endpoint | `patients/views.py:226` `extraction_samples` | extraction only |
| VCF FILTER at group level | `analysis_node.py:1520` `get_filter_codes`, `views.py:930` `extraction_vcf_locus_filters` | one node-level selection of filter *names* translated into every VCF — wrong, since only PASS is shared |
| `get_node_icon()` | `sample_node.py:416` | reads `self.sample.patient` only — group nodes get the plain square |
| `get_node_chips()` | `sample_node.py:425` | extraction chip (`fa-vial`, borrowed from specimen) + one chip per VCF |
| Add-node menu | `node_types.py:20`, `views_json.py:116` | one entry per class — a Patient node is "add Sample, change dropdown" |
| `analysis_templates_tag` | `related_analyses_tags.py:139` | `single_model_args = {sample, cohort, trio, quad, pedigree}` |

---

## Design

### 1. One picker, one tree — the editor

The level dropdown goes. The editor gets a single search box that finds any of the four kinds of thing
and a tree underneath showing the whole patient with the chosen row highlighted:

```
 Find patient, specimen, extraction or sample…         [ 🔍 ]

 ● Jane Citizen (F)                    PATIENT · 7 samples in 7 VCFs
   ├ ○ Liver biopsy FFPE   2024-03-01   SPECIMEN · 5 samples
   │   ├ ○ DNA  2024-03-04             EXTRACTION · 3 samples
   │   │   ├ TSO500_DNA_hard-filtered   93 variants   AD≥[ ] DP≥[ ] …
   │   │   ├ TSO500_DNA_cnv             16 variants   AD≥[ ] DP≥[ ] …
   │   │   └ TSO500_DNA_exon_cnv         4 variants   …
   │   └ ○ RNA  2024-03-04             EXTRACTION · 2 samples
   │       ├ TSO500_RNA_splice          12 variants
   │       └ TSO500_RNA_cnv          ⚠ GRCh37 (need GRCh38)
   └ ○ Blood  2023-11-12                SPECIMEN · 2 samples
       └ ○ DNA                          EXTRACTION · 2 samples
           ├ BloodWGS                    4.1M variants
           └ BloodWGS_sv             ⚠ VCF archived
```

- The radio on a patient / specimen / extraction row *is* the level control: picking a row sets
  `source_level` and the matching FK. Sample rows are picked the same way (level = SAMPLE).
- Rows outside the picked subtree stay visible but dimmed — that is the "what am I bringing in"
  answer, and it is how you move up or down a level without re-searching.
- Excluded samples show their reason inline (`⚠ VCF archived`, `⚠ GRCh37 (need GRCh38)`), the same
  strings `SampleGroup.warnings` already produces. Samples the user cannot see are one dimmed
  "2 samples hidden (no permission)" row per extraction.
- **Filters are global by default, per-VCF on demand.** Only `PASS` means the same thing in every VCF,
  so the global block under the tree is: `PASS only` checkbox, `AD`/`DP`/`GQ`/`PL`, zygosity, allele
  frequency. Sample rows are collapsed; expanding one shows that VCF's own FILTER codes (PASS greyed as
  the global setting) and its threshold overrides. Inputs show the global value as a placeholder and go
  bold once overridden, so a row with nothing typed inherits — exactly the semantics
  `get_sample_thresholds()` (`sample_node.py:170`) already has. A collapsed row with anything set shows
  a summary (`AD≥5 · DP≥20 · +LowDepth`). The separate per-sample table goes; the hidden
  `sample_thresholds` JSON field and `SampleThresholdsMixin` stay as they are.
- The summary line ("7 samples in 7 VCFs · 1 excluded") and the `View patient / specimen / …` cross-link
  come from the same payload.
- Locked-input templates (`lock_input_sources` on a non-template analysis) render the tree read-only:
  search box hidden, radios disabled; the global filters and per-row overrides stay editable.

A **grouped autocomplete** backs the search box. `dal_queryset_sequence` is not installed and is not
needed — select2 renders a result whose `children` array is a group natively, so a
`SampleSourceAutocompleteView(Select2QuerySetView)` overrides `get()` to return

```json
{"results": [
   {"text": "Patients",    "children": [{"id": "T:12",  "text": "Jane Citizen (F)"}]},
   {"text": "Specimens",   "children": [{"id": "P:40",  "text": "Liver biopsy FFPE — Jane Citizen"}]},
   {"text": "Extractions", "children": [{"id": "E:71",  "text": "DNA 2024-03-04 — Liver biopsy FFPE"}]},
   {"text": "Samples",     "children": [{"id": "S:903", "text": "TSO500_DNA_hard-filtered (TSO500_run7.vcf)"}]}
 ], "pagination": {"more": false}}
```

Each group reuses the corresponding existing autocomplete's `get_user_queryset` (so permissions, import
status, `exclude_archived` and the genome build forward all apply — see §3 for the two that need the
build filter added), capped at ~8 rows per group. Ids are `<level>:<pk>` using the enum values.

### 2. The card

Keep the pedigree badge at every level. The badge has always drawn the *patient* (square/circle,
struck through if deceased), and a patient is what every level resolves to, so
`get_node_icon()` derives the patient from whichever object is set:

```python
def _get_patient(self) -> Optional[Patient]:
    match self.source_level:
        case SampleNodeSourceLevel.SAMPLE:     return self.sample.patient if self.sample else None
        case SampleNodeSourceLevel.EXTRACTION: return self.extraction.specimen.patient if self.extraction else None
        case SampleNodeSourceLevel.SPECIMEN:   return self.specimen.patient if self.specimen else None
        case SampleNodeSourceLevel.PATIENT:    return self.patient
```

What the node actually *is* goes on the class strip, which is where a card says its kind. The strip text
currently comes from the classmethod `get_node_class_label_short()`; add an instance hook

```python
def get_node_strip_label(self) -> str:   # AnalysisNode: return self.get_node_class_label_short()
```

that `get_rendering_dict()` uses for `class_label_short`; the menu keeps using the classmethod.
`SampleNode` returns `SAMPLE` / `EXTRACTION` / `SPECIMEN` / `PATIENT`. The strip has ~45px of text at
7px bold, so `EXTRACTION` needs a `.node-klass-long` rule (6px, tighter tracking) that `analysis_nodes.js`
adds when the label is over 8 characters; the mockup's strip row shows it fitting.

**Chips are the hierarchy from the picked object down, nested the way the relations are** —
specimen ⊃ extraction ⊃ one `VCF ×N` leaf — and the picked object is the outermost chip, so an
extraction node carries the flask and a sample node the VCF icon:

| Level | Chips |
|---|---|
| Sample | `▮ VCF` — hover is the VCF it came from |
| Extraction | `⚗ DNA` group holding `▮ VCF ×3`; `⚠ ×N` |
| Specimen | `💧 FFPE` group holding one `⚗` group per extraction, each holding its `VCF ×N`; `⚠ ×N` |
| Patient | one `💧` group per specimen, as above; `⚠ ×N` |

VCFs are always a single `VCF ×N` chip whose hover lists the names. Group labels: tissue name, else
tissue status, else reference id for a specimen; nucleic acid for an extraction; hover is the full
`str()`. The `⚠ ×N` chip carries `css_class="node-chip-warning"` with the exclusion reasons as its
title — the "six months later, why fewer rows" question answered on the canvas. Above 3 extractions or
specimens in one group they collapse to a `×N` label — `MAX_VCF_CHIPS = 3` becomes the per-group cap.

Two additions to `NodeChip` (`node_display.py`) carry this:

```python
@dataclass(frozen=True)
class NodeChip:
    text: str
    icon: Optional[str] = None
    title: Optional[str] = None
    css_class: Optional[str] = None
    count: Optional[int] = None                  # small bubble after the text, eg ×3
    children: tuple["NodeChip", ...] = ()        # renders as a group wrapping its children
```

`renderChip()` in `analysis_nodes.js` draws `count` as `<span class="node-chip-count">` and a chip with
`children` as a `<span class="node-chip-group">` holding a label plus the child chips, recursively. New
rules in `analysis_nodes.css` beside `.node-chip`: `.node-chip-count`, `.node-chip-group` (full card
width, grey, children wrapping under the label; a nested group one shade lighter), `.node-chip-warning`
(amber) and `.node-chip-text:empty { display: none }` so the ⚠ chip is icon + count. `asdict()` already
serialises nested dataclasses, so `get_rendering_dict()` is unchanged.

Chip icons follow the models' own `preview_icon()`s, so chips, hover cards and search agree. Specimen
becomes `fa-solid fa-droplet` (`Specimen.preview_icon()`, `patients/models.py:364`) — biological material
against the extraction's lab glassware, and a blob against a stick at chip size; extraction keeps
`fa-solid fa-flask`. VCF keeps `fa-file-lines`. This also retires the current extraction chip's borrowed
vial (`sample_node.py:434`). Chip data comes from
`get_sample_group()`, already `cache_memoize`d via `get_source_vcf_names` — widen that memo to the tree
(`specimen → extraction → vcf names`, plus `excluded_reasons`) so the card costs one lookup; at patient
level it is the same shape the editor tree (§3) is built from.

`get_node_name()` (`sample_node.py:343`) already does `str(source_object)` + `(N samples)` for every
level; `_get_method_summary()` likewise. The editor's Patient tab (`samplenode_editor.html:120`)
switches from `sample.patient` to `node._get_patient()` so it appears at every level.

### 2b. VCF FILTER — PASS is global, everything else is per VCF

`NodeVCFFilter(node, vcf_filter)` (`analysis_node.py:1505`) already stores per-VCF selections: a
`VCFFilter` belongs to one VCF, and `vcf_filter=None` means PASS. What is wrong today is the group-level
reading of it — `NodeVCFFilter.get_filter_codes(node, vcf)` (`:1520`) translates the selected filter
*names* into every VCF's codes, as if `LowDepth` in a DRAGEN small-variant VCF meant the same as
`LowDepth` in a CNV VCF. The user-facing rule becomes:

> records read from a VCF = `PASS` (if the node's PASS row exists) ∪ that VCF's own ticked codes;
> no PASS row and nothing ticked for a VCF means that VCF is unfiltered.

So `get_filter_codes(node, vcf)` selects `NodeVCFFilter.objects.filter(node=node)` rows whose
`vcf_filter__vcf == vcf`, plus `None` if a PASS row exists — no name translation. The sample-level case
is unchanged, since a single-VCF node's rows all belong to that VCF. `SampleNode.get_filter_code()`
(`sample_node.py:294`) keeps its 0 / 1 / 2 meaning (none / PASS only / mixed) computed from the rows.

`VCFLocusFiltersMixin` (`forms_nodes.py:85`): the hidden `vcf_locus_filters` JSON becomes
`{"pass": true, "by_vcf": {"<vcf_id>": ["LowDepth", …]}}`; `save_vcf_locus_filters` writes one PASS row
when `pass` is set and one row per `(vcf, filter_id)` found in `VCFFilter.objects.filter(vcf_id=…,
filter_id=…)`, ignoring VCFs outside the node's group. At sample level the editor still shows the single
VCF's full filter list at the top (one VCF, so global and per-VCF coincide) and posts it under `by_vcf`.
The group branch of the `vcf_locus_filter` tag (`vcf_locus_filter_tags.py:6`) and the
`extraction_vcf_locus_filters` view (`analysis/views/views.py:930`, `analysis/urls.py:70`) are retired —
each sample row's filter list comes from the tree payload (§3). No migration.

### 3. Resolvers, endpoints, autocompletes

`patients/sample_grouping.py`:

```python
def get_sample_group(user, level: SampleNodeSourceLevel, source, genome_build=None) -> SampleGroup
```

with the existing `get_extraction_sample_group` becoming the `EXTRACTION` arm. The candidate querysets:

| Level | Candidate samples |
|---|---|
| Extraction | `Sample.objects.filter(extraction=e)` (as now) |
| Specimen | `Sample.objects.filter(extraction__specimen=s)` |
| Patient | `Sample.objects.filter(Q(patient=p) \| Q(extraction__specimen__patient=p)).distinct()` |

Patient is a union because the two links are set independently — `link_samples_and_vcfs_to_sequencing`
(`upload/vcf/vcf_import.py:~543`) carries `extraction` down and never sets `sample.patient`, while the
patient CSV import sets `patient` and may leave `extraction` null. `Patient.get_samples()`
(`patients/models.py:299`) only follows the direct FK; leave it alone (it drives the patient page) and
keep the union local to the grouping. Exclusion rules (`import_status`, archived, build, permission)
are shared and unchanged.

A **tree builder** for the editor, in the same module:

```python
def get_patient_sample_tree(user, patient, genome_build, selected: tuple[level, pk]) -> dict
```

walking `patient.specimen_set` → `extraction_set` → `sample_set`, plus an "Unlinked samples" pseudo
extraction for samples with `patient=p` and no extraction. Every sample row carries the same fields
`sample_group_as_json` emits today (`sample_id`, `sample`, `vcf_id`, `vcf`, `genome_build`,
`has_genotype`, `variant_count`) plus `included: bool`, `reason`, and
`vcf_filters: [{"filter_id", "description"}]` from `vcf.vcffilter_set` so the expanded row can draw that
VCF's FILTER checkboxes; every container row carries `level`, `id`, `label`, `sample_count`,
`in_selection: bool`. A specimen / extraction / sample id
resolves to its patient first so the tree is always the full patient.

Endpoints in `patients/views.py` / `urls.py`:

- `GET /patients/sample_group/<level>/<pk>/samples?genome_build=` — generalises `extraction_samples`
  (keep the old URL name pointing at the extraction arm; `view_extraction.html` links to it).
- `GET /patients/sample_group/<level>/<pk>/tree?genome_build=` — the editor's payload.
- `GET /patients/sample_source_autocomplete/` — the grouped search from §1.

Autocompletes (`patients/views_autocomplete.py`): `SpecimenAutocompleteView` and
`PatientAutocompleteView` become `GenomeBuildAutocompleteView`s filtering on
`extraction__sample__vcf__genome_build` / `Q(sample__vcf__genome_build) | Q(specimen__extraction__sample__vcf__genome_build)`
with `exclude_archived` on the same paths, `.distinct()`, matching what `ExtractionAutocompleteView`
does. They keep serving the specimen/extraction pages; the grouped view composes them.

### 4. Form

`SampleNodeForm` (`forms_nodes.py:837`):

- Drop `source_level`, `sample`, `extraction` from the rendered form; add
  `source = forms.CharField(widget=autocomplete.Select2(url="sample_source_autocomplete", forward=…))`
  carrying `"<level>:<pk>"`. `clean_source` splits it, checks the level is in
  `IMPLEMENTED_SOURCE_LEVELS`, loads the object through its model's `get_for_user`, and `save()` sets
  `source_level` plus the one FK, nulling the other three. `initial["source"]` is built from the
  instance so a saved node round-trips; `SampleNodeSourceLevel.SAMPLE` remains the model default.
- `genome_build_fields = ["source"]` so the build and `exclude_archived` forward to the grouped view.
- `SAMPLE_LEVEL_FIELDS` (`sample_gene_list`, `restrict_to_qc_gene_list`) stay hidden at group levels;
  `LOCKED_INPUT_FIELDS` becomes `["source", "restrict_to_qc_gene_list"]`.
- `Meta.exclude` drops `"specimen", "patient"` and gains `"sample", "extraction", "source_level"`.
- `IMPLEMENTED_SOURCE_LEVELS` = all four; the "not supported yet" branch in
  `_get_configuration_errors` goes.

Analysis variables: `AnalysisVariable` is keyed `(node, field)` and `populate_arguments`
(`models_analysis.py:646`) sets fields by name, so a template node saved at specimen level with
`specimen` marked variable keeps working unchanged. The `source` field itself is never the variable —
the variables UI keeps offering the four FKs.

### 5. Add-node menu

Four entries under Source — Sample, Extraction, Specimen, Patient — each creating a `SampleNode` with
`source_level` stamped, so a user looking for "Patient" finds it. Mechanism, per the earlier plan:

- `AnalysisNode.get_menu_entries()` classmethod returning `[(node_type_key, label, icon, initial_kwargs)]`,
  default one entry `(class_name, label, class_icon, {})`; `SampleNode` returns one per level with key
  `f"SampleNode:{level}"`, label the level's display, and a per-level `NodeIcon` for the menu only
  (pedigree square for Sample, `fa-flask`, `fa-droplet`, `fa-user-injured` — the models' own icons).
- `get_nodes_by_classification()` / `NODE_TYPES` display data and `_get_node_types_choices()` iterate
  entries rather than classes. `node_create` (`views_json.py:116`) splits the key on `:` and passes
  `initial_kwargs` to `objects.create`. The `SampleNode` entry list order is Sample, Extraction,
  Specimen, Patient; `initial = "SampleNode"` still resolves.
- A deployment without patients data can trim the list via a `settings.ANALYSIS_SAMPLE_NODE_LEVELS`
  tuple defaulting to all four — the "switched off" option from the earlier plan.

### 6. Launching templates from the patient, specimen and extraction pages

`analysis_templates_tag` (`related_analyses_tags.py:139`): `single_model_args` gains `extraction`,
`specimen`, `patient`; `_source_is_archived` (`:124`) treats those three as archived when every sample
in their group is archived. Then:

- `view_extraction.html` and `view_specimen.html` get an analysis templates section beside their Samples
  table, keyed on the analysis build (use the tag once per build present in the group, as
  `related_analyses_for_samples` does).
- `view_patient.html:176` switches `show_create_analyses=False` on and passes `patient=`.

`AnalysisTemplate` filtering is already by `class_name`, so a template whose variable is a `SampleNode`
`specimen` field lists on the specimen page with nothing else changing.

---

## Files

| Area | Files |
|---|---|
| Resolvers, tree, JSON | `patients/sample_grouping.py`, `patients/views.py`, `patients/urls.py` |
| Autocompletes | `patients/views_autocomplete.py`, `patients/urls.py` |
| Node | `analysis/models/nodes/sources/sample_node.py`, `analysis/models/nodes/analysis_node.py` (`get_node_strip_label`), `analysis/models/nodes/node_utils.py`, `analysis/models/nodes/node_types.py`, `analysis/views/views_json.py`, `analysis/forms/forms.py` |
| Form + editor | `analysis/forms/forms_nodes.py`, `analysis/views/nodes/sample_node_view.py`, `analysis/templates/analysis/node_editors/samplenode_editor.html`, `analysis/templatetags/vcf_locus_filter_tags.py`, `analysis/templates/analysis/tags/vcf_locus_filter_tag.html`, `analysis/views/views.py`, `analysis/urls.py` |
| VCF FILTER | `analysis/models/nodes/analysis_node.py` (`NodeVCFFilter.get_filter_codes`), `analysis/models/nodes/sources/sample_node.py` (`get_filter_code`, `has_filters`) |
| Canvas | `analysis/models/nodes/node_display.py` (`count`, `children`), `variantgrid/static_files/default_static/js/analysis_nodes.js` (`renderChip`), `variantgrid/static_files/default_static/css/analysis_nodes.css` (`.node-chip-count`, `.node-chip-group`, `.node-chip-warning`, `.node-klass-long`) |
| Templates tag + pages | `analysis/templatetags/related_analyses_tags.py`, `patients/templates/patients/view_extraction.html`, `view_specimen.html`, `view_patient.html` |
| Settings | `variantgrid/settings/components/default_settings.py` (`ANALYSIS_SAMPLE_NODE_LEVELS`) |

No migration: the columns and enum values already exist (`0111_sample_node_source_level`).

## Tests

Extend `analysis/tests/test_sample_node_extraction.py` (rename to `test_sample_node_levels.py`) with a
specimen fixture holding DNA + RNA extractions and a patient with a second specimen plus one
`patient`-only sample:

- specimen level resolves both arms; patient level resolves both specimens *and* the unlinked sample,
  once each (`distinct`)
- build / archived / permission exclusions surface at specimen and patient level
- `get_node_icon()` gives the pedigree glyph at every level; strip label per level
- chips: sample → one `VCF` chip; extraction → `DNA` group holding `VCF ×3`; specimen → `FFPE` ⊃ `DNA`/`RNA` ⊃ `VCF ×N`; patient → one group per specimen; `⚠ ×N` with reasons in title; `test_node_display` gains a check that `children` chips serialise
- form: `"P:40"` sets `source_level=SPECIMEN` + `specimen`, nulls the rest; unimplemented / unviewable
  ids are form errors; locked-input template renders read-only
- tree JSON: full patient from a sample id, `in_selection` flags, unlinked pseudo-extraction, each sample's `vcf_filters`
- VCF FILTER: global PASS + `LowDepth` ticked on VCF A → VCF A reads PASS ∪ LowDepth, VCF B reads PASS only; nothing set → both unfiltered; a `LowDepth` in VCF B is never selected by VCF A's row (replaces `test_node_level_filter_resolves_into_each_vcfs_codes` and `test_one_filter_selection_is_stored_for_a_node_spanning_vcfs`)
- grouped autocomplete: four groups, build and archive forwards honoured
- menu: four `SampleNode:*` entries, `node_create("SampleNode:T")` stamps `PATIENT`
- `analysis_templates_tag` accepts `patient=` / `specimen=` / `extraction=`

Update `test_chips_show_extraction_and_vcfs` for the nested `DNA [VCF ×3]` shape, and `patients` preview tests for the droplet, and
`test_node_display` still passes since `get_node_chips()` on an unsaved node returns `[]`.

## Done when

On the TSO 500 test run: add a Patient node from the menu, search the patient, see the full tree, pick
the liver specimen row — the card reads `SPECIMEN` under a female pedigree badge with an `FFPE` group holding `DNA [VCF ×3]` and
`RNA [VCF ×2]`, the grid shows all five callers' rows with a Source column, per-caller `min_ad` and a
`LowDepth` tick under one VCF's row apply to that VCF alone, and switching the radio to the DNA extraction row narrows it to three without
re-searching. The extraction page lists the analysis templates and launches one.
