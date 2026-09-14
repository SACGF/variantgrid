# Multi-variant case report - HTML, PDF, DOCX and JSON from a case's classifications (#444)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-14
Status: draft

Design for [#444](https://github.com/SACGF/variantgrid/issues/444) (multi-variant classification + reporting), the
reporting half of [sapath#431](https://github.com/SACGF/variantgrid_sapath/issues/431) (TSO 500), and the "#444
remainder" of Phase 8 in [`tso500_overall_plan.md`](tso500_overall_plan.md). It builds on what sapath#246 landed:
the Classify & Report tab (`analysis/classify_report.py`, `analysis/views/views_classify_report.py`) and
`ClassificationReport` (`classification/views/classification_export_report.py:ClassificationReport`) accepting a list of
classifications with `gene_groups` in its context.

Scope: the somatic (AMP) case report. The template and `CaseReport` models are generic so a germline case report is a
later template on the same machinery; this plan builds the somatic one.

Worked example: the three files in `../variantgrid_sapath/sapath/test_data/tso500/` - one de-identified printed TSO 500
report (`.txt`) and two JSON records in the shape the `../tso500_data` project extracts from the PDFs
(parse_pdf_to_json in TSO_report_parser1.py there). The report the lab prints has this structure:

```
header          patient, referrer, accession, specimen, tumour cell content, clinical indication
Results Summary table   small variants (tier, gene, protein, VAF) / copy number changes (tier, gene, copies) /
                        fusions / TMB / MSI
Summary Interpretation  one case-level paragraph
Variant Interpretation  Tier I -> Tier II -> Tier III; within a tier by gene; variants of one gene joined
                        with "and", each variant's narrative, then one gene-level paragraph
Comment / Therapy matching / Method / Genes tested / sign-off   mostly fixed text plus versions
```

## What the report run produces

Every run renders one server-side HTML document and derives the rest from it, all in pure Python in the request:

| Output | How | Who uses it |
|---|---|---|
| HTML | the report template (a Django template, same engine as the single-record report) rendered over the context | the scientist, as the preview in the browser |
| `.pdf` | [xhtml2pdf](https://xhtml2pdf.readthedocs.io/) converts that HTML - pure Python on reportlab, already in `requirements.in` | the read-only copy that goes with the case |
| `.docx` | [html2docx](https://pypi.org/project/html2docx/) converts the same HTML - pure Python on python-docx | medical scientists, when the LIS wants edited text |
| `.json` | a second Django template renders the same context to JSON; blank = the canonical context dump | downstream systems (the `tso500_data` variant database) |

One HTML template is the lab's whole report design, so the preview, the PDF and the Word file cannot disagree, and
the lab maintains the template the way it maintains `sapath/default_report.html` today. The report's structure is
headings, paragraphs and tables, which is what both converters do well; xhtml2pdf also gives `@page` headers, footers
and page breaks. The lab's final report is assembled in the LIS (the `.txt` example carries the LIS's NATA footer,
report ID and print date), so what VariantGrid produces is the content a scientist pastes or attaches, and the Word
file is there so they can touch it first. JSON comes from the same context so the structured record and the human
document describe the same classifications at the same versions.

The templates only loop: ordering, grouping and tier derivation live in Python (§Context). Rendering is server-side
because PDF and DOCX have no browser to run Vue in, so the case template reads `{{ v.evidence.c_hgvs.value }}` where
the single-record template reads `:evidence="c_hgvs"`.

---

## Data

### `ClassificationReportTemplate` - grows the multi-variant outputs

`classification/models/classification_report_models.py:ClassificationReportTemplate`. The existing `template` column
stays as the single-record HTML/Vue report (`view_template_report`); the new columns are the case report.

```python
class ClassificationReportTemplate(TimeStampedModel):
    name = models.TextField(primary_key=True)
    template = models.TextField(blank=True, default="")                       # existing: single-record HTML (Django + Vue)
    case_template = models.TextField(blank=True, default="")                  # Django template -> HTML, server-rendered; PDF and DOCX derive from it
    json_template = models.TextField(blank=True, default="")                  # Django template -> JSON; blank = canonical context dump
    case_fields = models.JSONField(default=list, blank=True)                  # case-level inputs the build form asks for
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices, null=True, blank=True)
```

- `case_template` present means the template can build a case report; the tab offers only those.
- `case_template` and `json_template` are validated on save: rendered against a fixture context, the JSON must
  `json.loads` and the HTML must convert through xhtml2pdf without error.
- `case_fields` is a list of `{"key", "label", "type": "text" | "bool" | "choice", "options": [...], "default"}`. It is
  how a deployment adds report-level inputs without a schema change - SA Path's assay-success and caveat flags
  (`Amplifications`, `Deletions`, `Variants`, `Fusions`, `TMB`, `MSI`, `Confirmed`; `Purity`, `Quality`, `Fail`) and
  the "Mutations Comment" are all `case_fields` on the TSO 500 template.
- `allele_origin_bucket` filters which templates a case is offered (`preferred_template_for` already says this is
  coming). Null = offered to every case. Templates stay global rather than per-lab, per the history on the issue.

### `CaseReport` - one run, pinned

New, in `classification/models/classification_report_models.py`. Classification may reference patients and snpdb,
and needs nothing from analysis, so it sits beside the template it renders. The case is a level of the
`Patient → Specimen → Extraction → Sample` hierarchy, spelt the way `SampleNode` spells it: `source_level` plus one
nullable FK per level (`patients/models_enums.py:SampleSourceLevel`, resolved through
`patients/sample_grouping.py:get_sample_group`).

```python
class CaseReportStatus(models.TextChoices):
    DRAFT = 'D', 'Draft'
    FINAL = 'F', 'Final'
    SUPERSEDED = 'S', 'Superseded'


class CaseReport(TimeStampedModel):
    template = models.ForeignKey(ClassificationReportTemplate, on_delete=PROTECT)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=PROTECT)
    status = models.CharField(max_length=1, choices=CaseReportStatus.choices, default=CaseReportStatus.DRAFT)
    supersedes = models.ForeignKey('self', null=True, blank=True, on_delete=SET_NULL)

    source_level = models.CharField(max_length=1, choices=SampleSourceLevel.choices)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    specimen = models.ForeignKey(Specimen, null=True, blank=True, on_delete=SET_NULL)
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)

    summary = models.TextField(blank=True)                    # the case-level "Summary Interpretation"
    case_values = models.JSONField(default=dict, blank=True)  # answers to template.case_fields
    external_report_id = models.TextField(null=True, blank=True)   # the LIS's report ID, entered when known
    report_date = models.DateField(null=True, blank=True)

    context_snapshot = models.JSONField()                     # the rendered context - the run is reproducible
    html = models.TextField()                                 # the rendered document; PDF and DOCX are derived from it
    pdf_file = models.FileField(upload_to=case_report_upload_path)
    docx_file = models.FileField(upload_to=case_report_upload_path)
    json_output = models.JSONField(null=True, blank=True)     # the rendered JSON, queryable

    class Meta:
        indexes = [models.Index(fields=["specimen"]), models.Index(fields=["patient"]), models.Index(fields=["sample"])]
        constraints = [models.CheckConstraint(name="case_report_one_source",
                                              check=exactly one of patient / specimen / extraction / sample is set)]


class CaseReportClassification(models.Model):
    case_report = models.ForeignKey(CaseReport, on_delete=CASCADE)
    classification_modification = models.ForeignKey(ClassificationModification, on_delete=PROTECT)  # the pinned version
    order = models.IntegerField()                             # position in the report, as computed at build time
    reported = models.BooleanField(default=True)              # "Report": Y/N - in the document vs listed as seen only

    class Meta:
        unique_together = ("case_report", "classification_modification")
        ordering = ("order",)
```

- A `CaseReport` pins **modifications**, so a later edit to a classification changes nothing about a report already
  issued; the tab's existing "unsubmitted changes" warning is what tells the scientist a rebuild is due.
- A `FINAL` report is immutable. Rebuilding makes a new `CaseReport` with `supersedes` set and marks the old one
  `SUPERSEDED`. Drafts can be deleted by their lab.
- `context_snapshot` holds everything the templates saw (measures included), so the same documents can be re-rendered
  after a template fix, and so a report's numbers are inspectable without re-deriving them from rows that may have
  changed (a `SpecimenMeasure` resend replaces the row).
- `reported` is the per-variant Report Y/N the JSON carries. It defaults from the `variant_reported` evidence key
  (`not_included` → False) and the scientist can flip it in the build form; the Tier section prints "No reportable
  variants detected" when a tier has variants but none reported.
- Files live under `MEDIA_ROOT/case_reports/<lab pk>/<report pk>/`; downloads go through a view that checks the
  user is in the report's lab or can see the case, never a direct media URL.

### The context - `ReportContext` and `ReportVariant`

Dataclasses in a new module, classification/report/case_report_context.py (the package holds the context builder, the
ordering rules and the three renderers). `ClassificationReport.context()` is refactored to build one of these too, so
the single-record HTML path and the case report share one definition of a variant's row.

```python
@dataclass
class ReportVariant:
    modification: ClassificationModification
    kind: str                     # 'small_variant' | 'copy_number' | 'fusion' | 'splice'
    gene_symbol: str              # sort key; fusions use the 5' partner
    gene_symbols: list[str]       # both fusion partners
    tier: Optional[str]           # somatic:clinical_significance value, eg 'tier_2'
    amp_tier: str                 # 'IA' 'IB' 'IIC' 'IID' 'III' 'IV' or '' - see §Tier
    tier_rank: int                # sort key derived from amp_tier then tier
    vaf: Optional[float]          # allele_frequency as a fraction
    copy_number: Optional[int]
    reported: bool
    sample: Optional[Sample]
    evidence: dict                # ClassificationReport.row_data(): every ekey as {value, note, formatted, label}
    warnings: list[str]           # eg tier / AMP level disagree


@dataclass
class GeneGroup:
    gene_symbol: str
    variants: list[ReportVariant]
    gene_summary: str             # h_summary - the one gene-level paragraph, see §Gene-level text
    gene_summary_source: Optional[int]   # pk of the modification it was taken from
    warnings: list[str]           # eg the gene's classifications disagree on h_summary


@dataclass
class TierGroup:
    tier: str                     # 'tier_1' ...
    label: str                    # 'Tier I - Variants of Strong Clinical Significance'
    genes: list[GeneGroup]


@dataclass
class KindGroup:
    kind: str
    label: str                    # 'Somatic Variants', 'Copy Number Changes', 'Gene Fusions'
    variants: list[ReportVariant]


@dataclass
class ReportContext:
    case_report: Optional[CaseReport]
    source_level: str
    patient: Optional[Patient]
    specimen: Optional[Specimen]
    extractions: list[Extraction]
    samples: list[Sample]
    sequencing_runs: list[str]
    measures: dict[str, SpecimenMeasure]   # keyed 'tmb', 'msi', 'gis', 'tumour_fraction', 'ploidy'
    variants: list[ReportVariant]          # every included classification in report order
    kind_groups: list[KindGroup]           # the Results Summary
    tier_groups: list[TierGroup]           # the Variant Interpretation
    gene_groups: list[GeneGroup]           # kept for templates written against the sapath#246 context
    summary: str
    case_values: dict
    lab: Lab
    user: User
    generated: datetime
    versions: dict                         # variantgrid, annotation version per build, DRAGEN from the VCF headers
```

Templates receive `asdict()` of this, with model instances replaced by small dicts (`pk`, `str`, the fields listed in
§Fields). That dict is also `context_snapshot`, and is the canonical JSON when `json_template` is blank.

---

## Ordering - the rules, in Python

The report order is fixed by the context builder; templates loop in the order given.

1. **Kind**: small variants, then copy number changes, then fusions, then splice variants. `kind` comes from the
   `Variant` the classification resolved to - a gene-level alt (`snpdb/gene_level_variants.py`, `GeneFusion` in
   `genes/models/models_gene_fusion.py`) says copy number or fusion; a classification with no resolved variant falls back
   to the `variant_class` evidence key (`copy_number_gain` / `copy_number_loss` → copy number, otherwise small variant).
2. **Tier**: `IA`, `IB`, `IIC`, `IID`, `III`, `IV`, then unclassified. Derived per §Tier.
3. **Gene symbol**: alphabetical.
4. **VAF**: descending, then copy number descending, then c.HGVS - so two variants in one gene print highest VAF
   first, as the example does.

`kind_groups` applies 1 then 2-4 within each kind (the Results Summary). `tier_groups` applies 2, then groups by gene
(3), then 4 within the gene, with kinds interleaved (a Tier IIC amplification prints under Tier II beside the Tier IIC
small variants, as the example's `GENE3 amplification (11 copies). Tier IIC.` does). A fusion is listed under its 5'
partner and its `gene_symbols` carries both, so a template can print `GENE1::GENE2`.

### Tier

`somatic:clinical_significance` holds the tier (`tier_1`, `tier_2`, `tier_3`, `tier_4`, `tier_1_or_2`;
`classification/enums/classification_enums.py:SomaticClinicalSignificance`) and `amp:level_a` … `amp:level_d` hold the
AMP evidence levels. The printed report needs the sub-tier (`IIC`), which is the two together:

| tier | AMP level populated | `amp_tier` |
|---|---|---|
| `tier_1` | A (any value) | `IA` |
| `tier_1` | B, and A empty | `IB` |
| `tier_2` | C (any value) | `IIC` |
| `tier_2` | D, and C empty | `IID` |
| `tier_3` | - | `III` |
| `tier_4` | - | `IV` |
| `tier_1` / `tier_2` with no matching level, or `tier_1_or_2` | - | the bare tier (`I`, `II`, `I/II`) plus a warning |

The warning surfaces on the build form beside the variant and in `ReportVariant.warnings`, so a report is never built
with a silently wrong sub-tier. The tier + level schema is the one settled for Shariant (labs supply whichever parts they
have), SA Path records its somatic work in it, and the report derives the sub-tier from it rather than adding a key.

### Gene-level text

The example prints one paragraph per gene after that gene's variants ("The GENE ABOVE does gene level things"). Gene
scope in the evidence keys (`EvidenceKey.copy_scope`, `classification/models/evidence_key.py`) currently covers
`h_summary`, `mode_of_inheritance`, `mechanism_of_disease` and the other gene/disease keys, all germline-facing.
`GeneGroup.gene_summary` is `h_summary`. Gene-level copy consensus (#1419) means the gene's classifications normally
carry the same text; where they differ, the builder takes the value from the most recently modified classification,
records which one in `gene_summary_source`, and adds a warning naming the classifications that disagree. The build form
shows the warning beside the gene with links to those records, whose form already offers the "newer gene content"
helper (`classification/views/views_gene_consensus.py`) to bring them into line before the report is finalised.

---

## Fields - what the report needs and where it comes from today

The TSO 500 report and JSON, field by field. "Have" means the value exists on a model or evidence key the builder can
read; the gap column is what this plan adds or what stays outside VariantGrid.

### Case header

| Report field | Source today | Gap |
|---|---|---|
| Patient name, DOB, sex | `Patient.first_name` / `last_name` / `date_of_birth` / `sex` | have; at SA Path also `HelixNGSOrder` |
| MRN, referring doctor, facility, location | - | LIS fields; the LIS wraps the content. Left out of VG deliberately |
| Ref / "Patient" (`C12345`, the patient ID Omico issues) | `Patient.patient_code`, posted by the client through the patient API | have; the JSON's `"Patient"` |
| Accession (`GM-26-…`) | `ExternalPK` on the specimen / `HelixAccession` at SA Path | have where the client posts it |
| Test requested / Panel (`TSO500`) | `SequencingRun.enrichment_kit` via `SequencingSample`; `FileUpload.metadata["source"]` | have when seqauto is linked; else a `case_field` |
| Seq Run (`TSO_26_024`) | `SequencingRun.name` via the extraction's `SequencingSample`s | have when linked; `sequencing_runs` in context |
| Specimen type, received, collected | `Specimen.tissue`, `received_date`, `collection_date` | have (tissue pending #1747) |
| FFPE block ID | `Specimen.reference_id` / `ExternalPK` | have |
| Tumour cell content % (purity) | `SpecimenMeasure` `TUMOUR_FRACTION`; also the per-record `somatic:tumor_cellularity` key | have at specimen level - the report reads the measure |
| Clinical indication | `sa_path:clinical_indication` per classification; `HelixClinicalIndication` at SA Path | case-level value wanted: a `case_field`, prefilled from the first classification's key |
| TMB value + Low/High | `SpecimenMeasure` `TMB` (`value`, `call`) | have |
| MSI % + Stable/Low/High | `SpecimenMeasure` `MSI` | have |
| HRD status / GIS | `SpecimenMeasure` `GIS` | have |
| Assay success flags ×7, caveats ×3 | - | `case_fields` on the template (`bool`) |
| Summary Interpretation (case level) | - (`somatic:summary_interpretation` is per variant) | `CaseReport.summary` |
| Mutations comment / "Note that…" | - | `case_field` (`text`) |
| NATA report ID, report date | `report_id` / `report_date` evidence keys are per classification | `CaseReport.external_report_id` / `report_date`, entered after the LIS issues them |
| Verified date, pathologist, scientist | `curation_verified_by` / `curation_verified_date` per classification | report-level sign-off: `FINAL` status carries `modified` + `user`; names are `case_fields` or LIS |
| Method versions: DRAGEN, VariantGrid, annotation sources | VCF header (`VCF.header`), `settings.VERSION`-style constant, `VariantAnnotationVersion` per build | have; `versions` in context |
| Genes tested | `EnrichmentKit` gene list | have when linked |

### Per variant

| Report / JSON field | Source | Gap |
|---|---|---|
| Gene | `gene_symbol` key; fusion partners from `GeneFusion` | have |
| Alteration (`var` / `amp` / `fusion`) | `ReportVariant.kind` | derived |
| Transcript, c., p. | `refseq_transcript_id`, `c_hgvs`, `p_hgvs` keys; both builds via `get_c_hgvs` | have |
| VAF | `allele_frequency` key (autopopulated from the sample genotype) | have |
| AMP tier incl. sub-tier | derived from `somatic:clinical_significance` + `amp:level_*` | §Tier |
| Copy number, fold change | `copy_number` key; `CohortGenotype.info["CN"]` / `FORMAT/SM` | have; fold change stays per observation |
| Report Y/N | `variant_reported` key → `CaseReportClassification.reported` | have, editable at build |
| Exon "5 of 11" | `exon` key | have as entered |
| COSMIC count, gnomAD | `cosmic_cnt`, `gnomad_af` keys | have |
| Variant narrative | `somatic:summary_interpretation` / `interpretation_summary` | have |
| Gene paragraph | `h_summary` by default | §Gene-level text |
| Fusion read support (`VAF: 29`) | `CohortGenotype` alt reads for the fusion row | have - surfaced as `vaf` for fusions |

The only new columns the report needs are the case-level ones on `CaseReport`; everything per variant already exists
as an evidence key or a variant fact.

---

## Rendering

All three renderers take a `ReportContext` and live in the new classification/report/ package.

- **HTML** - `engines['django'].from_string(template.case_template).render(context_dict)`, the way
  `ClassificationReport.get_template` renders the single-record template today. The template can read
  `v.evidence.<key>.formatted` / `.value` / `.note` for every evidence key, `:` in keys already turned to `_` by
  `row_data`. Stored in the `html` column of `CaseReport` and served as the preview.
- **PDF** - `xhtml2pdf.pisa.CreatePDF(html, dest=buffer)`. New dependency `xhtml2pdf` in `requirements.in`; it sits on
  `reportlab`, which is already there. The template's `<style>` carries the print CSS (`@page` size, margins, running
  header and footer, `page-break-before` on the Variant Interpretation section). Images, if any, are inline data URIs.
- **DOCX** - `html2docx(html, title).getvalue()`. New dependency `html2docx` (brings `python-docx`). Headings,
  paragraphs, tables, bold / italic / underline and lists carry across; the template is written to that vocabulary,
  which is all the example report uses.
- **JSON** - the Django engine renders `json_template` with `{% load js_tags %}` so values go through `jsonify`
  (`uicore/templatetags/js_tags.py`); the result is `json.loads`ed and stored in `json_output`. Blank template = the
  context dict itself.

All four render synchronously in the request; xhtml2pdf on a report of this size is well under a second. The build
form's template select has a `Preview` that renders the HTML without saving a `CaseReport`.

A build re-renders nothing it does not need to: `Rebuild documents` on a draft re-renders from `context_snapshot`
(a template fix), `New version` builds a new context from the current published versions.

### Template examples

Variant Interpretation section of `case_template` (Django template syntax; `vaf_percent`, `alteration` and
`gene_label` are convenience fields the context adds beside the raw values):

```
{% for tier in tier_groups %}
<h3 style="page-break-before: always">{{ tier.label }}</h3>
{% if not tier.genes %}<p>No variants detected.</p>{% endif %}
{% for gene in tier.genes %}
  <p>
  {% for v in gene.variants %}
    <b>{{ gene.gene_symbol }}.</b> {{ v.evidence.refseq_transcript_id.value }}:{{ v.evidence.c_hgvs.value }};
    {{ v.evidence.p_hgvs.value }}. Tier {{ v.amp_tier }}. VAF: {{ v.vaf_percent }}%
    {% if not forloop.last %}<br/>and<br/>{% endif %}
  {% endfor %}
  </p>
  {% for v in gene.variants %}<p>{{ v.evidence.somatic_summary_interpretation.value }}</p>{% endfor %}
  <p>{{ gene.gene_summary }}</p>
{% endfor %}
{% endfor %}
```

`json_template` producing the `tso500_data` shape:

```
{% load js_tags %}
{
  "Seq Run": {{ sequencing_runs.0|jsonify }},
  "Patient": {{ patient.patient_code|jsonify }},
  "Report Values": {
    "Panel": {{ case_values.panel|jsonify }},
    "TMB": {{ measures.tmb.value|stringformat:"s"|jsonify }},
    "MSI": {"Category": {{ measures.msi.call|jsonify }}},
    "Purity": {{ measures.tumour_fraction.value|jsonify }},
    "Assay Success": {{ case_values.assay_success|jsonify }},
    "Caveats": {{ case_values.caveats|jsonify }},
    "Nata Report": {"ID": {{ case_report.external_report_id|jsonify }}, "Date": {{ case_report.report_date|jsonify }}},
    "Comments": {"Summary Interpretation": {{ summary|jsonify }}, "Mutations Comment": {{ case_values.mutations_comment|jsonify }}}
  },
  "Variants": [
  {% for v in variants %}
    {"Alteration": {{ v.alteration|jsonify }},
     "Gene": {{ v.gene_label|jsonify }},
     "Description": {{ v.evidence.p_hgvs.value|jsonify }},
     "transcript": {{ v.evidence.refseq_transcript_id.value|jsonify }},
     "c": {{ v.evidence.c_hgvs.value|jsonify }},
     "VAF": {{ v.vaf|jsonify }},
     "CopyNo": {{ v.copy_number|jsonify }},
     "AMPTier": "Tier {{ v.amp_tier }}",
     "Report": "{% if v.reported %}Y{% else %}N{% endif %}"}{% if not forloop.last %},{% endif %}
  {% endfor %}
  ]
}
```

`case_values.assay_success` and `case_values.caveats` are `case_fields` of type `bool`, grouped by a `group` key on
the field definition so the form shows them as two rows of checkboxes and the context exposes each group as a dict.

The SA Path versions of both templates live in the sapath repo beside `sapath/default_report.html`, loaded into the
`ClassificationReportTemplate` row by a `ManualOperation` migration there.

---

## UI

All on the Classify & Report tab (`analysis/templates/analysis/classify_report_tab.html`), which already has the
tick-boxes and template select.

1. **Case types**: `ClassifyReportCase` (`analysis/classify_report.py:ClassifyReportCase`) gains `for_specimen` and
   `for_extraction`, resolving samples through `get_sample_group`, and `_get_case` in the views accepts
   `specimen` / `extraction`. The specimen and extraction pages get the tab, so a TSO 500 case (DNA + RNA arms) is
   reported from the specimen page, which is where its measures already show.
2. **Build report** replaces the current "Create report" POST-to-HTML with a modal: template (filtered by bucket),
   the template's `case_fields`, `summary` (prefilled from the latest draft for the case), and the ticked
   classifications listed in report order with their `amp_tier`, any `warnings`, and a Report Y/N toggle. Submitting
   creates the `CaseReport` and its rows, renders HTML, PDF, DOCX and JSON, and the modal shows the preview and the
   three downloads.
3. **Reports card** at the foot of the tab: every `CaseReport` for the case, newest first - status, template, who,
   when, the pinned classifications, downloads, and per row `Finalise`, `New version`, `Rebuild documents` (drafts),
   and the `external_report_id` / `report_date` fields to fill in once the LIS has issued them.
4. The existing HTML multi-variant view (`multi_classification_report`) stays as a preview for templates that only
   have an HTML `template`.

## Permissions

A `CaseReport` is visible to its lab's users and to anyone who can see the case (`Patient.filter_for_user` and the
sample's VCF permissions, the same rules the tab applies); only the lab's users build, finalise or delete. Every
classification pinned must be visible to the builder at build time, and the `CaseReportClassification` rows are the
record of what was shown.

---

## Order of work

1. **Models + migration + admin.** The template columns, `CaseReport`, `CaseReportClassification`; the admin gets
   `case_template` and `json_template` editors with the fixture-render check. `scripts/vg map` after.
2. **Context builder.** `ReportVariant` / `ReportContext`, the ordering, `amp_tier`, kinds, measures, versions.
   Refactor `ClassificationReport.context()` onto it so `record` / `classifications` / `gene_groups` keep their shape.
   Pure functions, tested over fake modifications.
3. **Renderers.** HTML, PDF (xhtml2pdf), DOCX (html2docx), JSON; `xhtml2pdf` and `html2docx` added to
   `requirements.in` and the lock recompiled with uv. A minimal generic case template and blank JSON template under a
   new classification/test_data/ directory for tests and as the shipped default.
4. **Tab.** Case types, the build modal, the reports card, download views, specimen and extraction page tabs.
   `vg page` on the sample and specimen pages before and after.
5. **SA Path templates** in the sapath repo: the TSO 500 `case_template` mirroring the `.txt` example, the `json_template`
   above, `case_fields` for the assay-success and caveat flags, panel and mutations comment; `allele_origin_bucket`
   somatic. Checked against the two example JSONs.
6. **Finalise.** Status transitions, supersede, and stamping the report onto each pinned classification:
   `report_date` (the finalise date), `variant_reported` (from `reported`) and `report_id` (when
   `external_report_id` is known; entering it later on a final report stamps that key alone). Each stamp is one
   `patch_value` with `save=True` followed by `publish_latest` at the record's current share level.
   `patch_value` writes nothing when the values already match, so re-finalising, rebuilding or re-entering the same
   LIS id creates no further modifications. A classification whose last edited version is newer than its published
   one is left alone and listed on the finalise dialog - publishing it would push someone's unsubmitted edits out
   under the report's name. After a successful stamp the `CaseReportClassification` re-points to the new published
   version, which differs from the rendered one only in those three keys, so the tab's stale check stays quiet and
   the pinned version is the one Shariant receives.

## Tests that earn their keep

- Ordering: kind → sub-tier → gene → VAF, with a fusion filed under its 5' partner and an amplification under Tier II
  beside small variants.
- `amp_tier` for every row of the table above, including the warning rows.
- `json_template` blank → canonical dump; a template that renders invalid JSON is rejected on save.
- The test template renders to HTML containing each included c.HGVS once, the PDF and DOCX derive from it without
  error, and a `reported=False` variant is absent from the tier text and present in the JSON with `"Report": "N"`.
- Specimen and extraction cases resolve to the same sample set as `get_sample_group`.
- A user outside the lab cannot download or finalise.

## Decisions

- **Sub-tier is derived**, from `somatic:clinical_significance` + `amp:level_*`, the schema settled for Shariant.
- **Gene paragraph is `h_summary`**, most recent wins where they differ, with a warning and links to fix.
- **Finalise stamps `report_date`, `variant_reported` and `report_id`** onto the pinned classifications, idempotently
  and only where nothing is unsubmitted (§Order of work 6).
- **`"Patient"` in the JSON is `Patient.patient_code`**, the Omico ID.
- **Somatic only for now**; the models leave room for a germline case template later.
