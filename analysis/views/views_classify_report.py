"""
The Classify & Report tab shown on the sample, patient, specimen and extraction pages.

It lives in the analysis app because it is built around VariantTag - analysis already depends on classification,
and classification must not depend on analysis.

The tab has two halves: the tag queue (@see analysis/classify_report.py) and the case report - the Build report
modal, which pins the ticked classifications and renders the four documents, and the Reports card listing what
the case has had built. A built report's own actions - the downloads, Finalise, New version, Rebuild documents
and the LIS fields - are classification/views/views_case_report.py, since they need the report, not the case.
"""
from typing import Optional

from django.core.exceptions import PermissionDenied
from django.http import Http404, HttpResponse
from django.http.response import HttpResponseBase, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST

from analysis.classify_report import (
    ClassifyReportCase,
    case_allele_origin_bucket,
    outstanding_tag_count,
    tag_summary,
)
from analysis.models import VariantTag
from analysis.models.nodes.analysis_node import AnalysisClassification
from analysis.variant_tag_operations import (
    get_sample_genotype_for_variant_tag,
    resolve_requires_classification_tags_for_samples,
    resolve_variant_tag,
)
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationModification,
    ClassificationReportTemplate,
)
from classification.report.case_report_builder import build_case_report, preview_case_report_html
from classification.report.case_report_context import build_report_variants
from classification.views.classification_export_report import ClassificationReport
from classification.views.views import classification_created_response, create_classification_object
from library.utils import first
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import Zygosity
from snpdb.forms import UserLabChoiceForm
from snpdb.models import Lab, Sample, UserSettings

CASE_TYPE_SAMPLE = "sample"
CASE_TYPE_PATIENT = "patient"
CASE_TYPE_SPECIMEN = "specimen"
CASE_TYPE_EXTRACTION = "extraction"


def _get_case(user, case_type: str, case_id: int) -> ClassifyReportCase:
    if case_type == CASE_TYPE_SAMPLE:
        return ClassifyReportCase.for_sample(user, Sample.get_for_user(user, case_id))
    if case_type == CASE_TYPE_PATIENT:
        return ClassifyReportCase.for_patient(user, Patient.get_for_user(user, case_id))
    if case_type == CASE_TYPE_SPECIMEN:
        return ClassifyReportCase.for_specimen(user, Specimen.get_for_user(user, case_id))
    if case_type == CASE_TYPE_EXTRACTION:
        return ClassifyReportCase.for_extraction(user, Extraction.get_for_user(user, case_id))
    raise Http404(f"Unknown case type '{case_type}'")


def _classify_report_context(case: ClassifyReportCase, case_type: str, case_id: int) -> dict:
    rows = case.queue_rows()
    modifications = list(case.classification_modifications())
    return {
        "case": case,
        "case_type": case_type,
        "case_id": case_id,
        "rows": rows,
        "outstanding": outstanding_tag_count(rows),
        "unclassified": sum(1 for row in rows if row.needs_classification),
        "tag_summary": tag_summary(rows),
        "classification_modifications": modifications,
        "report_templates": ClassificationReportTemplate.objects.exclude(template="").order_by("name"),
        "case_report_templates": ClassificationReportTemplate.case_templates_for_bucket(
            case_allele_origin_bucket(modifications)),
        "case_reports": list(case.case_reports()),
        "can_create_classifications": Classification.can_create_via_web_form(case.user),
    }


def classify_report_tab(request, case_type: str, case_id: int):
    case = _get_case(request.user, case_type, case_id)
    return render(request, 'analysis/classify_report_tab.html',
                  _classify_report_context(case, case_type, case_id))


def sample_classify_report_tab(request, sample_id: int):
    return classify_report_tab(request, CASE_TYPE_SAMPLE, sample_id)


def patient_classify_report_tab(request, patient_id: int):
    return classify_report_tab(request, CASE_TYPE_PATIENT, patient_id)


def specimen_classify_report_tab(request, specimen_id: int):
    return classify_report_tab(request, CASE_TYPE_SPECIMEN, specimen_id)


def extraction_classify_report_tab(request, extraction_id: int):
    return classify_report_tab(request, CASE_TYPE_EXTRACTION, extraction_id)


def classify_report_summary(request, case_type: str, case_id: int) -> JsonResponse:
    """ The counts the sample / patient page hangs off the Classify & Report tab label, so the page says
        whether there is anything to do before the tab is opened. Fetched after the page renders - working
        out which taggings are the case's walks every analysis its samples are in """
    case = _get_case(request.user, case_type, case_id)
    outstanding = sum(1 for variant_tag, _ in case.variant_tags() if not variant_tag.is_resolved)
    return JsonResponse({"outstanding": outstanding,
                         "classifications": case.classification_modifications().count()})


def classify_report_tag_dialog(request, case_type: str, case_id: int, variant_tag_id: int):
    """ The launcher for one tagged variant - the previous classifications of the allele, and where the new
        classification is going. Which (if any) previous classification applies is always the scientist's call """
    case = _get_case(request.user, case_type, case_id)
    variant_tag = VariantTag.get_for_user(request.user, variant_tag_id)
    row = case.queue_row(variant_tag)

    genome_build = variant_tag.genome_build
    vts = VariantTranscriptSelections(variant_tag.variant, genome_build)
    selected_transcript = next((td for td in vts.transcript_data if td.get("selected")), {})

    lab, lab_error = UserSettings.get_lab_and_error(request.user)
    lab_form = UserLabChoiceForm(user=request.user, default_lab=lab) if lab else None

    sample_genotype = get_sample_genotype_for_variant_tag(row.sample, variant_tag) if row.sample else None

    # Gene content is what's left to reuse when the lab has never seen this variant before - same
    # deduplicated rows as the create-from-variant page, so the two can't drift
    gene_groups = []
    if not row.copyable and row.gene_symbol:
        gene_groups = ClassificationConsensus.gene_consensus_groups(
            gene_symbol=row.gene_symbol, user=request.user, allele_origin_bucket=row.allele_origin_bucket)

    context = {
        "case": case,
        "case_type": case_type,
        "case_id": case_id,
        "row": row,
        "variant_tag": variant_tag,
        "genome_build": genome_build,
        "refseq_transcript_accession": selected_transcript.get(VariantTranscriptSelections.REFSEQ_TRANSCRIPT) or "",
        "ensembl_transcript_accession": selected_transcript.get(VariantTranscriptSelections.ENSEMBL_TRANSCRIPT) or "",
        "lab": lab,
        "lab_error": lab_error,
        "lab_form": lab_form,
        "sample_genotype": sample_genotype,
        "zygosity_display": Zygosity.display(sample_genotype.zygosity) if sample_genotype else None,
        "gene_groups": gene_groups,
    }
    return render(request, 'analysis/classify_report_tag_dialog.html', context)


@require_POST
def create_classification_for_case(request, case_type: str, case_id: int, variant_tag_id: int) -> HttpResponseBase:
    """ Create the classification a queue row is asking for - with a previous record's consensus copied in when
        the scientist picked one (copy_from_vcm_id), otherwise empty for them to fill in on the full form """
    case = _get_case(request.user, case_type, case_id)
    variant_tag = VariantTag.get_for_user(request.user, variant_tag_id)

    classification = create_classification_object(request)
    if analysis := variant_tag.analysis:
        if analysis.can_write(request.user):
            AnalysisClassification.objects.create(analysis=analysis, classification=classification)
    resolved = resolve_requires_classification_tags_for_samples(classification, case.samples, request.user)
    # The queue row shows the link either way - unresolved it also offers the button that says this is the person
    return classification_created_response(request, classification,
                                           {"resolved": variant_tag.pk in {vt.pk for vt in resolved}})


@require_POST
def resolve_variant_tag_for_case(request, case_type: str, case_id: int, variant_tag_id: int) -> HttpResponseBase:
    """ The scientist saying the case's classification is what this tagging was asking for - needed when the
        tagging never knew whose it was, so it couldn't be resolved automatically """
    case = _get_case(request.user, case_type, case_id)
    variant_tag = VariantTag.get_for_user(request.user, variant_tag_id)
    if not variant_tag.can_write(request.user):
        raise PermissionDenied(f"You have read-only access to VariantTag {variant_tag_id}")

    row = case.queue_row(variant_tag)
    if row.classification is None:
        raise Http404("This case has no classification for the tagged variant")
    resolve_variant_tag(variant_tag, row.classification, request.user)
    return JsonResponse({"resolved": True})


@require_POST
def multi_classification_report(request, case_type: str, case_id: int):
    """ One report over several of the case's classifications, grouped by gene """
    case = _get_case(request.user, case_type, case_id)
    selected_ids = request.POST.getlist("classification_modification_id")
    modifications = [cm for cm in case.classification_modifications() if str(cm.pk) in selected_ids]
    if not modifications:
        raise Http404("No classifications selected for this report")

    report_template = get_object_or_404(ClassificationReportTemplate, pk=request.POST["report_template"])
    # The header a multi-variant report opens with - who the case is, and what was sequenced.
    # Every level resolves its patient the same way, so an extraction case has one too
    extra_context = {
        "case": case.obj,
        "case_label": str(case.obj),
        "samples": case.samples,
        "patient": first(case.patients),
    }
    return ClassificationReport(modifications[0], user=request.user, classifications=modifications,
                                report_template=report_template, extra_context=extra_context).serve()


def _selected_modifications(case: ClassifyReportCase, request) -> list[ClassificationModification]:
    """ The ticked rows, resolved against what the case actually has - every classification a report
        pins has to be one this user can see, and the CaseReportClassification rows are that record """
    selected_ids = set(request.POST.getlist("classification_modification_id"))
    modifications = [cm for cm in case.classification_modifications() if str(cm.pk) in selected_ids]
    if not modifications:
        raise Http404("No classifications selected for this report")
    return modifications


def _case_report_template(request, modifications: list[ClassificationModification]) \
        -> Optional[ClassificationReportTemplate]:
    """ The template the form named, or the first the case is offered - only ever one of the bucket's,
        so a somatic template can't be picked for a germline case by posting its name """
    templates = ClassificationReportTemplate.case_templates_for_bucket(
        case_allele_origin_bucket(modifications))
    if name := request.POST.get("report_template"):
        return get_object_or_404(templates, pk=name)
    return templates.first()


def _reported_by_pk(request, modifications: list[ClassificationModification]) -> dict[int, bool]:
    """ The form's Report Y/N toggles, which start from the variant_reported evidence key """
    return {cm.pk: f"reported_{cm.pk}" in request.POST for cm in modifications}


def _flat_case_values(case_values: Optional[dict]) -> dict:
    """ A draft's answers keyed the way the form names them - the context holds a group as a dict """
    flat = {}
    for key, value in (case_values or {}).items():
        if isinstance(value, dict):
            flat.update(value)
        else:
            flat[key] = value
    return flat


def _report_lab(request, user) -> Lab:
    """ Whose report it is - the report is the lab's document, so the lab decides who may finalise it """
    if lab_id := request.POST.get("lab"):
        return get_object_or_404(Lab.valid_labs_qs(user), pk=lab_id)
    lab, lab_error = UserSettings.get_lab_and_error(user)
    if lab is None:
        raise PermissionDenied(lab_error or "You are not a member of a lab, so cannot build a report")
    return lab


@require_POST
def case_report_build_dialog(request, case_type: str, case_id: int):
    """ The Build report form: the template (filtered to the case's bucket), the template's own
        case_fields, the case level summary, and the ticked classifications in report order with
        their AMP sub-tier, any warnings and a Report Y/N toggle - so the order, the tiers and the
        document's own inputs are all checked before anything is built. Re-posted when the template
        changes, since a different template asks for different case_fields """
    case = _get_case(request.user, case_type, case_id)
    modifications = _selected_modifications(case, request)
    template = _case_report_template(request, modifications)

    draft = case.latest_draft_report()
    lab, lab_error = UserSettings.get_lab_and_error(request.user)
    context = {
        "case": case,
        "case_type": case_type,
        "case_id": case_id,
        "template": template,
        "case_report_templates": ClassificationReportTemplate.case_templates_for_bucket(
            case_allele_origin_bucket(modifications)),
        "variants": build_report_variants(modifications, request.user),
        "summary": draft.summary if draft else "",
        "case_values": _flat_case_values(draft.case_values if draft else None),
        "lab": lab,
        "lab_error": lab_error,
        "lab_form": UserLabChoiceForm(user=request.user, default_lab=lab) if lab else None,
    }
    return render(request, 'analysis/case_report_build_dialog.html', context)


@require_POST
def create_case_report(request, case_type: str, case_id: int):
    """ Build the case report the form describes, and show the preview and the downloads """
    case = _get_case(request.user, case_type, case_id)
    modifications = _selected_modifications(case, request)
    template = _case_report_template(request, modifications)
    if template is None:
        raise Http404("No case report template is configured for this case")

    case_report = build_case_report(
        request.user, template, _report_lab(request, request.user), case.source_level, case.obj,
        modifications, reported_by_pk=_reported_by_pk(request, modifications),
        summary=request.POST.get("summary", ""),
        case_values=template.case_values_from_form(request.POST))
    return render(request, 'analysis/case_report_built.html',
                  {"case_report": case_report, "case_type": case_type, "case_id": case_id})


@require_POST
def preview_case_report(request, case_type: str, case_id: int) -> HttpResponse:
    """ The document as it would read, rendered without saving a CaseReport """
    case = _get_case(request.user, case_type, case_id)
    modifications = _selected_modifications(case, request)
    template = _case_report_template(request, modifications)
    if template is None:
        raise Http404("No case report template is configured for this case")

    html = preview_case_report_html(
        request.user, template, case.source_level, case.obj, modifications,
        reported_by_pk=_reported_by_pk(request, modifications),
        summary=request.POST.get("summary", ""),
        case_values=template.case_values_from_form(request.POST))
    return HttpResponse(content=html, content_type="text/html")
