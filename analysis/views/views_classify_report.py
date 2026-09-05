"""
The Classify & Report tab shown on the sample and patient pages.

It lives in the analysis app because it is built around VariantTag - analysis already depends on classification,
and classification must not depend on analysis.

@see analysis/classify_report.py for the queue itself.
"""
from django.http import Http404
from django.http.response import HttpResponseBase
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST

from analysis.classify_report import ClassifyReportCase, outstanding_tag_count, tag_summary
from analysis.models import VariantTag
from analysis.models.nodes.analysis_node import AnalysisClassification
from analysis.variant_tag_operations import (
    get_sample_genotype_for_variant_tag,
    retire_requires_classification_tags_for_samples,
)
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from classification.models import Classification, ClassificationReportTemplate
from classification.views.classification_export_report import ClassificationReport
from classification.views.views import classification_created_response, create_classification_object
from patients.models import Patient
from patients.models_enums import Zygosity
from snpdb.forms import UserLabChoiceForm
from snpdb.models import Sample, UserSettings

CASE_TYPE_SAMPLE = "sample"
CASE_TYPE_PATIENT = "patient"


def _get_case(user, case_type: str, case_id: int) -> ClassifyReportCase:
    if case_type == CASE_TYPE_SAMPLE:
        return ClassifyReportCase.for_sample(user, Sample.get_for_user(user, case_id))
    if case_type == CASE_TYPE_PATIENT:
        return ClassifyReportCase.for_patient(user, Patient.get_for_user(user, case_id))
    raise Http404(f"Unknown case type '{case_type}'")


def _classify_report_context(case: ClassifyReportCase, case_type: str, case_id: int) -> dict:
    rows = case.queue_rows()
    return {
        "case": case,
        "case_type": case_type,
        "case_id": case_id,
        "rows": rows,
        "outstanding": outstanding_tag_count(rows),
        "tag_summary": tag_summary(rows),
        "classification_modifications": list(case.classification_modifications()),
        "report_templates": ClassificationReportTemplate.objects.exclude(template="").order_by("name"),
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
    retire_requires_classification_tags_for_samples(classification, case.samples, request.user)
    return classification_created_response(request, classification)


@require_POST
def multi_classification_report(request, case_type: str, case_id: int):
    """ One report over several of the case's classifications, grouped by gene """
    case = _get_case(request.user, case_type, case_id)
    selected_ids = request.POST.getlist("classification_modification_id")
    modifications = [cm for cm in case.classification_modifications() if str(cm.pk) in selected_ids]
    if not modifications:
        raise Http404("No classifications selected for this report")

    report_template = get_object_or_404(ClassificationReportTemplate, pk=request.POST["report_template"])
    # The header a multi-variant report opens with - who the case is, and what was sequenced
    extra_context = {
        "case": case.obj,
        "case_label": str(case.obj),
        "samples": case.samples,
        "patient": case.obj if case_type == CASE_TYPE_PATIENT else case.obj.patient,
    }
    return ClassificationReport(modifications[0], user=request.user, classifications=modifications,
                                report_template=report_template, extra_context=extra_context).serve()
