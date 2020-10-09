from datetime import datetime

import rest_framework
from crispy_forms.bootstrap import FieldWithButtons, StrictButton
from crispy_forms.layout import Layout, Field, Button, Submit
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.forms import formset_factory
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.views.decorators.http import require_POST
from global_login_required import login_not_required
from jfu.http import upload_receive, UploadResponse, JFUResponse
import json
import mimetypes
from requests.models import Response
from typing import Optional

from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from flags.models.models import FlagType
from genes.forms import GeneSymbolForm
from genes.hgvs import get_kind_and_transcript_accession_from_invalid_hgvs
from library.file_utils import rm_if_exists
from library.guardian_utils import is_superuser
from library.log_utils import log_traceback, report_event
from snpdb.forms import SampleChoiceForm, UserSelectForm, LabSelectForm
from snpdb.models import Variant, UserSettings, Sample, Allele, Lab
from snpdb.models.models_genome import GenomeBuild
from uicore.utils.form_helpers import form_helper_horizontal
from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import \
    create_classification_for_sample_and_variant_objects, generate_auto_populate_data
from classification.classification_stats import get_grouped_classification_counts, \
    get_classification_counts, get_criteria_counts
from classification.enums import SubmissionSource, SpecialEKeys
from classification.forms import EvidenceKeyForm
from classification.models import ClassificationAttachment, Classification, \
    ClassificationRef, ClassificationJsonParams, ClassificationConsensus
from classification.models.clinical_context_models import ClinicalContext
from classification.models.evidence_key import EvidenceKeyMap
from classification.models.flag_types import classification_flag_types
from classification.models.classification import ClassificationModification
from classification.classification_changes import ClassificationChanges
from classification.views.classification_datatables import ClassificationDatatableConfig
from classification.views.classification_export_csv import ExportFormatterCSV
from classification.views.classification_export_redcap import ExportFormatterRedcap
from variantopedia.forms import SearchForm, SearchAndClassifyForm


@user_passes_test(is_superuser)
def activity(request, latest_timestamp: Optional[str] = None):
    if latest_timestamp:
        latest_timestamp = datetime.fromtimestamp( float(latest_timestamp) )
    changes = ClassificationChanges.list_changes(latest_date=latest_timestamp)
    last_date = changes[len(changes)-1].date.timestamp() if changes else None
    context = {
        'changes': changes,
        'last_date': last_date,
        'can_create_classifications': settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN
    }
    return render(request, 'classification/activity.html', context)


def classifications(request):
    user_settings = UserSettings.get_for_user(request.user)

    initial = {'classify': True}
    search_and_classify_form = SearchAndClassifyForm(initial=initial)
    search_and_classify_form.fields['search'].label = "HGVS / dbSNP / VCF coordinate"
    helper = form_helper_horizontal()
    helper.layout = Layout(
        FieldWithButtons(Field('search', placeholder=""), Submit(name="action", value="Go"))
    )
    search_and_classify_form.helper = helper

    flag_types = FlagType.objects.filter(context=classification_flag_types.classification_flag_context)\
        .exclude(pk__in=[
            'classification_withdrawn',
            'classification_submitted',
            'classification_clinical_context_change']).order_by('label')
    flag_type_json = []
    for ft in flag_types:
        flag_type_json.append({'id': ft.pk, 'label': ft.label, 'description': ft.description})

    context = {
        "can_create_classification": Classification.can_create_via_web_form(request.user),
        "flag_types": flag_type_json,
        "gene_form": GeneSymbolForm(),
        "user_form": UserSelectForm(),
        "lab_form": LabSelectForm(),
        "labs": Lab.valid_labs_qs(request.user),
        "search_and_classify_form": search_and_classify_form,
        "genome_build": user_settings.default_genome_build,
        "VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME": settings.VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME,
        "VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN": settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN,
        "VARIANT_CLASSIFICATION_ID_FILTER": settings.VARIANT_CLASSIFICATION_ID_FILTER,
        "datatable_config": ClassificationDatatableConfig(request),
        "user_settings": user_settings,
    }
    return render(request, 'classification/classifications.html', context)


class AutopopulateView(APIView):

    def get(self, request):
        variant_id = request.GET.get("variant_id")
        genome_build_name = request.GET.get("genome_build_name")
        refseq_transcript_accession = request.GET.get("refseq_transcript_accession")
        ensembl_transcript_accession = request.GET.get("ensembl_transcript_accession")
        sample_id = request.GET.get("sample_id")
        copy_from_id = request.GET.get("copy_from_vcm_id")

        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        variant: Variant = get_object_or_404(Variant, pk=variant_id)
        sample: Optional[Sample] = None

        if sample_id:
            sample = Sample.get_for_user(request.user, sample_id)

        auto_data = generate_auto_populate_data(
            variant=variant,
            genome_build=genome_build,
            refseq_transcript_accession=refseq_transcript_accession,
            ensembl_transcript_accession=ensembl_transcript_accession,
            sample=sample,
            annotation_version=None
        )

        used_keys = set()
        used_keys.add(SpecialEKeys.AUTOPOPULATE)
        complete_values = list()

        for data_subset in auto_data.flatten():
            name = data_subset.name
            for key, value in data_subset.data_direct.items():
                if key not in used_keys:
                    used_keys.add(key)
                    if not isinstance(value, dict):
                        value = {'value': value}
                    complete_values.append({'key': key, 'blob': value, 'source': name})

        if copy_from_id:
            copy_from = ClassificationModification.objects.get(pk=copy_from_id)
            consensus_patch = ClassificationConsensus(variant, request.user, copy_from=copy_from).consensus_patch
            for key, blob in consensus_patch.items():
                if key not in used_keys:
                    used_keys.add(key)
                    complete_values.append({'key': key, 'blob': blob, 'source': 'copy from latest'})

        key_to_order = dict()
        index = 1
        for ekey in EvidenceKeyMap().all_keys:
            key_to_order[ekey.key] = index
            index = index + 1

        def key_order(tests_key: str):
            nonlocal key_order
            return key_to_order.get(tests_key, 0) * 2

        complete_values = sorted(complete_values, key=lambda entry: key_order(entry['key']))
        return rest_framework.response.Response(status=HTTP_200_OK, data={'data': complete_values})

@require_POST
def create_classification(request):
    if not Classification.can_create_via_web_form(request.user):
        raise PermissionDenied('User cannot create classifications via web form')

    variant_id = request.POST.get("variant_id")
    genome_build_name = request.POST["genome_build_name"]
    refseq_transcript_accession = request.POST.get("refseq_transcript_accession")
    ensembl_transcript_accession = request.POST.get("ensembl_transcript_accession")
    sample_id = request.POST.get("sample_id")
    copy_from_id = request.POST.get("copy_from_vcm_id")
    evidence_json = request.POST.get("evidence_json")

    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    if variant_id:
        variant = get_object_or_404(Variant, pk=variant_id)
    else:
        Classification.check_can_create_no_classification_via_web_form(request.user)
        variant = None

    if sample_id:
        sample = Sample.get_for_user(request.user, sample_id)
    else:
        sample = None

    classification = create_classification_for_sample_and_variant_objects(request.user, sample,
                                                                                          variant, genome_build,
                                                                                          refseq_transcript_accession=refseq_transcript_accession,
                                                                                          ensembl_transcript_accession=ensembl_transcript_accession)
    if evidence_json:
        evidence = json.loads(evidence_json)
        classification.patch_value(
            patch=evidence,
            clear_all_fields=False,
            user=request.user,
            source=SubmissionSource.FORM,
            leave_existing_values=True,
            save=True,
            make_patch_fields_immutable=False)

    classification.publish_latest(request.user)

    if copy_from_id:
        copy_from = ClassificationModification.objects.get(pk=copy_from_id)
        consensus_patch = ClassificationConsensus(variant, request.user, copy_from=copy_from).consensus_patch
        classification.patch_value(
            patch=consensus_patch,
            clear_all_fields=False,
            user=request.user,
            source=SubmissionSource.CONSENSUS,
            leave_existing_values=True,
            save=True,
            make_patch_fields_immutable=False)
        classification.publish_latest(request.user)

    return redirect(classification)


def classification_history(request, record_id):
    ref = ClassificationRef.init_from_str(request.user, str(record_id))
    ref.check_exists()
    #  changes page will show things before they were submitted and not hide values
    #  so need to restrict to people who've always had access to the record
    ref.check_security(must_be_writable=True)

    changes = ClassificationChanges.list_changes(classification=ref.record, limit=100)
    context = {
        'changes': changes,
        'can_create_classifications': Classification.can_create_via_web_form(request.user)
    }
    return render(request, 'classification/activity.html', context)


def view_classification(request, record_id):
    ref = ClassificationRef.init_from_str(request.user, str(record_id))
    ref.check_exists()
    ref.check_security()

    existing_files = get_classification_attachment_file_dicts(ref.record)
    other_classifications_summary = ref.record.get_other_classifications_summary_for_variant(request.user)

    record = ref.as_json(ClassificationJsonParams(current_user=request.user,
                                                         include_data=True,
                                                         include_lab_config=True))

    # default to the natural build of the classification
    genome_build = None
    variant = ref.record.variant
    if variant:
        genome_build = variant.genome_build
    if genome_build is None:
        genome_build = GenomeBuild.default_build()

    vc: Classification = ref.record

    context = {'vc': vc,
               'record': record,
               'genome_build': genome_build,
               'asterix_view': settings.VARIANT_CLASSIFICAITON_DEFAULT_ASTERIX_VIEW,
               'existing_files': existing_files,
               'other_classifications_summary': other_classifications_summary,
               'report_enabled': not not vc.lab.organization.classification_report_template,
               'attachments_enabled': settings.VARIANT_CLASSIFICATION_FILE_ATTACHMENTS
               }
    return render(request, 'classification/classification.html', context)


def view_classification_diff(request):
    records = []
    extra_data = dict()

    if request.GET.get('history'):
        vc_id = int(request.GET.get('history'))
        vc = Classification.objects.get(pk=vc_id)
        qs = ClassificationModification.objects.filter(classification=vc, published=True)
        qs = ClassificationModification.filter_for_user(request.user, qs).order_by('-created')
        records = list(qs.all())
        if request.GET.get('debug') != 'true' and len(records) > 1:
            consider = records
            records = []
            last_record: Optional[ClassificationModification] = None
            for record in consider:
                if last_record is None or not last_record.is_significantly_equal(record):
                    last_record = record
                    records.append(record)
                else:
                    extra_data[last_record.id] = {'first_seen': record.created.timestamp()}
                    # provide an earlier date value so if we updated the record with a minimal change 2 months after
                    # the original, we still acknowledge the fact that it was first uploaded 2 months ago

        if vc.can_write(request.user) and vc.last_published_version.id != vc.last_edited_version.id:
            records.insert(0, vc.last_edited_version)

    elif request.GET.get('clinical_context'):
        cc = ClinicalContext.objects.get(pk=request.GET.get('clinical_context'))
        records = cc.classification_modifications

    elif request.GET.get('variant_compare'):
        ref = ClassificationRef.init_from_str(user=request.user, id_str=request.GET.get('variant_compare'))
        compare_all = ClassificationModification.latest_for_user(user=request.user, allele=ref.record.variant.allele, published=True)
        # filter out variant we're comparing with, make it last in calculation
        latest_others_for_variant = [vcm for vcm in compare_all if vcm.classification.id != ref.record.id]
        compare_all = [ref.modification] + latest_others_for_variant
        records = compare_all

    elif request.GET.get('variant'):
        variant_id = int(request.GET.get('variant'))
        compare_all = ClassificationModification.latest_for_user(user=request.user, allele=Variant.objects.get(pk=variant_id).allele, published=True)
        records = compare_all

    elif request.GET.get('allele'):
        allele_id = int(request.GET.get('allele'))
        compare_all = ClassificationModification.latest_for_user(user=request.user, allele=Allele.objects.get(pk=allele_id), published=True)
        records = compare_all

    # filter out any Nones inserted by filtering on user permission etc
    records = [record for record in records if record]

    records_json = [vcm.as_json(ClassificationJsonParams(
        current_user=request.user,
        include_data=True,
        hardcode_extra_data=extra_data.get(vcm.id),
        fix_data_types=True
    )) for vcm in records]

    context = {
        'records_json': records_json
    }
    return render(request, 'classification/classification_diff.html', context)


def classification_qs(request):
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="classifications.csv"'

    extra_filters = request.GET.get("extra_filters")
    if extra_filters:
        extra_filters = json.loads(extra_filters)

    config = ClassificationDatatableConfig(request)
    qs = ClassificationModification.latest_for_user(user=request.user, published=True)
    qs = config.filter_queryset(qs)
    # TODO sort

    return qs


def classification_import_tool(request: HttpRequest) -> Response:
    EKeyFormSet = formset_factory(EvidenceKeyForm, extra=3)
    context = {
        "labs": Lab.valid_labs_qs(request.user)
    }
    return render(request, 'classification/classification_import_tool.html', context)


def export_classifications_grid(request):
    """
    CSV export of what is currently filtered into the classification grid
    """
    genome_build = UserSettings.get_for_user(request.user).default_genome_build
    qs = classification_qs(request)
    report_event(
        name='variant classification download',
        request=request,
        extra_data={
            'format': 'csv',
            'refer': 'classification listing',
            'approx_count': qs.count()
        }
    )
    return ExportFormatterCSV(user=request.user, genome_build=genome_build, qs=qs).export()


def export_classifications_grid_redcap(request):
    genome_build = UserSettings.get_for_user(request.user).default_genome_build
    qs = classification_qs(request)
    report_event(
        name='variant classification download',
        request=request,
        extra_data={
            'format': 'redcap',
            'refer': 'classification listing',
            'approx_count': qs.count()
        }
    )
    return ExportFormatterRedcap(user=request.user, genome_build=genome_build, qs=qs).export()


@require_POST
def classification_file_upload(request, classification_id):
    try:
        classification = get_object_or_404(Classification, pk=classification_id)
        if not classification.can_write(request.user):
            raise Exception('User can not edit this variant classification')
        uploaded_file = upload_receive(request)

        vc_attachment = ClassificationAttachment(classification=classification,
                                                        file=uploaded_file)
        vc_attachment.save()

        file_dict = vc_attachment.get_file_dict()
    except Exception as e:
        log_traceback()
        file_dict = {"error": str(e)}

    return UploadResponse(request, file_dict)


@require_POST
def classification_file_delete(request, pk):
    success = False
    try:
        vc_attachment = ClassificationAttachment.objects.get(pk=pk)
        if not vc_attachment.classification.can_write(request.user):
            raise Exception('User can not edit this variant classification')

        rm_if_exists(vc_attachment.file.path)
        if vc_attachment.thumbnail_path:
            rm_if_exists(vc_attachment.thumbnail_path)
        vc_attachment.delete()
        success = True
    except ClassificationAttachment.DoesNotExist:
        pass

    return JFUResponse(request, success)


def get_classification_attachment_file_dicts(classification):
    file_dicts = []
    for vca in classification.classificationattachment_set.all():
        file_dicts.append(vca.get_file_dict())

    file_dicts = list(reversed(file_dicts))  # JFU adds most recent at the end
    return file_dicts


def view_classification_file_attachment(request, pk, thumbnail=False):
    """ This is not done via static files, so we add security later """

    # TODO: Check security/access to Classification
    vc_attachment = get_object_or_404(ClassificationAttachment, pk=pk)

    if thumbnail and vc_attachment.thumbnail_path:
        filename = vc_attachment.thumbnail_path
    else:
        filename = vc_attachment.file.path

    content_type, _ = mimetypes.guess_type(filename)
    with open(filename, "rb") as f:
        image_data = f.read()
    return HttpResponse(image_data, content_type=content_type)


def view_classification_file_attachment_thumbnail(request, pk):
    return view_classification_file_attachment(request, pk, thumbnail=True)


def _get_lab_and_error(user: User):
    lab_error = None
    lab = None
    user_settings = UserSettings.get_for_user(user)
    try:
        lab = user_settings.get_lab()
    except ValueError as ve:
        lab_error = str(ve)
    return lab, lab_error


def create_classification_for_variant(request, variant_id, transcript_id=None):
    if not Classification.can_create_via_web_form(request.user):
        raise PermissionDenied('User cannot create classifications via web form')

    variant = Variant.objects.get(pk=variant_id)
    if variant.is_reference:
        msg = "We only classify variants - not reference alleles"
        raise ValueError(msg)

    variant_sample_autocomplete_form = SampleChoiceForm()
    variant_sample_autocomplete_form.fields['sample'].required = False

    vts = VariantTranscriptSelections(variant,
                                      variant.genome_build,
                                      initial_transcript_id=transcript_id,
                                      add_other_annotation_consortium_transcripts=True)
    lab, lab_error = _get_lab_and_error(request.user)

    consensus = ClassificationConsensus(variant, request.user)

    context = {'variant': variant,
               'variant_sample_autocomplete_form': variant_sample_autocomplete_form,
               "vts": vts,
               "lab": lab,
               "lab_error": lab_error,
               "initially_require_sample": settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE,
               "consensus": consensus}
    return render(request, 'classification/create_classification_for_variant.html', context)


def create_classification_from_hgvs(request, genome_build_name, hgvs_string):
    Classification.check_can_create_no_classification_via_web_form(request.user)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    sample_autocomplete_form = SampleChoiceForm(genome_build=genome_build)
    sample_autocomplete_form.fields['sample'].required = False
    lab, lab_error = _get_lab_and_error(request.user)
    refseq_transcript_accession = ""
    ensembl_transcript_accession = ""
    kind, transcript_accession = get_kind_and_transcript_accession_from_invalid_hgvs(hgvs_string)
    evidence = {}
    if transcript_accession:
        if transcript_accession.startswith("ENST"):
            ensembl_transcript_accession = transcript_accession
        else:
            refseq_transcript_accession = transcript_accession

        EKEYS_BY_KIND = {"c": SpecialEKeys.C_HGVS,
                         "p": SpecialEKeys.P_HGVS,
                         "g": SpecialEKeys.G_HGVS}
        ekey = EKEYS_BY_KIND.get(kind)
        if ekey:
            evidence[ekey] = hgvs_string

    warnings = [f"We could not understand HGVS '{hgvs_string}' - it may be invalid.",
                "This classification WILL NOT BE LINKED TO A VARIANT - no discordance checks can be performed."]
    for w in warnings:
        messages.add_message(request, messages.WARNING, w)

    context = {"genome_build": genome_build,
               "hgvs_string": hgvs_string,
               "evidence_json": json.dumps(evidence),
               "refseq_transcript_accession": refseq_transcript_accession,
               "ensembl_transcript_accession": ensembl_transcript_accession,
               'sample_autocomplete_form': sample_autocomplete_form,
               "lab": lab,
               "lab_error": lab_error,
               "initially_require_sample": settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE}
    return render(request, 'classification/create_classification_from_hgvs.html', context)


@login_not_required
def evidence_keys(request, external_page=True, max_share_level=None):
    """ public page to display EKey details """

    context = {'keys': EvidenceKeyMap.cached().all_keys}

    if external_page:
        context.update({
            "external_page": True,
            "evidence_keys_url": "evidence_keys",
            "evidence_keys_max_share_level_url": "evidence_keys_max_share_level",
            "base_template": 'external_base_content_as_submenu_page_content.html',
        })
    else:
        context.update({
            "evidence_keys_url": "evidence_keys_logged_in",
            "evidence_keys_max_share_level_url": "evidence_keys_logged_in_max_share_level",
            "base_template": 'snpdb/menu/menu_classifications_base.html',
        })
    return render(request, 'classification/evidence_keys.html', context)


def evidence_keys_logged_in(request, max_share_level=None):
    return evidence_keys(request, external_page=False, max_share_level=max_share_level)


def classification_graphs(request):
    if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
        show_unclassified = False
        visibility = "Shared"
        evidence_field = "published_evidence"
    else:
        show_unclassified = True
        visibility = f"Visible to user"
        evidence_field = "evidence"
    classification_counts = get_classification_counts(request.user, show_unclassified=show_unclassified)
    vc_gene_data = get_grouped_classification_counts(request.user, evidence_field, evidence_key="gene_symbol", max_groups=15, show_unclassified=show_unclassified)

    acmg_by_significance = get_criteria_counts(request.user, evidence_field)
    context = {
        "visibility": visibility,
        "show_unclassified": show_unclassified,
        "classification_labels": list(classification_counts.keys())[::-1],
        "classification_values": list(classification_counts.values())[::-1],
        "vc_gene_data": vc_gene_data,
        "acmg_by_significance": acmg_by_significance,
    }
    return render(request, 'classification/classification_graphs.html', context)
