import json
import mimetypes
import re
from datetime import datetime
from typing import Optional, List, Set, Dict

import rest_framework
from crispy_forms.bootstrap import FieldWithButtons
from crispy_forms.layout import Layout, Field, Submit
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.core.exceptions import PermissionDenied
from django.http import StreamingHttpResponse
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.urls import reverse
from django.utils.timezone import now
from django.views.decorators.http import require_POST
from django.views.generic import TemplateView
from global_login_required import login_not_required
from jfu.http import upload_receive, UploadResponse, JFUResponse
from requests.models import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import \
    create_classification_for_sample_and_variant_objects, generate_auto_populate_data
from classification.classification_changes import ClassificationChanges
from classification.classification_stats import get_grouped_classification_counts, \
    get_classification_counts, get_criteria_counts
from classification.enums import SubmissionSource, SpecialEKeys
from classification.models import ClassificationAttachment, Classification, \
    ClassificationRef, ClassificationJsonParams, ClassificationConsensus, ClassificationReportTemplate, ReportNames, \
    ConditionResolvedDict
from classification.models.classification import ClassificationModification
from classification.models.clinical_context_models import ClinicalContext
from classification.models.evidence_key import EvidenceKeyMap
from classification.models.flag_types import classification_flag_types
from classification.views.classification_datatables import ClassificationColumns
from classification.views.classification_export_csv import ExportFormatterCSV
from classification.views.classification_export_redcap import ExportFormatterRedcap
from flags.models import Flag, FlagComment
from flags.models.models import FlagType
from genes.forms import GeneSymbolForm
from genes.hgvs import get_kind_and_transcript_accession_from_invalid_hgvs
from library.django_utils import require_superuser, get_url_from_view_path
from library.file_utils import rm_if_exists
from library.guardian_utils import is_superuser
from library.log_utils import log_traceback
from library.utils import delimited_row
from snpdb.forms import SampleChoiceForm, UserSelectForm, LabSelectForm
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Variant, UserSettings, Sample, Lab, Allele
from snpdb.models.models_genome import GenomeBuild
from uicore.utils.form_helpers import form_helper_horizontal
from variantopedia.forms import SearchAndClassifyForm


@user_passes_test(is_superuser)
def activity(request, latest_timestamp: Optional[str] = None):
    if latest_timestamp:
        latest_timestamp = datetime.fromtimestamp(float(latest_timestamp))
    changes = ClassificationChanges.list_changes(latest_date=latest_timestamp)
    last_date = changes[len(changes)-1].date.timestamp() if changes else None
    context = {
        'changes': changes,
        'last_date': last_date,
        'can_create_classifications': settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN,
        'now': latest_timestamp if latest_timestamp else now()
    }
    return render(request, 'classification/activity.html', context)


def classifications(request):
    """
    Classification listing page
    """

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
        "VARIANT_CLASSIFICATION_REDCAP_EXPORT": settings.VARIANT_CLASSIFICATION_REDCAP_EXPORT,
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
        # TODO should probably know the lab with EvidenceKeys
        for ekey in EvidenceKeyMap.instance().all_keys:
            key_to_order[ekey.key] = index
            index = index + 1

        def key_order(tests_key: str):
            nonlocal key_order
            return key_to_order.get(tests_key, 0) * 2

        complete_values = sorted(complete_values, key=lambda entry: key_order(entry['key']))
        return rest_framework.response.Response(status=HTTP_200_OK, data={'data': complete_values})


@require_POST
def create_classification(request):
    return redirect(create_classification_object(request).get_absolute_url() + "?edit=true")


def create_classification_object(request) -> Classification:
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

    return classification


def classification_history(request, record_id):
    ref = ClassificationRef.init_from_str(request.user, str(record_id))
    ref.check_exists()
    #  changes page will show things before they were submitted and not hide values
    #  so need to restrict to people who've always had access to the record
    ref.check_security(must_be_writable=True)

    changes = ClassificationChanges.list_changes(classification=ref.record, limit=100)
    context = {
        'changes': changes,
        'can_create_classifications': Classification.can_create_via_web_form(request.user),
        'now': now()
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
        if len(variant.genome_builds) == 1:
            genome_build = next(iter(variant.genome_builds))

    if genome_build is None:
        user_settings = UserSettings.get_for_user(request.user)
        genome_build = user_settings.default_genome_build

    vc: Classification = ref.record

    context = {
        'vc': vc,
        'record': record,
        'genome_build': genome_build,
        'existing_files': existing_files,
        'other_classifications_summary': other_classifications_summary,
        'report_enabled': ClassificationReportTemplate.objects.filter(name=ReportNames.DEFAULT_REPORT).exclude(template__iexact='').exists(),
        'attachments_enabled': settings.VARIANT_CLASSIFICATION_FILE_ATTACHMENTS
    }
    return render(request, 'classification/classification.html', context)


def view_classification_diff(request):
    extra_data = dict()

    if history_str := request.GET.get('history'):
        vc_id = int(history_str)
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

    elif clinical_context_str := request.GET.get('clinical_context'):
        cc = ClinicalContext.objects.get(pk=clinical_context_str)
        records = cc.classification_modifications
        records.sort(key=lambda cm: cm.curated_date_check, reverse=True)

    elif variant_id_str := request.GET.get('variant_compare'):
        ref = ClassificationRef.init_from_str(user=request.user, id_str=variant_id_str)
        compare_all = ClassificationModification.latest_for_user(user=request.user, allele=ref.record.variant.allele, published=True)
        # filter out variant we're comparing with, make it last in calculation
        latest_others_for_variant = [vcm for vcm in compare_all if vcm.classification.id != ref.record.id]
        latest_others_for_variant.sort(key=lambda cm: cm.curated_date_check, reverse=True)
        compare_all = [ref.modification] + latest_others_for_variant
        records = compare_all

    elif allele_id_str := request.GET.get('allele'):
        allele_id = int(allele_id_str)
        compare_all = list(ClassificationModification.latest_for_user(user=request.user, allele=Allele.objects.get(pk=allele_id), published=True))
        compare_all.sort(key=lambda cm: cm.curated_date_check, reverse=True)
        records = compare_all

    elif cids := request.GET.get('cids'):
        records = [ClassificationModification.latest_for_user(user=request.user, classification=cid, published=True).first() for cid in [cid.strip() for cid in cids.split(',')]]
        records = [record for record in records if record]

    else:
        raise ValueError("Diff not given valid diff type")

    # filter out any Nones inserted by filtering on user permission etc
    records = [record for record in records if record]

    records_json = [vcm.as_json(ClassificationJsonParams(
        current_user=request.user,
        include_data=True,
        hardcode_extra_data=extra_data.get(vcm.id),
        fix_data_types=True,
        api_version=1,
    )) for vcm in records]

    for record in records_json:
        condition_dict: ConditionResolvedDict
        if condition_dict := record.get('resolved_condition'):
            condition_ekey = record.get('data').get(SpecialEKeys.CONDITION, {})
            condition_ekey["resolved"] = condition_dict
            record.get('data')[SpecialEKeys.CONDITION] = condition_ekey

    context = {
        'records_json': records_json
    }
    return render(request, 'classification/classification_diff.html', context)


def classification_import_tool(request: HttpRequest) -> Response:
    """
    This page allows you to upload a file, the fact that the file is probably to do with classifications
    doesn't really matter, there's not logic to process it here, just to upload it
    """
    all_labs = list(Lab.valid_labs_qs(request.user, admin_check=True))
    selected_lab = None
    if len(all_labs) == 1:
        selected_lab = all_labs[0]
    else:
        # choose the first test lab for import to process on
        for lab in all_labs:
            if 'test' in lab.name.lower():
                selected_lab = lab
                break

    context = {
        "labs": all_labs,
        "selected_lab": selected_lab
    }
    return render(request, 'classification/classification_import_tool.html', context)


def classification_qs(request):
    config = ClassificationColumns(request)
    qs = ClassificationModification.latest_for_user(user=request.user, published=True)
    qs = config.filter_queryset(qs)
    # TODO sort
    return qs


def export_classifications_grid(request):
    """
    CSV export of what is currently filtered into the classification grid
    """
    genome_build = UserSettings.get_for_user(request.user).default_genome_build
    qs = classification_qs(request)
    return ExportFormatterCSV(user=request.user, genome_build=genome_build, qs=qs).export()


def export_classifications_grid_redcap(request):
    genome_build = UserSettings.get_for_user(request.user).default_genome_build
    qs = classification_qs(request)
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


def get_classification_attachment_file_dicts(classification) -> List[Dict]:
    file_dicts: List[Dict] = list()
    vca: ClassificationAttachment
    for vca in classification.classificationattachment_set.all():
        file_dicts.append(vca.get_file_dict())

    file_dicts = list(reversed(file_dicts))  # JFU adds most recent at the end
    return file_dicts


def view_classification_file_attachment(request, pk, thumbnail=False) -> HttpResponse:
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


class CreateClassificationForVariantView(TemplateView):
    template_name = 'classification/create_classification_for_variant.html'

    def _get_variant(self) -> Variant:
        return Variant.objects.get(pk=self.kwargs["variant_id"])

    def _get_genome_build(self) -> GenomeBuild:
        return GenomeBuild.get_name_or_alias(self.kwargs["genome_build_name"])

    def _get_form_post_url(self):
        return reverse("create_classification")

    def _get_sample_form(self):
        variant_sample_autocomplete_form = SampleChoiceForm()
        variant_sample_autocomplete_form.fields['sample'].required = False
        return variant_sample_autocomplete_form

    def get_context_data(self, *args, **kwargs):
        if not Classification.can_create_via_web_form(self.request.user):
            raise PermissionDenied('User cannot create classifications via web form')

        variant = self._get_variant()
        if variant.is_reference:
            msg = "We only classify variants - not reference alleles"
            raise ValueError(msg)

        genome_build = self._get_genome_build()
        vts = VariantTranscriptSelections(variant, genome_build,
                                          add_other_annotation_consortium_transcripts=True)
        lab, lab_error = UserSettings.get_lab_and_error(self.request.user)

        consensus = ClassificationConsensus(variant, self.request.user)
        return {'variant': variant,
                "genome_build": genome_build,
                "form_post_url": self._get_form_post_url(),
                'variant_sample_autocomplete_form': self._get_sample_form(),
                "vts": vts,
                "lab": lab,
                "lab_error": lab_error,
                "initially_require_sample": settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_INITIALLY_REQUIRE_SAMPLE,
                "consensus": consensus}


def create_classification_from_hgvs(request, genome_build_name, hgvs_string):
    Classification.check_can_create_no_classification_via_web_form(request.user)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    sample_autocomplete_form = SampleChoiceForm(genome_build=genome_build)
    sample_autocomplete_form.fields['sample'].required = False
    lab, lab_error = UserSettings.get_lab_and_error(request.user)
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
def evidence_keys(request: HttpRequest) -> HttpResponse:
    """ public page to display EKey details """

    context = {
        'keys': EvidenceKeyMap.instance().all_keys
    }
    return render(request, 'classification/evidence_keys.html', context)


def classification_graphs(request):
    if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
        show_unclassified = False
        visibility = "shared & matched to an allele"
    else:
        show_unclassified = True
        visibility = f"visible to {request.user.username}"
    classification_counts = get_classification_counts(request.user, show_unclassified=show_unclassified, unique_alleles=True)

    evidence_field = "published_evidence"
    vc_gene_data = get_grouped_classification_counts(
        user=request.user,
        field=evidence_field,
        evidence_key="gene_symbol",
        max_groups=15,
        show_unclassified=show_unclassified,
        allele_level=True)

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


def lab_gene_classification_counts(request):
    if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
        visibility = "Shared"
    else:
        visibility = f"Visible to user"

    data_columns_whitelist = {

    }
    gene_list_categories_whitelist = {
        "Panel App Panel",
        "Lab Gene Classification Counts",
    }

    context = {"visibility": visibility,
               "data_columns_whitelist": data_columns_whitelist,
               "gene_list_categories_whitelist": gene_list_categories_whitelist,
               "default_enrichment_kits": []}
    return render(request, 'classification/lab_gene_classification_counts.html', context)


@require_superuser
def clin_sig_change_data(request):

    def yield_data():

        genome_build = GenomeBuildManager.get_current_genome_build()
        yield delimited_row(['Classification First Submitted (UTC)', 'Date Changed (UTC)', 'Org', 'Lab', f'c.hgvs {genome_build}', 'URL', 'From', 'To', 'Status', 'Comments', 'Discordance(s)', 'Other Labs for Allele at Time'], '\t')

        flag_changed_re = re.compile(r"^Classification (has )?changed from (?P<from>.*?) to (?P<to>.*?)$")
        de_number_re = re.compile(r"(.*?) [(].*?[)]")

        def de_number(clin_sig: str) -> str:
            if clin_sig:
                clin_sig = clin_sig.strip()
                if match := de_number_re.match(clin_sig):
                    return match.group(1)
                return clin_sig

        all_flags_qs = Flag.objects.filter(flag_type=classification_flag_types.significance_change)
        flag: Flag
        for flag in all_flags_qs.order_by('-created'):
            flag_comment: FlagComment
            cs_from: str = 'ERROR'
            cs_to: str = 'ERROR'
            comments: List[str] = list()
            resolution: Optional[str] = 'In Progress'
            classification_created: datetime
            date_raised: datetime = flag.created
            discordance_dates: List[datetime] = list()
            other_labs: Set[Lab] = set()

            flag_collection = flag.collection
            if source := flag_collection.source_object:
                if isinstance(source, Classification):
                    c_hgvs = source.get_c_hgvs(genome_build=GenomeBuildManager.get_current_genome_build())
                    org = source.lab.organization.name
                    lab = source.lab.name
                    url = get_url_from_view_path(source.get_absolute_url())
                    classification_created = source.created

                    # get earliest discordance
                    first_discordance: Flag
                    for discordance in Flag.objects.filter(collection=source.flag_collection_safe,
                                                           flag_type=classification_flag_types.discordant).order_by(
                            'created'):
                        discordance_dates.append(discordance.created)

                    if allele := source.allele:
                        cl: Classification
                        for cl in Classification.objects.filter(variant__in=allele.variants):
                            if cl.lab != source.lab and cl.created < flag.created:
                                other_labs.add(cl.lab)

                else:
                    continue
            else:
                # the classification has been deleted
                continue

            for index, flag_comment in enumerate(flag.flagcomment_set.order_by('created')):
                if text := flag_comment.text:
                    if index == 0:
                        if match := flag_changed_re.match(text):
                            cs_from = de_number(match.group('from'))
                            cs_to = de_number(match.group('to'))
                        else:
                            print(text)
                    else:
                        comments.append(text)
                if flag_comment and flag_comment.resolution:
                    resolution = flag_comment.resolution.label

            other_lab_list = list(other_labs)
            other_lab_list.sort()
            other_lab_str = ", ".join(other_lab.name for other_lab in other_lab_list)
            discordance_date_str = " ".join([dd.strftime('%Y-%m-%d %H:%M:%S') for dd in discordance_dates])
            yield delimited_row([classification_created.strftime('%Y-%m-%d %H:%M:%S'), date_raised.strftime('%Y-%m-%d %H:%M:%S'), org, lab, c_hgvs, url, cs_from, cs_to, resolution, '\n'.join(comments), discordance_date_str, other_lab_str], '\t')

    response = StreamingHttpResponse(yield_data(), content_type='text/tsv')
    # modified_str = now().strftime("%a, %d %b %Y %H:%M:%S GMT")  # e.g. 'Wed, 21 Oct 2015 07:28:00 GMT'
    # response['Last-Modified'] = modified_str
    response['Content-Disposition'] = f'attachment; filename="clin_sig_changes.tsv"'
    return response
