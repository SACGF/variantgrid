import io
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Dict, Any, Optional, Iterable, List, Set

from django.contrib import messages
from django.db.models import QuerySet, When, Value, Case, IntegerField, Count, Q, TextField
from django.db.models.functions import Cast
from django.http import HttpResponse, StreamingHttpResponse, HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render, redirect, get_object_or_404
from django.urls import reverse
from django.views import View
from django.views.decorators.http import require_POST

from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, ClinVarExportBatch, ClinVarExportBatchStatus, \
    EvidenceKeyMap, ClinVarExportStatus, ClinVarExportSubmission, ConditionResolved
from classification.models.clinvar_export_prepare import ClinvarExportPrepare
from classification.utils.clinvar_matcher import ClinVarLegacyRow, ClinVarLegacyMatches, ClinVarLegacyExportMatchType
from classification.views.classification_dashboard_view import ClassificationDashboard
from genes.hgvs import CHGVS
from library.cache import timed_cache
from library.django_utils import add_save_message, get_url_from_view_path, require_superuser, RequireSuperUserView
from library.utils import html_to_text, export_column, ExportRow, local_date_string, ExportDataType, JsonDataType
from library.utils.django_utils import render_ajax_view
from ontology.models import OntologyTerm, AncestorCalculator, OntologyTermRelation
from snpdb.lab_picker import LabPickerData
from snpdb.models import ClinVarKey, Lab, Allele, GenomeBuild
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, CellData


@timed_cache(size_limit=30, ttl=60)
def allele_for(allele_id: int) -> Allele:
    return Allele.objects.select_related('clingen_allele').get(pk=allele_id)


class ClinVarExportBatchColumns(DatatableConfig):

    def render_status(self, row: CellData):
        return ClinVarExportBatchStatus(row['status']).label

    def batch_id(self, row: CellData):
        return row['id']

    def __init__(self, request):
        super().__init__(request)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('clinvar_export_batch_detail', expected_height=120)
        self.rich_columns = [
            RichColumn("id", label="ID", orderable=True, renderer=self.batch_id, client_renderer='batchId', default_sort=SortOrder.DESC),
            RichColumn("clinvar_key", label="ClinVar Key", orderable=True, enabled=False),
            RichColumn("allele_origin_bucket", label="Allele Origin", orderable=True, client_renderer="render_allele_origin_bucket"),
            RichColumn("created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn("record_count", label="Records In Submission", orderable=True),
            RichColumn("status", renderer=self.render_status, orderable=True)
        ]

    def get_initial_query_params(self, clinvar_key: Optional[str] = None):
        clinvar_keys = ClinVarKey.clinvar_keys_for_user(self.user)
        if clinvar_key:
            clinvar_keys = clinvar_keys.filter(pk=clinvar_key)
        cve = ClinVarExportBatch.objects.filter(clinvar_key__in=clinvar_keys)
        return cve

    def get_initial_queryset(self):
        initial_qs = self.get_initial_query_params(
            clinvar_key=self.get_query_param('clinvar_key')
        )
        initial_qs = initial_qs.annotate(record_count=Count('clinvarexportsubmission'))
        return initial_qs


def clinvar_export_batch_detail(request, clinvar_export_batch_id: int):
    clinvar_export_batch: ClinVarExportBatch = ClinVarExportBatch.objects.get(pk=clinvar_export_batch_id)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)
    submissions = clinvar_export_batch.clinvarexportsubmission_set.order_by('created')

    return render_ajax_view(request, 'classification/clinvar_export_batch_detail.html', {
        "batch": clinvar_export_batch,
        "submissions": submissions
    })


def _export_id_to_batch_ids(qs: QuerySet[ClinVarExport]) -> Dict[int, List[int]]:
    export_to_batches = defaultdict(list)
    for export_id, batch_id in ClinVarExportSubmission.objects.filter(
        clinvar_export__in=qs.values_list("id", flat=True)).order_by(
        "-submission_batch").values_list("clinvar_export", "submission_batch"):
        export_to_batches[export_id].append(batch_id)
    return export_to_batches


class ClinVarExportColumns(DatatableConfig[ClinVarExport]):

    def render_allele(self, row: Dict[str, Any]) -> str:
        allele = allele_for(row["clinvar_allele__allele"])
        return f"{allele:CA}"

    def pre_render(self, qs: QuerySet[ClinVarExport]):
        # find all the batches these records are in
        # do this once rather than per row
        self.export_to_batches = _export_id_to_batch_ids(qs)

    def power_search(self, qs: QuerySet[ClinVarExport], search_string: str) -> QuerySet[ClinVarExport]:
        try:
            # if searching on a number could be batch ID or submission ID
            if search_string.isnumeric():
                batch_id = int(search_string)
                submissions = ClinVarExportSubmission.objects.filter(submission_batch_id=batch_id).values_list("clinvar_export_id", flat=True)
                # filter rather than just return submissions to make sure user has permission to see items that belong to the batch
                # also make sure an individual export doesn't have the batch ID
                return qs.filter(Q(pk__in=submissions) | Q(pk=batch_id))
            elif search_string.upper().startswith("CE_"):
                batch_id = search_string[3:]
                return qs.filter(Q(pk=batch_id))

            # if searcing a ClinGenAlleleID
            elif search_string.upper().startswith("CA"):
                clingen_number_str = search_string[2:]
                # make sure we're actually searching on a number and not "cat disease"
                if clingen_number_str:  # if user is just typing 0s, don't filter anything yet
                    if clingen_number := int(clingen_number_str):
                        return qs.annotate(
                            clingen_allele_id_str=(Cast('clinvar_allele__allele__clingen_allele_id', output_field=TextField()))
                        ).filter(clingen_allele_id_str__startswith=str(clingen_number))
                    else:
                        return qs
        except ValueError:
            pass

        # if searching on a c.hgvs, MONDO ID
        return super().power_search(qs, search_string)

    def batches(self, row: CellData) -> JsonDataType:
        return self.export_to_batches.get(row.get("id"))

    def render_c_hgvs(self, row: CellData) -> JsonDataType:
        if row["classification_based_on__classification__allele_info__allele"]:
            genome_build = row["classification_based_on__published_evidence__genome_build__value"]
            c_hgvs_str: str
            if "h37" in genome_build:
                c_hgvs_str = row["classification_based_on__classification__allele_info__grch37__c_hgvs"]
            else:
                c_hgvs_str = row["classification_based_on__classification__allele_info__grch38__c_hgvs"]

            data: Dict[str, Any]
            c_hgvs = CHGVS(c_hgvs_str)
            if c_hgvs.raw_c != c_hgvs.full_c_hgvs:
                data = {
                    "allele_origin_bucket": row["allele_origin_bucket"],
                    "genome_build": genome_build,
                    "transcript": c_hgvs.transcript,
                    "gene_symbol": c_hgvs.gene_symbol,
                    "variant": c_hgvs.raw_c,
                    "always_show_genome_build": True
                    # below row makes this link to the variant, but probably not desired action
                    # "variantId": row["classification_based_on__classification__variant"]
                }
            else:
                data = {"full": c_hgvs_str}

            if allele_id := row["clinvar_allele__allele"]:
                allele = allele_for(allele_id)
                data["allele"] = f"{allele:CA}"

            return data
        else:
            return None

    def render_id(self, row: CellData) -> JsonDataType:
        return {
            "id": row["id"],
            "allele_origin_bucket": row["allele_origin_bucket"]
        }

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.export_to_batches = {}

        self.search_box_enabled = True
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('clinvar_export_detail')
        self.rich_columns = [
            RichColumn("id",
                       label="ID",
                       orderable=True,
                       default_sort=SortOrder.DESC,
                       client_renderer='renderId',
                       renderer=self.render_id,
                       extra_columns=["allele_origin_bucket"]
            ),
            RichColumn(name="c_hgvs", label='Allele',
                    sort_keys=["classification_based_on__classification__allele_info__grch38__c_hgvs"],
                    extra_columns=[
                        "clinvar_allele__allele",
                        "classification_based_on__classification__allele_info__allele",
                        "classification_based_on__published_evidence__genome_build__value",
                        "classification_based_on__classification__allele_info__grch37__c_hgvs",
                        "classification_based_on__classification__allele_info__grch38__c_hgvs",
                    ],
                    renderer=self.render_c_hgvs, client_renderer='VCTable.hgvs',
                    search=[
                        # "clinvar_allele__allele__clingen_allele__id",  # need a string
                        "classification_based_on__classification__allele_info__grch37__c_hgvs",
                        "classification_based_on__classification__allele_info__grch38__c_hgvs"
                    ]
            ),
            RichColumn("condition",
                       label="Condition Umbrella",
                       client_renderer='VCTable.condition',
                       search=["condition__display_text"],
                       sort_keys=["condition__sort_text"]
            ),
            RichColumn(name="batches", label="In Batch IDs", renderer=self.batches, client_renderer='renderBatches', orderable=False, search=False, extra_columns=["id"]),
            RichColumn("status", label="Sync Status", client_renderer='renderStatus', sort_keys=["status_sort"], orderable=True, search=False),
            RichColumn("scv", label="SCV", orderable=True),
        ]

    def get_initial_query_params(self, clinvar_key: Optional[str] = None, status: Optional[str] = None) -> QuerySet[ClinVarExport]:
        clinvar_keys = ClinVarKey.clinvar_keys_for_user(self.user)
        if clinvar_key:
            clinvar_keys = clinvar_keys.filter(pk=clinvar_key)

        cve = ClinVarExport.objects.filter(clinvar_allele__clinvar_key__in=clinvar_keys)
        if status:
            statuses = status.split(',')
            cve = cve.filter(status__in=statuses)

        return cve

    def get_initial_queryset(self) -> QuerySet[ClinVarExport]:
        initial_qs = self.get_initial_query_params(
            clinvar_key=self.get_query_param('clinvar_key'),
            status=self.get_query_param('status')
        )
        # add sorting for status so important changes show first
        whens = [
            When(status=ClinVarExportStatus.NEW_SUBMISSION, then=Value(1)),
            When(status=ClinVarExportStatus.CHANGES_PENDING, then=Value(2)),
            When(status=ClinVarExportStatus.UP_TO_DATE, then=Value(3)),
            When(status=ClinVarExportStatus.IN_ERROR, then=Value(4)),
            When(status=ClinVarExportStatus.EXCLUDE, then=Value(5))
        ]
        case = Case(*whens, default=Value(0), output_field=IntegerField())
        initial_qs = initial_qs.annotate(status_sort=case)

        return initial_qs


def clinvar_export_review(request: HttpRequest, clinvar_export_id: int) -> HttpResponseBase:
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=clinvar_export_id)  # fixme get or 404
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    if request.method == "POST":
        clinvar_export.scv = request.POST.get("scv") or ""
        clinvar_export.save()
        add_save_message(request, valid=True, name="ClinVarExport")
        return redirect(clinvar_export)

    common_condition: Optional[OntologyTerm] = None
    clinvar_exports_for_allele = list(ClinVarExport.objects.filter(
        clinvar_allele_id=clinvar_export.clinvar_allele_id
    ).exclude(pk=clinvar_export_id).all())
    if clinvar_exports_for_allele:
        all_condition_terms: Set[OntologyTerm] = set()
        clinvar_exports_for_allele.append(clinvar_export)
        for clinvar_allele in clinvar_exports_for_allele:
            for term in clinvar_allele.condition_resolved.terms:
                if mondo := OntologyTermRelation.as_mondo(term):
                    all_condition_terms.add(mondo)

        clinvar_exports_for_allele = sorted(clinvar_exports_for_allele, key=lambda cve: cve.pk)
        common_condition = AncestorCalculator.common_ancestor(all_condition_terms)

    return render(request, 'classification/clinvar_export.html', context={
        "clinvar_export": clinvar_export,
        "cm": clinvar_export.classification_based_on,
        "clinvar_exports_for_allele": clinvar_exports_for_allele,
        "common_condition": common_condition
    })


def clinvar_export_history(request: HttpRequest, clinvar_export_id: int) -> HttpResponseBase:
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=clinvar_export_id)
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    history = clinvar_export.clinvarexportsubmission_set.order_by('-created')

    return render(request, 'classification/clinvar_export_history.html', context={
        'clinvar_export': clinvar_export,
        'history': history
    })


@dataclass(frozen=True)
class ClinVarExportSummaryData:
    clinvar_export: ClinVarExport
    batch_ids: Optional[List[int]]


@dataclass(frozen=True)
class ClinVarExportSummary(ExportRow):

    clinvar_export: ClinVarExport
    batch_ids: Optional[List[int]]
    has_duplicates: bool = False

    @property
    def classification(self):
        return self.clinvar_export.classification_based_on

    @cached_property
    def genome_build_row(self) -> GenomeBuild:
        if classification := self.classification:
            try:
                return classification.get_genome_build()
            except KeyError:
                pass

    @export_column("ID")
    def _id(self):
        return self.clinvar_export.clinvar_export_id

    @export_column("ClinVar Export URL")
    def _clinvar_export_url(self):
        return get_url_from_view_path(self.clinvar_export.get_absolute_url())

    @export_column("Allele URL")
    def _allele_url(self):
        if modification := self.classification:
            classification = modification.classification
            if allele_id := classification.allele_id:
                return get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": allele_id}))

    @export_column("Classification URL")
    def _classification_url(self):
        if modification := self.classification:
            return get_url_from_view_path(reverse('view_classification', kwargs={"classification_id": modification.id_str}))

    @export_column("Record ID")
    def _record_id(self):
        if modification := self.classification:
            classification = modification.classification
            return classification.lab.group_name + "/" + classification.lab_record_id

    @export_column("Allele Origin")
    def _allele_origin(self):
        return self.clinvar_export.get_allele_origin_bucket_display()

    @export_column("Genome Build")
    def _genome_build(self):
        return str(self.genome_build_row) if self.genome_build_row else None

    @export_column("ClinGenAllele")
    def _clingen_allele(self):
        if modification := self.classification:
            if allele := modification.classification.allele_object:
                return str(allele.clingen_allele)

    @export_column("c.HGVS")
    def _c_hgvs(self):
        if genome_build := self.genome_build_row:
            return self.classification.classification.get_c_hgvs(genome_build)

    @export_column("Condition Umbrella")
    def _condition_umbrella(self):
        return self.clinvar_export.condition_resolved.as_plain_text

    @export_column("Condition")
    def _condition(self):
        if modification := self.classification:
            return modification.condition_text

    @export_column("Affected Status")
    def _affected_status(self):
        if modification := self.classification:
            return EvidenceKeyMap.pretty_value_for(modification, SpecialEKeys.AFFECTED_STATUS)

    @export_column("Mode of Inheritance")
    def _mode_of_inheritance(self):
        if modification := self.classification:
            return EvidenceKeyMap.pretty_value_for(modification, SpecialEKeys.MODE_OF_INHERITANCE)

    @export_column("Clinical Significance")
    def _clinical_significance(self):
        if modification := self.classification:
            return EvidenceKeyMap.pretty_value_for(modification, SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @export_column("Interpretation Summary")
    def _interpretation_summary(self):
        if classification := self.classification:
            return html_to_text(classification.get(SpecialEKeys.INTERPRETATION_SUMMARY))

    @export_column("Curation Date")
    def _curation_date(self):
        if classification := self.classification:
            return EvidenceKeyMap.pretty_value_for(classification, SpecialEKeys.CURATION_DATE)

    @export_column("Classification Submitted to $site_name", data_type=ExportDataType.datetime)
    def _classification_imported_created(self):
        if modification := self.classification:
            return modification.created

    @export_column("Sync Status")
    def _sync_status(self):
        if self.clinvar_export.status == ClinVarExportStatus.UP_TO_DATE:
            if clinvar_error := self.clinvar_export.last_submission_error:
                return "ClinVar Error"
        return self.clinvar_export.get_status_display()

    @export_column("SCV")
    def _scv(self):
        return self.clinvar_export.scv

    @export_column("Latest Batch ID")
    def _batch_id(self):
        if submission := self.clinvar_export.last_submission:
            return submission.submission_batch_id

    @export_column("Latest Batch Status")
    def _batch_status(self):
        if submission := self.clinvar_export.last_submission:
            submission.submission_batch.get_status_display()

    @export_column("All Batch IDs")
    def _batch_ids_all(self):
        if batch_ids := self.batch_ids:
            return ",".join(str(batch_id) for batch_id in sorted(batch_ids))

    @export_column("Others for Allele")
    def _possible_duplicate(self):
        if self.has_duplicates:
            duplicates = ClinVarExport.objects.filter(clinvar_allele=self.clinvar_export.clinvar_allele).exclude(
                pk=self.clinvar_export.pk).order_by('pk')
            return ", ".join(str(duplicate) for duplicate in duplicates)

    @export_column("Common Condition w Others")
    def _common_condition(self):
        if self.has_duplicates:
            all_for_allele = ClinVarExport.objects.filter(clinvar_allele=self.clinvar_export.clinvar_allele).order_by('pk')
            all_conditions = list()
            for other in all_for_allele:
                all_conditions += other.condition_resolved.terms
            all_conditions = list(sorted(set(all_conditions)))
            all_conditions_str = ", ".join(f"{term.pk} {term.name}" for term in all_conditions)
            try:
                if common_condition := AncestorCalculator.common_ancestor(all_conditions):
                    return str(common_condition) + f" from ({all_conditions_str})"
            except ValueError:
                return f"No common condition found from ({all_conditions_str})"
        return ""

    @export_column("Messages")
    def _messages(self):
        all_messages: List[str] = []
        if clinvar_error := self.clinvar_export.last_submission_error:
            all_messages.append(f"(CLINVAR ERROR) {clinvar_error}")

        if json_header := self.clinvar_export.submission_grouping:
            if j_messages := json_header.all_messages:
                all_messages += [f"({message.severity.upper()}) {message.text}" for message in j_messages if
                                 message.severity != "info"]

        if json_body := self.clinvar_export.submission_body:
            if j_messages := json_body.all_messages:
                all_messages += [f"({message.severity.upper()}) {message.text}" for message in j_messages if message.severity != "info"]

        return "\n".join(all_messages)

    ## This column is a bit much
    # @export_column("JSON")
    # def json(self):
    #     if json_version := self.clinvar_export.submission_full:
    #         return json.dumps(json_version.pure_json(), indent=4)

def clinvar_export_download(request: HttpRequest, clinvar_key_id: str) -> HttpResponseBase:
    # get all the individual rows
    clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)
    qs = ClinVarExport.objects.filter(clinvar_allele__clinvar_key=clinvar_key).order_by('-clinvar_allele__pk')

    # get all map of export_id to all the batch ids that export id was in
    export_id_to_batch_ids = _export_id_to_batch_ids(qs)

    # get a list of duplicates
    duplicates = ClinVarExport.objects.filter(
        # classification_based_on__isnull=False,
        clinvar_allele__clinvar_key_id=clinvar_key_id
    ).values('clinvar_allele').annotate(total=Count('clinvar_allele')).order_by('-total')

    has_possible_duplicate_clinvar_allele_id: set[int] = set()

    for duplicate in duplicates:
        total = duplicate.get("total")
        if total <= 1:
            break
        else:
            has_possible_duplicate_clinvar_allele_id.add(duplicate.get('clinvar_allele'))

    data_iterator = (
        ClinVarExportSummary(
            clinvar_export=ce,
            batch_ids=export_id_to_batch_ids.get(ce.pk),
            has_duplicates=ce.clinvar_allele_id in has_possible_duplicate_clinvar_allele_id
        ) for ce in qs.iterator()
    )
    csv_generator = ClinVarExportSummary.csv_generator(data_iterator)

    date_str = local_date_string()
    response = StreamingHttpResponse(csv_generator, content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="export_batches_{date_str}.csv"'
    return response


def clinvar_export_batch_download(request: HttpRequest, clinvar_export_batch_id: int) -> HttpResponseBase:
    clinvar_export_batch: ClinVarExportBatch = ClinVarExportBatch.objects.get(pk=clinvar_export_batch_id)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)

    # code duplicated from admin, but don't feel this should go into models
    batch_json = clinvar_export_batch.to_json()
    batch_json_str = json.dumps(batch_json)
    response = HttpResponse(batch_json_str, content_type='application/json')
    response['Content-Disposition'] = f'attachment; filename=clinvar_export_preview_{clinvar_export_batch.pk}.json'
    return response


def clinvar_export_summary(request: HttpRequest, clinvar_key_id: Optional[str] = None) -> HttpResponseBase:
    clinvar_key: ClinVarKey
    if not clinvar_key_id:
        if clinvar_key := ClinVarKey.clinvar_keys_for_user(request.user).first():
            return redirect(reverse('clinvar_key_summary', kwargs={'clinvar_key_id': clinvar_key.pk}))
        else:
            # page has special support if no clinvar key is available to the user
            return render(request, 'classification/clinvar_key_summary_none.html')

    clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    labs = Lab.objects.filter(clinvar_key=clinvar_key).order_by('name')
    dashboard = ClassificationDashboard(LabPickerData.for_labs(labs, user=request.user))

    export_columns = ClinVarExportColumns(request)
    export_batch_columns = ClinVarExportBatchColumns(request)

    return render(request, 'classification/clinvar_key_summary.html', {
        'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
        'clinvar_key': clinvar_key,
        'labs': labs,
        'missing_condition_count': dashboard.classifications_wout_standard_text,
        'count_records': export_columns.get_initial_query_params(clinvar_key=clinvar_key_id).count(),
        'count_batch': export_batch_columns.get_initial_query_params(clinvar_key=clinvar_key_id).count()
    })


@require_POST
def clinvar_export_refresh(request: HttpRequest, clinvar_key_id: str) -> HttpResponseBase:
    clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)
    logs = ClinvarExportPrepare.update_export_records_for_keys(clinvar_keys={clinvar_key, })
    print("\n".join(logs))
    messages.add_message(request, level=messages.INFO, message="Prepare complete")
    # TODO, actually store or display the logs somewhere - against the ClinVarAlleles?
    return redirect(reverse('clinvar_key_summary', kwargs={'clinvar_key_id': clinvar_key_id}))


@require_POST
def clinvar_export_create_batch(request: HttpRequest, clinvar_key_id: str) -> HttpResponseBase:
    clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    clinvar_export_ids_str = request.POST.get('clinvar_export_ids')
    all_ids = set(int(cid) for cid in re.compile(r"\d+").findall(clinvar_export_ids_str))
    clinvar_export_qs = ClinVarExport.objects.filter(clinvar_allele__clinvar_key=clinvar_key, pk__in=all_ids)
    id_count = clinvar_export_qs.count()
    if clinvar_export_qs.count() != len(all_ids):
        # find the IDs that didn't exist for the lab
        found_ids = set(clinvar_export_qs.values_list('pk', flat=True).all())
        missing = sorted(all_ids - found_ids)
        missing_str = ", ".join(str(m) for m in missing)
        messages.add_message(request, level=messages.ERROR, message=f"Could not find some of the IDs for this ClinVarKey - not creating a batch. Missing IDs : {missing_str}")
    elif id_count == 0:
        messages.add_message(request, level=messages.ERROR, message=f"No IDs provided")
    else:
        batches = ClinVarExportBatch.create_batches(clinvar_export_qs, force_update=True)
        for batch in batches:
            batch_records_count = batch.clinvarexportsubmission_set.count()
            messages.add_message(request, level=messages.SUCCESS, message=f"Created Export Batch {batch} with {batch_records_count} records")
            if missing := len(all_ids) - batch_records_count:
                batch_record_ids = set(batch.clinvarexportsubmission_set.values_list('clinvar_export_id', flat=True).all())
                missing_ids = sorted(all_ids - batch_record_ids)
                missing_ids_list = ', '.join(str(i) for i in missing_ids)
                messages.add_message(request, level=messages.WARNING, message=f"Some records were not added to the batch, already up to date or in error : {missing_ids_list}")
        if not batches:
            lab_record_ids = ', '.join(str(i) for i in all_ids)
            messages.add_message(request, level=messages.ERROR, message=f"{lab_record_ids} IDs were already up to date or in error.")
    return redirect(reverse('clinvar_key_summary', kwargs={'clinvar_key_id': clinvar_key_id}))


def clinvar_export_detail(request: HttpRequest, clinvar_export_id: int) -> HttpResponseBase:
    clinvar_export = get_object_or_404(ClinVarExport, pk=clinvar_export_id)
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    interpretation_summary: Optional[str] = None
    if cm := clinvar_export.classification_based_on:
        if interpret_summary_html := cm.get(SpecialEKeys.INTERPRETATION_SUMMARY):
            interpretation_summary = html_to_text(interpret_summary_html)

    return render(request, 'classification/clinvar_export_detail.html', {
        'clinvar_export': clinvar_export,
        'classification': clinvar_export.classification_based_on,
        'interpretation_summary': interpretation_summary
    })


class ClinVarMatchView(RequireSuperUserView, View):

    def get(self, request, **kwargs):
        clinvar_key_id = kwargs.get('clinvar_key_id')
        clinvar_key: ClinVarKey
        if not clinvar_key_id:
            if clinvar_key := ClinVarKey.clinvar_keys_for_user(request.user).first():
                return redirect(reverse('clinvar_match', kwargs={'clinvar_key_id': clinvar_key.pk}))
            else:
                # page has special support if no clinvar key is available to the user
                return render(request, 'classification/clinvar_key_summary_none.html')

        clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
        clinvar_key.check_user_can_access(request.user)

        return render(request, 'classification/clinvar_match.html', {
            'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
            'clinvar_key': clinvar_key
        })

    def post(self, request, **kwargs):
        clinvar_key_id = kwargs.get('clinvar_key_id')
        clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
        clinvar_key.check_user_can_access(request.user)

        clinvar_legacy_rows: List[ClinVarLegacyRow] = []
        try:
            file_obj = io.StringIO(request.FILES.get('file').read().decode("utf-8"))
            clinvar_legacy_rows = list(ClinVarLegacyRow.load_file(file_obj, clinvar_key))
        except AttributeError as ve:
            messages.error(request, f"File is missing or in the wrong format ({ve})")

        return render(request, 'classification/clinvar_match.html', {
            'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
            'clinvar_key': clinvar_key,
            'rows': clinvar_legacy_rows
        })


@require_superuser
def clinvar_match_detail(request, clinvar_key_id: str):
    clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    data_str = request.GET.get('data_str')

    legacy_row = ClinVarLegacyRow.from_data_str(clinvar_key, data_str)
    matches: List[ClinVarLegacyMatches] = legacy_row.find_variant_grid_allele() if data_str else []

    # see if we match to a ClinVarExport, but not using SCV
    # implying that the SCV might need to be copied over
    has_clinvar_export = False
    has_scv_match = False
    for legacy_match in matches:
        for clinvar_export_match in legacy_match.clinvar_export_matches:
            has_clinvar_export = has_clinvar_export or bool(clinvar_export_match.clinvar_export)
            has_scv_match = has_scv_match or ClinVarLegacyExportMatchType.SCV_MATCHES in clinvar_export_match.match_types

    action_required = has_clinvar_export and not has_scv_match

    return render(request, 'classification/clinvar_match_detail.html', {
        'matches': matches,
        'action_required': action_required
    })
