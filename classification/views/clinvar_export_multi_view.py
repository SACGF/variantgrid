from dataclasses import dataclass
from enum import IntEnum
from functools import cached_property, reduce
from typing import Iterable, Optional, Any

from django.contrib import messages
from django.db import transaction
from django.db.models import Count, Q, QuerySet
from django.http import HttpRequest
from django.shortcuts import render, redirect
from django.urls import reverse
from more_itertools import first

from classification.models import ClinVarAllele, ClinVarExport, ConditionResolved, ClassificationModification, \
    ClinVarExportDeleteStatus, ClinVarAlleleMultiStatus
from classification.models.clinvar_export_convertor import ClinVarExportConverter
from library.django_utils import require_superuser
from library.utils import JsonObjType
from ontology.models import AncestorCalculator, OntologyTerm
from snpdb.views.datatable_view import DatatableConfig, DC, DatatableConfigQuerySetMode, RichColumn, CellData, SortOrder


class ClinVarAlleleMultiMergeStatus(IntEnum):
    CAN_MERGE = 1
    DELETE_ALL = 2
    MULTIPLE_SCVS = 3
    MULTIPLE_SUBMISSIONS = 4
    NO_COMMON_CONDITION = 5

    @property
    def can_action(self) -> bool:
        return self in {ClinVarAlleleMultiMergeStatus.CAN_MERGE, ClinVarAlleleMultiMergeStatus.DELETE_ALL}

    @property
    def display(self) -> str:
        match self:
            case ClinVarAlleleMultiMergeStatus.CAN_MERGE: return "Can Merge"
            case ClinVarAlleleMultiMergeStatus.DELETE_ALL: return "Delete All Records (none are worth keeping)"
            case ClinVarAlleleMultiMergeStatus.MULTIPLE_SCVS: return "Cannot merge, multiple unique SCVs"
            case ClinVarAlleleMultiMergeStatus.MULTIPLE_SUBMISSIONS: return "Cannot merge, multiple records have already been submitted to ClinVar"
            case ClinVarAlleleMultiMergeStatus.NO_COMMON_CONDITION: return "No Common Condition"
            case _: return self.name


@dataclass
class ClinVarAlleleMultiMergePreview:
    pk: int
    scv: Optional[str]
    merge_submissions: bool
    condition: ConditionResolved
    classification_based_on: ClassificationModification


@dataclass
class ClinVarAlleleMultiMergeOutput:
    clinvar_allele: ClinVarAllele
    preview: Optional[ClinVarAlleleMultiMergePreview]
    status: ClinVarAlleleMultiMergeStatus

    @property
    def can_action(self) -> bool:
        return self.status.can_action

    @transaction.atomic()
    def apply(self) -> list[str]:
        messages = []
        delete_us: list[ClinVarExport]
        if self.status == ClinVarAlleleMultiMergeStatus.DELETE_ALL:
            delete_us = list(self.clinvar_allele.clinvarexport_set.exclude(delete_status=ClinVarExportDeleteStatus.DELETED))
        elif self.status == ClinVarAlleleMultiMergeStatus.CAN_MERGE:
            merge_into = self.clinvar_allele.clinvarexport_set.filter(pk=self.preview.pk).get()
            merge_into.scv = self.preview.scv or ""
            merge_into.condition = self.preview.condition.to_json()
            merge_into.classification_based_on = self.preview.classification_based_on
            ClinVarExportConverter(merge_into).convert().apply()
            merge_into.save()
            messages.append(f"Updated CE_{merge_into.pk}")
            delete_us = list(self.clinvar_allele.clinvarexport_set.exclude(pk=self.preview.pk).exclude(delete_status=ClinVarExportDeleteStatus.DELETED))

        for delete_me in delete_us:
            messages.append(f"Deleted CE_{delete_me.pk}")
            delete_me.delete()

        return messages


class ClinVarAlleleMultiExport:

    def __init__(self, clinvar_allele: ClinVarAllele):
        self.clinvar_allele = clinvar_allele

    @cached_property
    def clinvar_exports(self) -> list[ClinVarExport]:
        return list(sorted(self.clinvar_allele.clinvarexport_set.exclude(delete_status=ClinVarExportDeleteStatus.DELETED)))

    @cached_property
    def clinvar_exports_with_classification(self) -> list[ClinVarExport]:
        return list(sorted(self.clinvar_allele.clinvarexport_set.exclude(delete_status=ClinVarExportDeleteStatus.DELETED).exclude(classification_based_on__isnull=True)))

    @cached_property
    def common_ancestor(self) -> OntologyTerm:
        common_ancestor: OntologyTerm
        try:
            common_ancestor = AncestorCalculator.common_ancestor(terms=reduce(lambda x,y: x+y, (export.condition_resolved.terms for export in self.clinvar_exports_with_classification), []))
        except ValueError:
            pass
        return OntologyTerm.get_or_stub("MONDO:0000001")

    @property
    def common_ancestor_json(self) -> JsonObjType:
        return ConditionResolved(terms=[self.common_ancestor]).to_json(include_join=False)

    def approve(self):
        all_conditions = set()
        for clinvar_export in self.clinvar_exports_with_classification:
            if resolved := clinvar_export.condition_resolved:
                all_conditions.add(resolved)
        self.clinvar_allele.accepted_condition_splits = list(sorted(all_conditions))
        self.clinvar_allele.multiple_review_status = ClinVarAlleleMultiStatus.MULTI_CONDITION_APPROVED
        self.clinvar_allele.save()

    def unapprove(self):
        self.clinvar_allele.accepted_condition_splits = []
        self.clinvar_allele.multiple_review_status = ClinVarAlleleMultiStatus.MULTI_CONDITION_REQUIRES_REVIEW
        self.clinvar_allele.save()


    @cached_property
    def merge_preview(self) -> ClinVarAlleleMultiMergeOutput:
        consider_exports = list(self.clinvar_exports)

        # filter out records that are completely pointless
        def is_something(clinvar_export: ClinVarExport):
            return bool(clinvar_export.scv) or bool(clinvar_export.classification_based_on) \
                or clinvar_export.has_submission

        consider_exports = [c for c in consider_exports if is_something(c)]
        if not consider_exports:
            return ClinVarAlleleMultiMergeOutput(clinvar_allele=self.clinvar_allele, preview=None, status=ClinVarAlleleMultiMergeStatus.DELETE_ALL)

        all_scvs = {c.scv for c in consider_exports if c.scv}
        if len(all_scvs) > 1:
            return ClinVarAlleleMultiMergeOutput(clinvar_allele=self.clinvar_allele, preview=None, status=ClinVarAlleleMultiMergeStatus.MULTIPLE_SCVS)

        exports_with_submissions = [c for c in consider_exports if c.has_submission]
        exports_with_submission_count = len(exports_with_submissions)
        if exports_with_submission_count > 1:
            return ClinVarAlleleMultiMergeOutput(clinvar_allele=self.clinvar_allele, preview=None, status=ClinVarAlleleMultiMergeStatus.MULTIPLE_SUBMISSIONS)

        conditions_that_matter = [c.condition_resolved for c in consider_exports if bool(c.classification_based_on)]
        resulting_condition: Optional[ConditionResolved] = None

        if len(conditions_that_matter) == 1:
            resulting_condition = conditions_that_matter[0]
        elif len(conditions_that_matter) > 1:
            resulting_condition = ConditionResolved(terms=[self.common_ancestor])

        if resulting_condition is None:
            return ClinVarAlleleMultiMergeOutput(clinvar_allele=self.clinvar_allele, preview=None, status=ClinVarAlleleMultiMergeStatus.NO_COMMON_CONDITION)

        use_scv = None
        if len(all_scvs) == 1:
            use_scv = first(all_scvs)

        merge_submissions = exports_with_submission_count > 0

        currently_linked_classification_modifications = [c.classification_based_on for c in consider_exports if bool(c.classification_based_on)]
        link_classification: ClassificationModification = None
        if len(currently_linked_classification_modifications) > 0:
            link_classification = max(currently_linked_classification_modifications, key=lambda c: c.curated_date_check)

        pk: int
        if exports_with_submissions:
            pk = exports_with_submissions[0].pk
        else:
            pk = min(c.pk for c in consider_exports)

        return ClinVarAlleleMultiMergeOutput(
            clinvar_allele=self.clinvar_allele,
            preview=ClinVarAlleleMultiMergePreview(
                pk=pk,
                scv=use_scv,
                merge_submissions=merge_submissions,
                condition=resulting_condition,
                classification_based_on=link_classification
            ),
            status=ClinVarAlleleMultiMergeStatus.CAN_MERGE
        )


class ClinVarAlleleMultiExportColumns(DatatableConfig[ClinVarAllele]):

    def get_initial_queryset(self) -> QuerySet[DC]:
        return ClinVarAllele.objects.filter(multiple_review_status__in=(ClinVarAlleleMultiStatus.MULTI_CONDITION_REQUIRES_REVIEW, ClinVarAlleleMultiStatus.MULTI_CONDITION_APPROVED))

    def map_object(self, obj: DC) -> Any:
        return ClinVarAlleleMultiExport(obj)

    def pre_render(self, qs: QuerySet[DC]):
        return qs.select_related('clinvar_key', 'allele').prefetch_related('clinvarexport_set')

    def render_exports(self, row: CellData):
        multi: ClinVarAlleleMultiExport = row.obj
        data = []
        for ce in multi.clinvar_exports_with_classification:
            data.append({"pk": ce.pk, "status": ce.get_status_display(), "condition": ce.condition})
        return data

    def __init__(self, request: HttpRequest):
        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="pk", label="ID", renderer=lambda cell: {"text": cell.obj.clinvar_allele.pk, "url": reverse('clinvar_export_multi', kwargs={"clinvar_allele_pk": cell.obj.clinvar_allele.pk})}, client_renderer='TableFormat.linkUrl', orderable=True),
            RichColumn(key="clinvar_key", label="ClinVar Key", renderer=lambda cell: str(cell.obj.clinvar_allele.clinvar_key)),
            RichColumn(key="allele", label="Allele", renderer=lambda cell: str(cell.obj.clinvar_allele.allele), sort_keys=['allele']),
            RichColumn(key="export_bucket", label="Export Type", renderer=lambda cell: str(cell.obj.clinvar_allele.clinvar_export_bucket), sort_keys=['clinvar_allele__clinvar_export_bucket'], client_renderer='VCTable.allele_origin_cell'),
            RichColumn(name="common-condition", label="Common Condition", renderer=lambda cell: cell.obj.common_ancestor_json, client_renderer='VCTable.condition'),
            RichColumn(name="exports", label="ClinVar Exports", renderer=self.render_exports, client_renderer='renderClinVarExports'),
            RichColumn(name="status", label="Duplicate Status", renderer=lambda cell: cell.obj.clinvar_allele.get_multiple_review_status_display(), sort_keys=['multiple_review_status'], orderable=True, default_sort=SortOrder.DESC)
        ]



@require_superuser
def view_multi_clinvar_exports_listing(request: HttpRequest):
    return render(
        request=request,
        template_name="classification/clinvar_export_multi_listing.html"
    )


@require_superuser
def view_multi_clinvar_exports(request: HttpRequest, clinvar_allele_pk: int):
    if request.method == "POST":
        action = request.POST.get("action")
        merge_me = ClinVarAlleleMultiExport(ClinVarAllele.objects.get(pk=clinvar_allele_pk))
        if action == "merge":
            message_list = merge_me.merge_preview.apply()
            for message in message_list:
                messages.add_message(request, messages.SUCCESS, message)
        elif action == "approve":
            merge_me.approve()
            messages.add_message(request, messages.SUCCESS, f"ClinVarAllele {clinvar_allele_pk} marked as accepting provided conditions")
        elif action == "unapprove":
            merge_me.unapprove()
            messages.add_message(request, messages.SUCCESS,f"ClinVarAllele {clinvar_allele_pk} marked as accepting provided conditions")
        else:
            raise ValueError(f"Unknown action {action}")

        return redirect(reverse('clinvar_export_multi', kwargs={"clinvar_allele_pk": clinvar_allele_pk}))
    else:
        return render(
            request=request,
            context={"clinvar_allele_mult": ClinVarAlleleMultiExport(ClinVarAllele.objects.get(pk=clinvar_allele_pk))},
            template_name="classification/clinvar_export_multi.html"
        )