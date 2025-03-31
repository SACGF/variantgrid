from dataclasses import dataclass
from enum import IntEnum
from functools import cached_property, reduce
from typing import Iterable, Optional

from django.db.models import Count
from django.http import HttpRequest
from django.shortcuts import render
from more_itertools import first

from classification.models import ClinVarAllele, ClinVarExport, ConditionResolved, ClassificationModification
from library.django_utils import require_superuser
from ontology.models import AncestorCalculator, OntologyTerm


class ClinVarAlleleMultiMergeStatus(IntEnum):
    CAN_MERGE = 1
    DELETE_ALL = 2
    MULTIPLE_SCVS = 3
    MULTIPLE_SUBMISSIONS = 4
    NO_COMMON_CONDITION = 5

    @property
    def can_action(self) -> bool:
        return self in {ClinVarAlleleMultiMergeStatus, ClinVarAlleleMultiMergeStatus.DELETE_ALL}

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
    preview: Optional[ClinVarAlleleMultiMergePreview]
    status: ClinVarAlleleMultiMergeStatus


class ClinVarAlleleMultiExport:

    def __init__(self, clinvar_allele: ClinVarAllele):
        self.clinvar_allele = clinvar_allele

    @cached_property
    def clinvar_exports(self) -> list[ClinVarExport]:
        return list(sorted(self.clinvar_allele.clinvarexport_set.all()))

    @cached_property
    def common_ancestor(self) -> Optional[OntologyTerm]:
        try:
            return AncestorCalculator.common_ancestor(terms=reduce(lambda x,y: x+y, (export.condition_resolved.terms for export in self.clinvar_exports), []))
        except ValueError:
            return None

    @cached_property
    def merge_preview(self) -> ClinVarAlleleMultiMergeOutput:
        consider_exports = list(self.clinvar_exports)

        # filter out records that are completely pointless
        def is_something(clinvar_export: ClinVarExport):
            return bool(clinvar_export.scv) or bool(clinvar_export.classification_based_on) \
                or clinvar_export.has_submission

        consider_exports = [c for c in consider_exports if is_something(c)]
        if not consider_exports:
            return ClinVarAlleleMultiMergeOutput(preview=None, status=ClinVarAlleleMultiMergeStatus.DELETE_ALL)

        all_scvs = {c.scv for c in consider_exports if c.scv}
        if len(all_scvs) > 1:
            return ClinVarAlleleMultiMergeOutput(preview=None, status=ClinVarAlleleMultiMergeStatus.MULTIPLE_SCVS)

        exports_with_submission_count = len([c for c in consider_exports if c.has_submission])
        if exports_with_submission_count > 1:
            return ClinVarAlleleMultiMergeOutput(preview=None, status=ClinVarAlleleMultiMergeStatus.MULTIPLE_SUBMISSIONS)

        conditions_that_matter = [c.condition_resolved for c in consider_exports if bool(c.classification_based_on)]
        resulting_condition: Optional[ConditionResolved] = None

        if len(conditions_that_matter) == 1:
            resulting_condition = conditions_that_matter[0]
        elif len(conditions_that_matter) > 1:
            resulting_condition = self.common_ancestor

        if resulting_condition is None:
            return ClinVarAlleleMultiMergeOutput(preview=None, status=ClinVarAlleleMultiMergeStatus.NO_COMMON_CONDITION)

        use_scv = None
        if len(all_scvs) == 1:
            use_scv = first(all_scvs)

        merge_submissions = exports_with_submission_count > 0

        currently_linked_classification_modifications = [c.classification_based_on for c in consider_exports if bool(c.classification_based_on)]
        link_classification: ClassificationModification = None
        if len(currently_linked_classification_modifications) > 0:
            link_classification = max(currently_linked_classification_modifications, key=lambda c: c.curated_date_check)

        pk = min(c.pk for c in consider_exports)

        return ClinVarAlleleMultiMergeOutput(
            preview=ClinVarAlleleMultiMergePreview(
                pk=pk,
                scv=use_scv,
                merge_submissions=merge_submissions,
                condition=resulting_condition,
                classification_based_on=link_classification
            ),
            status=ClinVarAlleleMultiMergeStatus.CAN_MERGE
        )


def multi_export_clinvar_alleles() -> Iterable[ClinVarAlleleMultiExport]:
    return [ClinVarAlleleMultiExport(ca) for ca in (ClinVarAllele.objects.annotate(clinvar_export_count=Count('clinvarexport')).filter(clinvar_export_count__gte=2).iterator())]


@require_superuser
def view_multi_clinvar_exports_listing(request: HttpRequest):
    return render(
        request=request,
        context={"clinvar_allele_mults": multi_export_clinvar_alleles()},
        template_name="classification/clinvar_export_multi_listing.html"
    )


@require_superuser
def view_multi_clinvar_exports(request: HttpRequest, clinvar_allele_pk: int):
    return render(
        request=request,
        context={"clinvar_allele_mult": ClinVarAlleleMultiExport(ClinVarAllele.objects.get(pk=clinvar_allele_pk))},
        template_name="classification/clinvar_export_multi.html"
    )