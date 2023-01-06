from functools import cached_property
from typing import Optional, List, Union

from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.shortcuts import render, get_object_or_404
from guardian.shortcuts import get_objects_for_user
from rest_framework.response import Response
from rest_framework.views import APIView

from classification.models import ConditionTextMatch, ConditionText, update_condition_text_match_counts, MultiCondition, \
    ConditionMatchingSuggestion, condition_matching_suggestions
from classification.views.classification_dashboard_view import ClassificationDashboard
from library.utils import empty_to_none
from ontology.models import OntologyTerm
from ontology.ontology_matching import normalize_condition_text
from snpdb.lab_picker import LabPickerData
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class ConditionTextColumns(DatatableConfig):

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="normalized_text", label='Text', orderable=True, client_renderer="idRenderer", extra_columns=["id"], sort_keys=['normalized_text', 'lab__name']),
            RichColumn(key="classifications_count", label="Classification Count", orderable=True, default_sort=SortOrder.DESC, sort_keys=['classifications_count', '-normalized_text']),
            RichColumn(key="classifications_count_outstanding", label="Classifications Outstanding", default_sort=SortOrder.DESC, orderable=True, sort_keys=['classifications_count_outstanding', '-normalized_text'])
        ]
        if self.lab_picker.multi_labs_selected:
            self.rich_columns.insert(
                0,
                RichColumn(key="lab__name", label='Lab', orderable=True, sort_keys=['lab__name', 'normalized_text'])
            )

    @cached_property
    def lab_picker(self) -> LabPickerData:
        return LabPickerData.for_user(user=self.user, selection=self.get_query_param("lab_id"))

    def get_initial_queryset(self):
        # exclude where we've auto matched and have 0 outstanding left
        cts_qs: QuerySet[ConditionText] = get_objects_for_user(self.user, ConditionText.get_read_perm(), klass=ConditionText, accept_global_perms=True)
        cts_qs = cts_qs.filter(lab_id__in=self.lab_picker.lab_ids)
        return cts_qs

    def filter_queryset(self, qs: QuerySet) -> QuerySet:
        qs = qs.filter(classifications_count__gt=0)
        if filter_text := empty_to_none(self.get_query_param("text_filter")):
            normal_text = normalize_condition_text(filter_text)
            words = [word.strip() for word in normal_text.split(" ")]
            for word in words:
                qs = qs.filter(normalized_text__icontains=word)
        if empty_to_none(self.get_query_param("filter_outstanding")) == "true":
            qs = qs.exclude(classifications_count_outstanding=0)

        return qs


def condition_matchings_view(request, lab_id: Optional[Union[str,int]] = None):
    lab_picker = LabPickerData.from_request(request, lab_id, 'condition_matchings_lab')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, 'classification/condition_matchings.html', context={
        'dlab': ClassificationDashboard(lab_picker),
    })


def next_condition_text(current: ConditionText, user: User) -> Optional[ConditionText]:
    qs = ConditionText.objects.filter(classifications_count__lte=current.classifications_count, classifications_count_outstanding__gte=1).order_by('-classifications_count', 'normalized_text')
    qs = get_objects_for_user(user, ConditionText.get_read_perm(), qs, accept_global_perms=True)
    for ct in qs:
        if ct == current:
            continue
        if ct.classifications_count == current.classifications_count and ct.normalized_text < current.normalized_text:
            continue
        return ct
    return None


def condition_matching_view(request, pk: int):
    ct: ConditionText = get_object_or_404(ConditionText.objects.filter(pk=pk))
    ct.check_can_view(request.user)

    root_match = get_object_or_404(ConditionTextMatch.objects.filter(condition_text=ct, gene_symbol__isnull=True))

    return render(request, 'classification/condition_matching.html', context={
        'condition_text': ct,
        'condition_match': root_match,
        'next_condition_text': next_condition_text(current=ct, user=request.user)
    })


class ConditionTextMatchingAPI(APIView):

    def get(self, request, **kwargs) -> Response:
        ct_id = kwargs['pk']
        ct: ConditionText = ConditionText.objects.get(pk=ct_id)
        user: User = request.user
        ct.check_can_view(user)

        suggestions = condition_matching_suggestions(ct)
        return Response({
            "count_outstanding": ct.classifications_count_outstanding,
            "suggestions": [cms.as_json() for cms in suggestions]
        })

    def post(self, request, **kwargs) -> Response:
        ct_id = kwargs['pk']
        ct: ConditionText = ConditionText.objects.get(pk=ct_id)
        user: User = request.user
        ct.check_can_write(user)

        data = request.data.get("changes")
        ctm: Optional[ConditionTextMatch] = None
        errors: List[str] = []
        for update in data:
            terms = update.get('terms')
            valid_terms: List[str] = []
            if terms:
                terms = [term.strip() for term in terms]
                for term in terms:
                    if term:
                        try:
                            om = OntologyTerm.get_or_stub(term)
                            valid_terms.append(om.id)
                        except ValueError:
                            errors.append(f'"{term}" is not a valid ontology term.')

        if not errors:
            ctm = ct.conditiontextmatch_set.get(pk=update.get('ctm_id'))
            ctm.condition_xrefs = valid_terms
            ctm.condition_multi_operation = MultiCondition(update.get('joiner') or MultiCondition.NOT_DECIDED)
            ctm.last_edited_by = request.user
            ctm.save()

            update_condition_text_match_counts(ct)
            ct.save()

        if len(data) == 1 and ctm and not ctm.is_root and not (ctm.is_gene_level and not ctm.condition_xrefs):
            # if we only have 1 suggestion, and it's not the root, and isn't a blanked out record with empty conditions
            cms = ConditionMatchingSuggestion(ctm)
            cms.validate()
            return Response({
                "errors": errors,
                "suggestions": [cms.as_json()],
                "count_outstanding": ct.classifications_count_outstanding
            })
        else:
            suggestions = condition_matching_suggestions(ct)
            return Response({
                "errors": errors,
                "suggestions": [cms.as_json() for cms in suggestions],
                "count_outstanding": ct.classifications_count_outstanding
            })
