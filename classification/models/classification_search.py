import operator
from functools import reduce
from typing import Any, Optional, Union
from django.db.models import Q
from django.dispatch import receiver
from django.utils.safestring import SafeString

from classification.enums import SpecialEKeys
from classification.models import Classification, ClassificationModification, EvidenceKeyMap
from snpdb.models import Lab
from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse


class SearchResponseClassification(SearchResponseRecordAbstract[Classification]):

    @classmethod
    def search_type(cls) -> str:
        return "Classification"

    def display(self) -> Union[str, SafeString]:
        last_published = self.record.last_published_version
        classification: Classification = self.record
        clin_sig = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).\
            pretty_value(last_published.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)) or 'Unclassified'
        return f"{clin_sig} {classification.lab.name} / {classification.lab_record_id}"


@receiver(search_signal, sender=SearchInput)
def search_ontology(sender: Any, search_input: SearchInput, **kwargs) -> Optional[SearchResponse]:
    response: SearchResponse[SearchResponseClassification] = SearchResponse(SearchResponseClassification)

    search_string = search_input.search_string
    """ Search for LabId which can be either:
        "vc1080" or "Molecular Genetics, Frome Road / vc1080" (as it appears in classification) """

    filters = [Q(classification__lab_record_id=search_string)]  # exact match
    slash_index = search_string.find("/")
    if slash_index > 0:
        lab_name = search_string[:slash_index].strip()
        lab_record_id = search_string[slash_index + 1:].strip()
        if lab_name and lab_record_id:
            q_lab_name = Q(classification__lab__name=lab_name)
            q_lab_record = Q(classification__lab_record_id=lab_record_id)
            filters.append(q_lab_name & q_lab_record)
    q_vcm = reduce(operator.or_, filters)
    vcm_qs = ClassificationModification.filter_for_user(search_input.user).filter(is_last_published=True)
    vcm_ids = vcm_qs.filter(q_vcm).values('classification')

    # check for source ID but only for labs that the user belongs to
    vcm_source_ids = vcm_qs.filter(classification__lab__in=Lab.valid_labs_qs(search_input.user, admin_check=True)).filter(
        published_evidence__source_id__value=search_string).values('classification')

    # convert from modifications back to Classification so absolute_url returns the editable link
    qs = Classification.objects.filter(Q(pk__in=vcm_ids) | Q(pk__in=vcm_source_ids))
    response.extend(SearchResponseClassification.from_iterable(qs))

    return response
