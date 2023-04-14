import operator
from functools import reduce
from typing import Any, Optional, Type

from django.conf import settings
from django.db.models import Q
from django.dispatch import receiver
from pyhgvs import HGVSName, InvalidHGVSName

from annotation.models import VariantAnnotationVersion
from classification.models import Classification, ClassificationModification
from snpdb.models import Lab, Organization
from snpdb.search2 import search_signal, SearchInput, SearchResponse


@receiver(search_signal, sender=SearchInput)
def classification_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response = SearchResponse("Classification")

    search_string = search_input.search_string
    """ Search for LabId which can be either:
        "vc1080" or "Molecular Genetics, Frome Road / vc1080" (as it appears in classification) """

    filters = [Q(classification__lab_record_id=search_string)]  # exact match
    slash_index = search_string.find("/")
    if slash_index > 0:
        parts = [p.strip() for p in search_string.split("/")]
        org_name: Optional[str] = None
        lab_name: Optional[str] = None
        lab_record_id: str = parts[-1]

        if len(parts) >= 2:
            lab_name = parts[-2]
        if len(parts) >= 3:
            org_name = parts[-3]

        if lab_name:
            org: Optional[Organization] = None
            if org_name:
                org = Organization.objects.filter(Q(name=org_name) | Q(short_name=org_name) | Q(group_name=org_name)).first()

            lab_qs = Lab.objects.filter(Q(name=lab_name) | Q(group_name__endswith="/" + lab_name))
            if org:
                lab_qs = lab_qs.filter(organization=org)

            if lab_qs:
                filters.append(Q(classification__lab_record_id=lab_record_id) & Q(classification__lab__in=lab_qs))

    # We can also filter for c.HGVS in classifications
    if settings.SEARCH_CLASSIFICATION_HGVS_SUFFIX:
        if search_string.startswith("c.") or search_string.startswith("n."):
            cdna = search_string[2:]
            try:
                hgvs_name = HGVSName()
                hgvs_name.parse_cdna(cdna)  # Valid HGVS cdna section

                # Look for direct string in evidence
                filters.append(Q(classification__evidence__c_hgvs__icontains=search_string))

                # Look for HGVS in transcript annotation pointing to classified allele
                vav = VariantAnnotationVersion.latest(search_input.genome_build_preferred)
                va_path = "classification__allele__variantallele"
                vta_path = f"{va_path}__variant__varianttranscriptannotation"
                filters.append(Q(**{f"{va_path}__genome_build": search_input.genome_build_preferred,
                                    f"{vta_path}__version": vav,
                                    f"{vta_path}__hgvs_c__endswith": search_string}))

            except InvalidHGVSName:
                pass

    q_cm = reduce(operator.or_, filters)
    cm_qs = ClassificationModification.filter_for_user(search_input.user).filter(is_last_published=True)
    cm_ids = cm_qs.filter(q_cm).values('classification')

    # check for source ID but only for labs that the user belongs to
    cm_source_ids = cm_qs.filter(classification__lab__in=Lab.valid_labs_qs(search_input.user, admin_check=True)).filter(
        classification__last_source_id=search_string).values('classification')

    # convert from modifications back to Classification so absolute_url returns the editable link
    qs = Classification.objects.filter(Q(pk__in=cm_ids) | Q(pk__in=cm_source_ids))
    response.extend(qs)

    return response
