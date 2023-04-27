import operator
from functools import reduce
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver
from pyhgvs import HGVSName, InvalidHGVSName
from threadlocals.threadlocals import get_current_user

from annotation.models import VariantAnnotationVersion
from classification.models import Classification, ClassificationModification, ImportedAlleleInfo
from library.preview_request import preview_extra_signal, PreviewKeyValue
from ontology.models import OntologyTerm
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Lab, Organization, Allele
from snpdb.search import search_receiver, SearchInputInstance, SearchExample
from snpdb.user_settings_manager import UserSettingsManager


@search_receiver(
    search_type=Classification,
    example=SearchExample(
        note="The lab record ID",
        examples=["vc1545"]
    )
)
def classification_search(search_input: SearchInputInstance):

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

    # FIXME put this in its own search, or just remove it and the setting
    # maybe have variants report how many classifcations they have
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
    yield Classification.objects.filter(Q(pk__in=cm_ids) | Q(pk__in=cm_source_ids))


@receiver(preview_extra_signal, sender=Allele)
def allele_preview_classifications_extra(sender, user: User, obj: Allele, **kwargs):
    cms = ClassificationModification.latest_for_user(user=user, allele=obj)
    count = cms.count()
    extras = [PreviewKeyValue("Classification Count", count)]
    if count:
        genome_build = GenomeBuildManager.get_current_genome_build()
        column = ClassificationModification.column_name_for_build(genome_build)
        if c_hgvs := sorted(c_hgvs for c_hgvs in cms.order_by(column).values_list(column, flat=True).distinct().all() if c_hgvs):
            for hgvs in c_hgvs:
                extras.append(PreviewKeyValue(genome_build.name, hgvs, important=True))

    return extras


@receiver(preview_extra_signal, sender=OntologyTerm)
def ontology_preview_classifications_extra(sender, user: User, obj: OntologyTerm, **kwargs):
    terms = [{"term_id": obj.pk}]
    qs = Classification.filter_for_user(user).filter(condition_resolution__resolved_terms__contains=terms)
    if num_classifications := qs.count():
        return [PreviewKeyValue("Classification count", num_classifications)]
