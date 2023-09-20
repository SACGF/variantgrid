import operator
from functools import reduce
from typing import Optional, List

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver

from annotation.models import VariantAnnotationVersion
from classification.models import Classification, ClassificationModification
from library.preview_request import preview_extra_signal, PreviewKeyValue
from ontology.models import OntologyTerm
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Lab, Organization, Allele, Variant, GenomeBuild
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Classification,
    example=SearchExample(
        note="The lab record ID",
        examples=["CR_1545" if settings.VARIANT_CLASSIFICATION_ID_OVERRIDE_PREFIX else "vc1545"]
    )
)
def classification_search(search_input: SearchInputInstance):

    search_string = search_input.search_string
    """ Search for LabId which can be either:
        "vc1080" or "Molecular Genetics/ Frome Road / vc1080" (as it appears in classification) or CR_1080 or Molecular Genetics/ Frome Road / CR_1080 """
    if search_string.startswith("CR_"):
        search_string = search_string[3:]
        filters = [Q(classification__id=search_string)]
    else:
        filters = [Q(classification__lab_record_id__iexact=search_string)]  # exact match
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
                if lab_record_id.startswith("CR_"):
                    lab_record_id = lab_record_id[3:]
                    filters.append(Q(classification__id=lab_record_id) & Q(classification__lab__in=lab_qs))
                else:
                    filters.append(Q(classification__lab_record_id=lab_record_id) & Q(classification__lab__in=lab_qs))

    q_cm = reduce(operator.or_, filters)
    cm_qs = ClassificationModification.filter_for_user(search_input.user).filter(is_last_published=True)
    cm_ids = cm_qs.filter(q_cm).values('classification')

    # check for source ID but only for labs that the user belongs to
    cm_source_ids = cm_qs.filter(classification__lab__in=Lab.valid_labs_qs(search_input.user, admin_check=True)).filter(
        classification__last_source_id=search_string).values('classification')

    # convert from modifications back to Classification so absolute_url returns the editable link
    yield Classification.objects.filter(Q(pk__in=cm_ids) | Q(pk__in=cm_source_ids))


def _allele_preview_classifications_extra(user: User, obj: Allele, genome_build: GenomeBuild) -> List[PreviewKeyValue]:
    cms = ClassificationModification.latest_for_user(user=user, allele=obj)
    extras = []
    hgvs_extras = []
    if count := cms.count():
        extras += [PreviewKeyValue.count(Classification, count)]
        column = ClassificationModification.column_name_for_build(genome_build)
        # provide the c.HGVS for alleles
        if c_hgvs := sorted(
                c_hgvs for c_hgvs in cms.order_by(column).values_list(column, flat=True).distinct().all() if c_hgvs):
            for hgvs in c_hgvs:
                hgvs_extras.append(PreviewKeyValue(None, hgvs, dedicated_row=True))

    if not hgvs_extras:
        try:
            v = obj.variant_for_build(genome_build)
            hgvs_extras = _variant_hgvs_extra(v, genome_build)
        except ValueError:
            pass

    return extras + hgvs_extras


def _variant_hgvs_extra(variant: Variant, genome_build: GenomeBuild) -> List[PreviewKeyValue]:
    hgvs_extras = []
    if c_hgvs := variant.get_canonical_c_hgvs(genome_build):
        hgvs_extras.append(PreviewKeyValue(None, c_hgvs, dedicated_row=True))
    else:
        vav = VariantAnnotationVersion.latest(genome_build)
        if va := variant.variantannotation_set.filter(version=vav).first():
            variant_summary = "intergenic"
            if va.distance:
                for direction in ["upstream", "downstream"]:
                    if direction in va.consequence:
                        variant_summary = f"{va.distance}bp {direction} of {va.symbol}"
                        break
            hgvs_extras.append(PreviewKeyValue(None, variant_summary, dedicated_row=True))

    return hgvs_extras


@receiver(preview_extra_signal, sender=Allele)
def allele_preview_classifications_extra(sender, user: User, obj: Allele, **kwargs):
    genome_build = GenomeBuildManager.get_current_genome_build()
    return _allele_preview_classifications_extra(user, obj, genome_build)


@receiver(preview_extra_signal, sender=Variant)
def variant_preview_classifications_extra(sender, user: User, obj: Variant, **kwargs):
    genome_build = next(iter(obj.genome_builds))
    if allele := obj.allele:
        return _allele_preview_classifications_extra(user, allele, genome_build)
    else:
        return _variant_hgvs_extra(obj, genome_build)


@receiver(preview_extra_signal, sender=OntologyTerm)
def ontology_preview_classifications_extra(sender, user: User, obj: OntologyTerm, **kwargs):
    terms = [{"term_id": obj.pk}]
    qs = ClassificationModification.latest_for_user(user=user, published=True, exclude_withdrawn=True, classification__condition_resolution__resolved_terms__contains=terms)
    if num_classifications := qs.count():
        return [PreviewKeyValue.count(Classification, num_classifications)]
