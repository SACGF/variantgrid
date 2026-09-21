from collections import defaultdict
from functools import reduce
from operator import itemgetter, or_

from django.db.models import Q
from django.template import Library

from classification.models import Classification
from classification.views.classification_datatables import ClassificationColumns
from pedigree.models import CohortSamplePedFileRecord, Pedigree
from snpdb.models import CohortSample, Duo, Quad, Trio
from snpdb.models.models_enums import ImportStatus

register = Library()

TRIO_SAMPLES_SELECT_RELATED = ("mother__sample", "father__sample", "proband__sample")


def related_data_context(context, samples):
    tag_context = {
        "samples": samples,
        "user": context["user"]
    }
    classifications = Classification.objects.filter(sample__in=samples)
    if classifications.exists():
        tag_context["has_classifications"] = True
        tag_context["sample_ids_list"] = [s.pk for s in samples]
        tag_context["datatable_config"] = ClassificationColumns(context["request"])
    return tag_context


@register.inclusion_tag("snpdb/tags/related_data_for_patient.html", takes_context=True)
def related_data_for_patient(context, patient):
    context = related_data_context(context, patient.get_samples())
    context["cases"] = patient.case_set.all()
    return context


def summarise_samples(x_and_samples_list, show_sample_info):
    x_and_samples = []
    for obj, samples_list in sorted(x_and_samples_list.items(), key=itemgetter(0)):
        samples_info = None
        if show_sample_info:
            samples_info = ','.join(sorted(samples_list))
        x_and_samples.append((obj, samples_info))
    return x_and_samples


@register.inclusion_tag("snpdb/tags/related_data_for_samples.html", takes_context=True)
def related_data_for_samples(context, samples, show_sample_info=True):
    cohorts_and_samples_list = defaultdict(list)
    trios_and_samples_list = defaultdict(list)
    quads_and_samples_list = defaultdict(list)
    duos_and_samples_list = defaultdict(list)
    pedigrees_and_samples_list = defaultdict(list)

    cohort_samples = list(CohortSample.objects.filter(sample__in=samples).select_related("cohort", "sample"))
    successful_cohort_samples = [cs for cs in cohort_samples
                                 if cs.cohort.import_status == ImportStatus.SUCCESS]
    for cs in successful_cohort_samples:
        cohorts_and_samples_list[cs.cohort].append(cs.sample.name)

    if successful_cohort_samples:
        cs_ids = {cs.pk for cs in successful_cohort_samples}
        for family_class, family_and_samples_list in [(Trio, trios_and_samples_list),
                                                      (Quad, quads_and_samples_list),
                                                      (Duo, duos_and_samples_list)]:
            member_q = reduce(or_, [Q(**{f"{field}__in": cs_ids}) for field in family_class.MEMBER_FIELDS])
            member_samples = [f"{field}__sample" for field in family_class.MEMBER_FIELDS]
            for family in family_class.objects.filter(member_q).select_related(*member_samples):
                for cs in family.get_cohort_samples():
                    if cs.pk in cs_ids:
                        family_and_samples_list[family].append(cs.sample.name)

    if cohort_samples:
        seen_pedigree_cohort_samples = set()
        record_qs = CohortSamplePedFileRecord.objects.filter(cohort_sample__in=cohort_samples)
        for record in record_qs.select_related("pedigree", "cohort_sample__sample"):
            pair = (record.pedigree_id, record.cohort_sample_id)
            if pair not in seen_pedigree_cohort_samples:
                seen_pedigree_cohort_samples.add(pair)
                pedigrees_and_samples_list[record.pedigree].append(record.cohort_sample.sample.name)

    cohorts_and_samples = summarise_samples(cohorts_and_samples_list, show_sample_info)
    trios_and_samples = summarise_samples(trios_and_samples_list, show_sample_info)
    quads_and_samples = summarise_samples(quads_and_samples_list, show_sample_info)
    duos_and_samples = summarise_samples(duos_and_samples_list, show_sample_info)
    pedigrees_and_samples = summarise_samples(pedigrees_and_samples_list, show_sample_info)

    context = related_data_context(context, samples)
    context.update({"show_sample_info": show_sample_info,
                    "cohorts_and_samples": cohorts_and_samples,
                    "trios_and_samples": trios_and_samples,
                    "quads_and_samples": quads_and_samples,
                    "duos_and_samples": duos_and_samples,
                    "pedigrees_and_samples": pedigrees_and_samples})
    return context


@register.inclusion_tag("snpdb/tags/related_data_for_cohort.html", takes_context=True)
def related_data_for_cohort(context, cohort):
    context = related_data_context(context, cohort.get_samples())
    # Trios/pedigrees are often built off a sub cohort (@see VCF page) - show them on the parent too
    sub_cohorts = list(cohort.sub_cohort_set.all())
    cohorts = [cohort] + sub_cohorts
    context["cohort"] = cohort
    context["sub_cohorts"] = sub_cohorts
    context["trios"] = list(Trio.objects.filter(cohort__in=cohorts).select_related(*TRIO_SAMPLES_SELECT_RELATED))
    context["quads"] = list(Quad.objects.filter(cohort__in=cohorts))
    context["duos"] = list(Duo.objects.filter(cohort__in=cohorts))
    context["pedigrees"] = list(Pedigree.objects.filter(cohort__in=cohorts))
    return context


@register.inclusion_tag("snpdb/tags/related_data_for_family.html", takes_context=True)
def related_data_for_family(context, family):
    """ Duo/Trio/Quad - anything with FamilyGroupMixin.get_samples() """
    return related_data_context(context, family.get_samples())


@register.inclusion_tag("snpdb/tags/related_data_for_pedigree.html")
def related_data_for_pedigree(pedigree):
    return {"pedigree": pedigree}
