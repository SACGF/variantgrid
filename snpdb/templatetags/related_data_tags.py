from collections import defaultdict
from django.template import Library
from operator import itemgetter

from snpdb.models.models_enums import ImportStatus
from classification.models import VariantClassification
from classification.views.variant_classification_datatables import VariantClassificationDatatableConfig

register = Library()


def related_data_context(context, samples):
    tag_context = {
        "samples": samples,
        "user": context["user"]
    }
    variant_classifications = VariantClassification.objects.filter(sample__in=samples)
    if variant_classifications.exists():
        tag_context["has_variant_classifications"] = True
        tag_context["sample_ids_list"] = [s.pk for s in samples]
        tag_context["datatable_config"] = VariantClassificationDatatableConfig(context["request"])
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
    pedigrees_and_samples_list = defaultdict(list)

    for sample in samples:
        for cs in sample.cohortsample_set.all():
            cohort = cs.cohort
            if cohort.import_status == ImportStatus.SUCCESS:
                cohorts_and_samples_list[cohort].append(sample.name)

                for trio in cohort.trio_set.all():
                    trios_and_samples_list[trio].append(sample.name)

            for pedigree in cohort.pedigree_set.all():
                pedigrees_and_samples_list[pedigree].append(sample.name)

    cohorts_and_samples = summarise_samples(cohorts_and_samples_list, show_sample_info)
    trios_and_samples = summarise_samples(trios_and_samples_list, show_sample_info)
    pedigrees_and_samples = summarise_samples(pedigrees_and_samples_list, show_sample_info)

    context = related_data_context(context, samples)
    context.update({"show_sample_info": show_sample_info,
                    "cohorts_and_samples": cohorts_and_samples,
                    "trios_and_samples": trios_and_samples,
                    "pedigrees_and_samples": pedigrees_and_samples})
    return context


@register.inclusion_tag("snpdb/tags/related_data_for_cohort.html", takes_context=True)
def related_data_for_cohort(context, cohort):
    context = related_data_context(context, cohort.get_samples())
    context["cohort"] = cohort
    return context


@register.inclusion_tag("snpdb/tags/related_data_for_trio.html", takes_context=True)
def related_data_for_trio(context, trio):
    return related_data_context(context, trio.cohort.get_samples())


@register.inclusion_tag("snpdb/tags/related_data_for_pedigree.html")
def related_data_for_pedigree(pedigree):
    return {"pedigree": pedigree}
