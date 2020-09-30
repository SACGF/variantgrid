from django.template.library import Library

from snpdb.templatetags.related_data_tags import related_data_context

register = Library()


@register.inclusion_tag("pathtests/tags/related_data_for_case.html", takes_context=True)
def related_data_for_case(context, case):
    samples = case.patient.get_samples()
    context.update(related_data_context(context, samples))
    context["pathology_test_orders"] = case.pathologytestorder_set.all()
    context["samples"] = samples
    return context


@register.inclusion_tag("pathtests/tags/related_data_for_pathology_test_order.html", takes_context=True)
def related_data_for_pathology_test_order(context, pathology_test_order):
    samples = pathology_test_order.case.patient.get_samples()
    context.update(related_data_context(context, samples))
    context["samples"] = samples
    try:
        context["sapath_link"] = pathology_test_order.sapathologyrequestgenelistpathologytestorderlink
    except:  # May not have SA Path app installed
        pass
    return context
