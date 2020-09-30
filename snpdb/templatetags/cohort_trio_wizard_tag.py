from django.template import Library, loader

register = Library()


@register.simple_tag(takes_context=True)
def cohort_trio_wizard(context, cohort):
    cohort.cohort_genotype_collection  # Ensure it's loaded

    context = {"url_name_visible": context["url_name_visible"],
               'cohort': cohort,
               'sample_count': cohort.cohortsample_set.count()}
    t = loader.get_template("snpdb/tags/cohort_trio_wizard_tag.html")
    return t.render(context)
