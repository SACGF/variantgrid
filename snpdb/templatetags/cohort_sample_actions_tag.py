from django.template import Library, loader
from django.urls import reverse

from snpdb.archive import DataArchivedError
from snpdb.models import CohortGenotypeCollection, ImportStatus

register = Library()

# What a selection of samples can be turned into, by how many samples that takes
FAMILY_ANALYSES = [
    {"num_samples": 2, "label": "Duo", "url_name": "duo_wizard", "model": "snpdb.Duo"},
    {"num_samples": 3, "label": "Trio", "url_name": "trio_wizard", "model": "snpdb.Trio"},
    {"num_samples": 4, "label": "Quad", "url_name": "quad_wizard", "model": "snpdb.Quad"},
]


@register.simple_tag(takes_context=True)
def cohort_sample_actions(context, cohort):
    """ Bar for the samples ticked on the VCF / Cohort page - the count, and the one thing that many
        samples can launch @see js/sample_selection_actions.js """
    url_name_visible = context["url_name_visible"]
    sample_count = cohort.cohortsample_set.count()
    if cohort.import_status != ImportStatus.SUCCESS or sample_count < 2:
        return ""

    try:
        _ = cohort.cohort_genotype_collection  # Everything on offer here needs the packed genotypes
    except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        return ""

    family_analyses = []
    if url_name_visible["analysis"]:
        family_analyses = [fa for fa in FAMILY_ANALYSES
                           if url_name_visible[fa["url_name"]] and sample_count >= fa["num_samples"]]

    tag_context = {
        "cohort": cohort,
        "family_analyses": family_analyses,
        "sub_cohort_url": reverse("create_sub_cohort", kwargs={"cohort_id": cohort.pk}),
    }
    t = loader.get_template("snpdb/tags/cohort_sample_actions_tag.html")
    return t.render(tag_context)
