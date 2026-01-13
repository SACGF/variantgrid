import operator
import uuid
from collections import defaultdict
from functools import reduce

from django.conf import settings
from django.db.models import Q, Model
from django.template import Library

from analysis.forms import get_analysis_template_form_for_variables_only_of_class
from analysis.models import MutationalSignature, Analysis, AnalysisTemplate
from analysis.models.models_karyomapping import KaryomappingAnalysis
from analysis.related_analyses import get_related_analysis_details_for_samples, \
    get_related_analysis_details_for_cohort, \
    get_related_analysis_details_for_pedigree, get_related_analysis_details_for_trio
from pedigree.models import Pedigree
from snpdb.models import Cohort, Trio

register = Library()


def get_all_analyses_for_user(user, samples, cohorts=None, trios=None, pedigrees=None):
    analysis_details: dict[Analysis, list] = defaultdict(list)
    if pedigrees:
        for analysis, details in get_related_analysis_details_for_pedigree(user, pedigrees):
            analysis_details[analysis].append(f'Pedigree: {details}')

    if cohorts:
        for analysis, details in get_related_analysis_details_for_cohort(user, cohorts):
            analysis_details[analysis].append(f'Cohort {details}')

    if trios:
        for analysis, details in get_related_analysis_details_for_trio(user, trios):
            analysis_details[analysis].append(f'Trio: {details}')

    for analysis, details in get_related_analysis_details_for_samples(user, samples):
        analysis_details[analysis].append(f"Sample: {details}")

    return [(analysis, ", ".join(details)) for analysis, details in analysis_details.items()]


def update_context_with_related_analysis(context, samples, cohorts=None, trios=None,
                                         pedigrees=None, show_sample_info=True):
    user = context["user"]
    analysis_details = get_all_analyses_for_user(user, samples, cohorts=cohorts, trios=trios, pedigrees=pedigrees)
    karyomapping_analyses = KaryomappingAnalysis.filter_for_user(user).filter(trio__in=trios)
    analyses_list = [i[0].pk for i in analysis_details]
    analyses_with_tags = Analysis.objects.filter(pk__in=analyses_list, varianttag__isnull=False)
    variant_tag_genome_build_names = analyses_with_tags.order_by("genome_build").values_list("genome_build", flat=True).distinct()

    context.update({"analysis_details": analysis_details,
                    "karyomapping_analyses": karyomapping_analyses,
                    "show_output_nodes": settings.ANALYSIS_RELATED_DOWNLOAD_OUTPUT_NODES,
                    "show_sample_info": show_sample_info,
                    "variant_tag_genome_build_names": variant_tag_genome_build_names,
                    "analysis_ids_list": analyses_list})


@register.inclusion_tag("analysis/tags/related_analyses_for_samples.html", takes_context=True)
def related_analyses_for_samples(context, samples, show_sample_info):
    cohorts = Cohort.objects.filter(cohortsample__sample__in=samples).distinct()
    trio_sample_q_list = [Q(**{f"{trio_field}__sample__in": samples}) for trio_field in ["mother", "father", "proband"]]
    trios = Trio.objects.filter(cohort__in=cohorts).filter(reduce(operator.or_, trio_sample_q_list)).distinct()
    pedigrees = Pedigree.objects.filter(cohort__in=cohorts).filter(cohortsamplepedfilerecord__cohort_sample__sample__in=samples).distinct()
    sample_mutational_signatures = MutationalSignature.objects.filter(sample__in=samples).distinct().select_related("sample")

    update_context_with_related_analysis(context, samples, cohorts, trios, pedigrees, show_sample_info=show_sample_info)
    context["samples"] = samples
    context["sample_mutational_signatures"] = sample_mutational_signatures
    return context


@register.inclusion_tag("analysis/tags/related_analyses_for_cohort.html", takes_context=True)
def related_analyses_for_cohort(context, cohort):
    pedigrees = cohort.pedigree_set.all()
    trios = cohort.trio_set.all()
    cohorts = [cohort] + list(cohort.sub_cohort_set.all())

    update_context_with_related_analysis(context, cohort.get_samples(), cohorts, trios, pedigrees)
    context["cohort"] = cohort
    return context


@register.inclusion_tag("analysis/tags/related_analyses_for_trio.html", takes_context=True)
def related_analyses_for_trio(context, trio):
    pedigrees = trio.cohort.pedigree_set.all()
    update_context_with_related_analysis(context, trio.get_samples(), [trio.cohort], [trio], pedigrees)
    context["trio"] = trio
    return context


@register.inclusion_tag("analysis/tags/related_analyses_for_pedigree.html", takes_context=True)
def related_analyses_for_pedigree(context, pedigree):
    trios = pedigree.cohort.trio_set.all()
    update_context_with_related_analysis(context, pedigree.get_samples(), [pedigree.cohort], trios, [pedigree])
    context["pedigree"] = pedigree
    return context


@register.inclusion_tag("analysis/tags/analysis_templates_tag.html", takes_context=True)
def analysis_templates_tag(context, genome_build, autocomplete_field=True, has_somatic_sample=False, has_sample_gene_list=False, requires_sample_gene_list=None,
                           **kwargs):
    user = context["user"]
    single_model_args = {"sample", "cohort", "trio", "pedigree"}
    params_error_message = f"analysis_templates_tag should be passed dict with exactly one Model value for {','.join(single_model_args)}. Args: {kwargs}"

    hidden_inputs = {}  # values should be primary keys
    klass = None
    for k, v in kwargs.items():
        if isinstance(v, Model):
            if k in single_model_args:
                if klass:
                    raise ValueError(params_error_message)
                klass = type(v)
            v = v.pk

        hidden_inputs[k] = v
    if klass is None:
        raise ValueError(params_error_message)

    # Hack to add proband as sample
    if trio := kwargs.get("trio"):
        hidden_inputs["sample"] = trio.proband.sample_id

    class_name = klass._meta.label
    # Show/Hide AnalysisTemplateVersions based on requires_sample_gene_list

    missing_templates = []
    requires_sample_somatic = None
    if not has_somatic_sample:
        requires_sample_somatic = False
        missing_templates.append("somatic sample")

    if not has_sample_gene_list:
        requires_sample_gene_list = False
        missing_templates.append("sample gene list")

    AnalysisTemplateForm = get_analysis_template_form_for_variables_only_of_class(class_name,
                                                                                  autocomplete_field=autocomplete_field,
                                                                                  requires_sample_somatic=requires_sample_somatic,
                                                                                  requires_sample_gene_list=requires_sample_gene_list)

    analysis_template_links = AnalysisTemplate.filter(user, class_name=class_name,
                                                      requires_sample_somatic=requires_sample_somatic,
                                                      requires_sample_gene_list=requires_sample_gene_list,
                                                      atv_kwargs={"appears_in_links": True})

    flattened_uuid = str(uuid.uuid4()).replace("-", "_")
    return {
        "genome_build": genome_build,
        "flattened_uuid": flattened_uuid,
        "autocomplete_field": autocomplete_field,
        "analysis_template_form": AnalysisTemplateForm(prefix=flattened_uuid),
        "analysis_template_links": analysis_template_links,
        "hidden_inputs": hidden_inputs,
        "missing_templates": ", ".join(missing_templates),
    }


@register.inclusion_tag("analysis/tags/analysis_output_node_downloads.html", takes_context=True)
def analysis_output_node_downloads(context, analysis):
    or_filters = [
        Q(output_node=True)
    ]
    # Always add tag node
    # q_tag_node = Q(name=TagNode.ANALYSIS_TAGS_NAME, visible=False)
    # or_filters.append(q_tag_node)
    q = reduce(operator.or_, or_filters)
    qs_output_nodes = analysis.analysisnode_set.all().filter(q).order_by("pk").select_subclasses()
    # output_node_values = list(qs_output_nodes.values_list("name", "status", "count"))

    return {
        "analysis": analysis,
        "output_nodes": qs_output_nodes,
    }
