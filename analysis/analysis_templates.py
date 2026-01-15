import logging
import re

from django.conf import settings
from django.contrib.auth.models import User

from analysis.models import Analysis, AnalysisNode, AnalysisTemplate, AnalysisTemplateRun, \
    AnalysisTemplateRunArgument, SampleAnalysisTemplateRun
from analysis.models.nodes.node_utils import get_toposorted_nodes, reload_analysis_nodes
from analysis.related_analyses import get_related_analysis_details_for_samples
from genes.models import ActiveSampleGeneList
from library.guardian_utils import add_public_group_read_permission
from snpdb.models import Sample, GenomeBuild


def run_analysis_template(analysis_template: AnalysisTemplate,
                          genome_build: GenomeBuild,
                          user: User = None,
                          **kwargs) -> AnalysisTemplateRun:

    template_run = AnalysisTemplateRun.create(analysis_template, genome_build, user=user)
    template_run.populate_arguments(kwargs)
    populate_analysis_from_template_run(template_run)
    return template_run


def populate_analysis_from_template_run(template_run):
    """ Relies on populate_arguments being called to create AnalysisTemplateRunArguments for AnalysisVariables """

    template_run.populate_analysis_name()
    node_field_values = template_run.get_node_field_values()

    nodes_qs = AnalysisNode.objects.filter(pk__in=node_field_values).select_subclasses()
    nodes_with_expected_errors = []
    for group in get_toposorted_nodes(nodes_qs):
        for node in group:
            # Set the fields (want to do all at once before save)
            for field, value in node_field_values[node.pk].items():
                setattr(node, field, value)
                node.queryset_dirty = True

            try:
                node.save()
                error_message = node.get_errors(include_parent_errors=False, flat=True)
            except Exception as e:
                error_message = str(e)

            if error_message:
                AnalysisTemplateRunArgument.objects.filter(variable__node=node).update(error=error_message)
                if node.hide_node_and_descendants_upon_template_configuration_error:
                    nodes_with_expected_errors.append(node)

    if nodes_with_expected_errors:
        descendants_node_ids_to_hide = set()
        for node in nodes_with_expected_errors:
            descendants_node_ids_to_hide.update([n.pk for n in node.descendants_set()])

        template_args = AnalysisTemplateRunArgument.objects.filter(variable__node_id__in=descendants_node_ids_to_hide)
        template_args.update(error="Hidden due to ancestor error")
        nodes_to_hide = {n.pk for n in nodes_with_expected_errors} | descendants_node_ids_to_hide
        AnalysisNode.objects.filter(pk__in=nodes_to_hide).update(visible=False)

    print("Everything fine - recalculating everything!")
    reload_analysis_nodes(template_run.analysis.pk)


def get_sample_analysis(sample: Sample, analysis_template: AnalysisTemplate) -> Analysis:
    """ Can only be 1 per sample/template run combo """

    try:
        satr = SampleAnalysisTemplateRun.objects.get(sample=sample,
                                                     analysis_template_run__template_version__template=analysis_template,
                                                     analysis_template_run__template_version__version=analysis_template.latest_version())
        return satr.analysis_template_run.analysis
    except SampleAnalysisTemplateRun.DoesNotExist:
        at_run = run_analysis_template(analysis_template, sample.genome_build, sample=sample)
        at_run.analysis.visible = False
        at_run.analysis.save()
        add_public_group_read_permission(at_run.analysis)
        SampleAnalysisTemplateRun.objects.create(sample=sample, analysis_template_run=at_run)
        return at_run.analysis


def _get_auto_launch_analysis_templates_for_sample(user, sample, skip_already_analysed=False):
    if skip_already_analysed:
        if get_related_analysis_details_for_samples(user, [sample]):
            return []

    sample_enrichment_kit_name = None
    if sample.enrichment_kit:
        sample_enrichment_kit_name = sample.enrichment_kit.name

    analysis_templates = []
    for enrichment_kit_name, sample_regex, analysis_template_name in settings.ANALYSIS_AUTO_LAUNCH_SAMPLE_ANALYSIS_TEMPLATE_CONFIG:
        if sample_enrichment_kit_name == enrichment_kit_name:
            if sample_regex and not re.match(sample_regex, sample.name):
                continue
            analysis_templates.append(AnalysisTemplate.objects.get(name=analysis_template_name))
    return analysis_templates

def auto_launch_analysis_templates_for_sample(user, sample, analysis_description=None, skip_already_analysed=False):
    for analysis_template in _get_auto_launch_analysis_templates_for_sample(user, sample,
                                                                            skip_already_analysed=skip_already_analysed):
        template_version = analysis_template.active
        template_arguments = {"sample": sample}
        if template_version.requires_sample_gene_list:
            try:
                template_arguments["sample_gene_list"] = sample.activesamplegenelist
            except ActiveSampleGeneList.DoesNotExist:
                logging.warning("Skipping auto analysis '%s' for sample: %s. Will try again if QC Gene Lists created", analysis_template, sample)
                continue
        template_run = AnalysisTemplateRun.create(analysis_template, sample.genome_build, user=user)
        if analysis_description:
            template_run.analysis.description = analysis_description
            template_run.analysis.save()
        template_run.populate_arguments(template_arguments)
        populate_analysis_from_template_run(template_run)
