import logging
import re
from dataclasses import dataclass
from typing import Optional, List

from auditlog.context import disable_auditlog
from django.conf import settings
from django.contrib.auth.models import User

from analysis.models import Analysis, AnalysisNode, AnalysisTemplate, AnalysisTemplateRun, \
    AnalysisTemplateRunArgument, SampleAnalysisTemplateRun, CohortAnalysisTemplateRun, AutoLaunchAnalysisTemplate
from analysis.models.nodes.node_utils import get_toposorted_nodes, reload_analysis_nodes
from analysis.related_analyses import get_related_analysis_details_for_samples
from genes.models import ActiveSampleGeneList
from library.guardian_utils import add_public_group_read_permission
from snpdb.models import Sample, GenomeBuild, Cohort


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

    with disable_auditlog():
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


def _get_single_template_run_analysis(klass, analysis_template: AnalysisTemplate,
                                      obj, field_name) -> Analysis:
    """ Can only be 1 per sample/template run combo """
    kwargs = {
        field_name: obj,
    }

    try:
        satr = klass.objects.get(analysis_template_run__template_version__template=analysis_template,
                                 analysis_template_run__template_version__version=analysis_template.latest_version(),
                                 **kwargs)
        return satr.analysis_template_run.analysis
    except klass.DoesNotExist:
        at_run = run_analysis_template(analysis_template, obj.genome_build, **kwargs)
        at_run.analysis.visible = False
        at_run.analysis.save()
        add_public_group_read_permission(at_run.analysis)
        klass.objects.create(analysis_template_run=at_run, **kwargs)
        return at_run.analysis


def get_sample_analysis(sample: Sample, analysis_template: AnalysisTemplate) -> Analysis:
    """ Can only be 1 per sample/template run combo """

    return _get_single_template_run_analysis(SampleAnalysisTemplateRun, analysis_template,
                                             sample, "sample")


def get_cohort_analysis(cohort: Cohort, analysis_template: AnalysisTemplate) -> Analysis:
    return _get_single_template_run_analysis(CohortAnalysisTemplateRun, analysis_template,
                                             cohort, "cohort")


def _get_auto_launch_analysis_templates_for_sample(user, sample, skip_already_analysed=False):
    if skip_already_analysed:
        if get_related_analysis_details_for_samples(user, [sample]):
            return []

    sample_enrichment_kit_name = None
    if sample.enrichment_kit:
        sample_enrichment_kit_name = sample.enrichment_kit.name

    matches = get_auto_launch_analysis_template_matches(user, sample_enrichment_kit_name, sample.name)
    return [m.analysis_template for m in matches if m.match]

@dataclass
class AutoLaunchAnalysisTemplateMatch:
    enrichment_kit_name: Optional[str]
    enrichment_kit_match: bool
    sample_regex: Optional[str]
    sample_regex_match: bool
    analysis_template: AnalysisTemplate

    @property
    def match(self) -> bool:
        return self.enrichment_kit_match and self.sample_regex_match


def get_auto_launch_analysis_template_matches(user, sample_enrichment_kit_name, sample_name) -> List[AutoLaunchAnalysisTemplateMatch]:
    matches = []
    templates_qs = AnalysisTemplate.filter_for_user(user)
    for auto_launch in AutoLaunchAnalysisTemplate.objects.filter(template__in=templates_qs):
        enrichment_kit_name = None
        if enrichment_kit := auto_launch.enrichment_kit:
            enrichment_kit_name = enrichment_kit.name

        sample_regex_match = True  # blank means anything
        if auto_launch.sample_regex:
            sample_regex_match = bool(re.match(auto_launch.sample_regex, sample_name))
        match = AutoLaunchAnalysisTemplateMatch(enrichment_kit_name=enrichment_kit_name,
                                                enrichment_kit_match=sample_enrichment_kit_name == enrichment_kit_name,
                                                sample_regex=auto_launch.sample_regex,
                                                sample_regex_match=sample_regex_match,
                                                analysis_template=auto_launch.template)
        matches.append(match)
    return matches

def auto_launch_analysis_templates_for_sample(user, sample, analysis_description=None, skip_already_analysed=False):
    for analysis_template in _get_auto_launch_analysis_templates_for_sample(user, sample,
                                                                            skip_already_analysed=skip_already_analysed):
        template_version = analysis_template.active
        template_arguments = {"sample": sample}
        if template_version.requires_sample_gene_list:
            try:
                template_arguments["sample_gene_list"] = sample.activesamplegenelist.sample_gene_list
            except ActiveSampleGeneList.DoesNotExist:
                logging.warning("Skipping auto analysis '%s' for sample: %s. Will try again if QC Gene Lists created", analysis_template, sample)
                continue
        template_run = AnalysisTemplateRun.create(analysis_template, sample.genome_build, user=user)
        if analysis_description:
            template_run.analysis.description = analysis_description
            template_run.analysis.save()
        template_run.populate_arguments(template_arguments)
        populate_analysis_from_template_run(template_run)
