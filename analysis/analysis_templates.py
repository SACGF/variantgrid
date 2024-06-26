from auditlog.context import disable_auditlog
from django.contrib.auth.models import User

from analysis.models import Analysis, AnalysisNode, AnalysisTemplate, AnalysisTemplateRun, \
    AnalysisTemplateRunArgument, SampleAnalysisTemplateRun, CohortAnalysisTemplateRun
from analysis.models.nodes.node_utils import get_toposorted_nodes, reload_analysis_nodes
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
