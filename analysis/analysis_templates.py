from typing import Optional

from auditlog.context import disable_auditlog
from django.contrib.auth.models import User

from analysis.enums import AnalysisTemplateRunTrigger
from analysis.models import (
    Analysis,
    AnalysisNode,
    AnalysisTemplate,
    AnalysisTemplateRun,
    AnalysisTemplateRunArgument,
    AnalysisTemplateVersion,
    CohortAnalysisTemplateRun,
    SampleAnalysisTemplateRun,
)
from analysis.models.nodes.analysis_node import AnalysisDag
from analysis.models.nodes.node_utils import get_toposorted_nodes, reload_analysis_nodes
from annotation.models import AnnotationVersion
from library.django_utils.database_utils import advisory_lock
from library.guardian_utils import add_public_group_read_permission
from snpdb.models import Cohort, GenomeBuild, Sample


def run_analysis_template(analysis_template: AnalysisTemplate,
                          genome_build: GenomeBuild,
                          trigger: AnalysisTemplateRunTrigger,
                          user: User = None,
                          template_version: AnalysisTemplateVersion = None,
                          annotation_version: AnnotationVersion = None,
                          **kwargs) -> AnalysisTemplateRun:
    template_run = AnalysisTemplateRun.create(analysis_template, genome_build, trigger, user=user,
                                              template_version=template_version,
                                              annotation_version=annotation_version)
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
            dag = AnalysisDag(template_run.analysis.pk)
            descendants_node_ids_to_hide = set()
            for node in nodes_with_expected_errors:
                descendants_node_ids_to_hide.update(dag.descendant_ids(node.pk))

            template_args = AnalysisTemplateRunArgument.objects.filter(variable__node_id__in=descendants_node_ids_to_hide)
            template_args.update(error="Hidden due to ancestor error")
            nodes_to_hide = {n.pk for n in nodes_with_expected_errors} | descendants_node_ids_to_hide
            AnalysisNode.objects.filter(pk__in=nodes_to_hide).update(visible=False)

        reload_analysis_nodes(template_run.analysis.pk)


def _get_single_template_run_analysis(klass, analysis_template: AnalysisTemplate,
                                      obj, field_name, trigger: AnalysisTemplateRunTrigger,
                                      annotation_version: Optional[AnnotationVersion] = None) -> Analysis:
    """ One per sample (or cohort), template version and annotation version - a run left on a superseded
        annotation version is passed over for a new one, so the sample tab and the annotated exports show
        the current annotation (#2209). annotation_version defaults to the build's latest """
    kwargs = {
        field_name: obj,
    }

    template_version = analysis_template.active
    if template_version is None:
        raise ValueError(f"{analysis_template} has no active version")
    if annotation_version is None:
        annotation_version = AnnotationVersion.latest(obj.genome_build)

    # Concurrent requests (a page and an export) would each create a run - the second waits and finds the first's
    with advisory_lock(klass._meta.db_table, obj.pk):
        satr = klass.objects.filter(analysis_template_run__template_version=template_version,
                                    analysis_template_run__analysis__annotation_version=annotation_version,
                                    **kwargs).order_by("-pk").first()
        if satr:
            return satr.analysis_template_run.analysis

        at_run = run_analysis_template(analysis_template, obj.genome_build, trigger, template_version=template_version,
                                       annotation_version=annotation_version, **kwargs)
        at_run.analysis.visible = False
        at_run.analysis.save()
        add_public_group_read_permission(at_run.analysis)
        klass.objects.create(analysis_template_run=at_run, **kwargs)
        return at_run.analysis


def get_sample_analysis(sample: Sample, analysis_template: AnalysisTemplate,
                        annotation_version: Optional[AnnotationVersion] = None) -> Analysis:
    return _get_single_template_run_analysis(SampleAnalysisTemplateRun, analysis_template,
                                             sample, "sample", AnalysisTemplateRunTrigger.SAMPLE,
                                             annotation_version=annotation_version)


def get_cohort_analysis(cohort: Cohort, analysis_template: AnalysisTemplate,
                        annotation_version: Optional[AnnotationVersion] = None) -> Analysis:
    return _get_single_template_run_analysis(CohortAnalysisTemplateRun, analysis_template,
                                             cohort, "cohort", AnalysisTemplateRunTrigger.COHORT,
                                             annotation_version=annotation_version)
