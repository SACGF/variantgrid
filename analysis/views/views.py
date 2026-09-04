from collections import Counter
from functools import cached_property

import numpy as np
from django.apps import apps
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.core.exceptions import PermissionDenied
from django.http.response import JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls.base import reverse
from django.utils.text import slugify
from django.views.decorators.http import require_POST

from analysis import forms
from analysis.analysis_import_export import analysis_export_to_dict
from analysis.analysis_templates import (
    get_auto_launch_analysis_template_matches,
    populate_analysis_from_template_run,
)
from analysis.forms import (
    AnalysisChoiceForm,
    AnalysisTemplateAutoLaunchForm,
    AnalysisTemplateForm,
    AnalysisTemplateTypeChoiceForm,
    AnalysisTemplateVersionForm,
    AutoLaunchFormSet,
    SelectGridColumnForm,
)
from analysis.models import (
    Analysis,
    AnalysisNode,
    AnalysisTemplate,
    AnalysisTemplateRun,
    AnalysisTemplateVersion,
    AnalysisVariable,
    TagNode,
    VariantTag,
)
from analysis.models.enums import MinimisationResultType, NodeStatus
from analysis.models.mutational_signatures import MutationalSignature
from analysis.models.nodes.analysis_node import AnalysisClassification
from analysis.models.nodes.node_counts import get_node_count_colors, get_tag_node_count_colors
from analysis.models.nodes.node_types import get_node_display_data_by_class_name, get_node_types_hash
from analysis.variant_tag_operations import retire_requires_classification_tags
from analysis.views.analysis_permissions import get_analysis_or_404
from annotation.models.models import MutationalSignatureInfo
from classification.views.views import (
    CreateClassificationForVariantView,
    create_classification_object,
)
from library.django_utils import add_save_message, get_field_counts, set_form_read_only
from library.guardian_utils import is_superuser
from seqauto.models import EnrichmentKit
from snpdb.forms import SampleChoiceForm
from snpdb.models import (
    GenomeBuild,
    JobsControl,
    Sample,
    UserSettings,
    get_igv_data,
)
from variantgrid.tips import get_tips


def analysis_list(request):
    form = forms.CreateAnalysisForm(request.POST or None, user=request.user)
    if request.method == "POST":
        if form.is_valid():
            analysis = form.save()
            return redirect(analysis)
        add_save_message(request, False, "Analysis")

    context = {"create_analysis_form": form,
               "analysis_choice_form": AnalysisChoiceForm()}
    return render(request, 'analysis/analyses.html', context)


def analysis_templates(request):
    form = forms.CreateAnalysisTemplateForm(request.POST or None, user=request.user)
    if request.method == "POST":
        if form.is_valid():
            analysis_template = form.save()
            return redirect(analysis_template.analysis)

    context = {"create_analysis_template_form": form,
               "analysis_template_choice_form": AnalysisTemplateTypeChoiceForm()}
    return render(request, 'analysis/analysis_templates.html', context)


def analysis_templates_auto_launch(request):
    """ Shows how the auto launch will work to users """

    template_auto_launch_form = AnalysisTemplateAutoLaunchForm(request.POST or None)

    sample_enrichment_kit_name = None
    sample_name = ""
    if request.method == "POST":
        if sample_enrichment_kit_id := request.POST.get("enrichment_kit"):
            if enrichment_kit := EnrichmentKit.objects.filter(pk=sample_enrichment_kit_id).first():
                sample_enrichment_kit_name = enrichment_kit.name
        sample_name = request.POST["example_sample_name"]

    auto_launch_analysis_template_matches = get_auto_launch_analysis_template_matches(request.user, sample_enrichment_kit_name, sample_name)
    context = {
        "auto_launch_analysis_template_matches": auto_launch_analysis_template_matches,
        "template_auto_launch_form": template_auto_launch_form,
        "sample_enrichment_kit_name": sample_enrichment_kit_name,
        "sample_name": sample_name,
    }
    return render(request, 'analysis/analysis_templates_auto_launch.html', context)


def get_analysis_settings(user, analysis):
    user_settings = UserSettings.get_for_user(user)
    igv_data = get_igv_data(user, genome_build=analysis.genome_build)
    if analysis.canonical_transcript_collection:
        canonical_transcript_collection = str(analysis.canonical_transcript_collection)
    else:
        canonical_transcript_collection = ""

    analysis_settings = {
        "annotation_version": analysis.annotation_version_id,
        "node_count_types": analysis.get_node_count_types(),
        "canonical_transcript_collection": canonical_transcript_collection,
        "grid_sample_label_template": analysis.grid_sample_label_template,
        "variant_tag_stale_days": analysis.variant_tag_stale_days,
        "show_igv_links": analysis.show_igv_links,
        "igv_data": igv_data,
        "open_variant_details_in_new_window": user_settings.variant_link_in_analysis_opens_new_tab,
        "genome_build": str(analysis.genome_build),
    }
    return analysis_settings


def view_analysis(request, analysis_id, active_node_id=0):
    analysis = get_analysis_or_404(request.user, analysis_id)

    nodes = analysis.analysisnode_set.filter(visible=True).select_subclasses()
    node_classes_kwargs = {}
    if analysis.lock_input_sources:
        node_classes_kwargs = {"source_nodes": False}

    node_help_dict = {}
    user_settings = UserSettings.get_for_user(request.user)
    if user_settings.tool_tips:
        node_help_dict = {label: subclass.get_help_text() for label, subclass in get_node_types_hash().items()}

    analysis_variables = [[av.node_id, av.field] for av in AnalysisVariable.objects.filter(node__analysis=analysis)]
    analysis_tags_node = TagNode.get_analysis_tags_node(analysis)

    context = {
        "node_classes_form": forms.AnalysisNodeClassesForm(**node_classes_kwargs),
        "nodes": nodes,
        "node_count_colors": get_node_count_colors("color") + get_tag_node_count_colors(request.user, "color"),
        "analysis": analysis,
        "analysis_settings": get_analysis_settings(request.user, analysis),
        "analysis_tags_node": analysis_tags_node,
        "active_node_id": active_node_id,
        "node_help": node_help_dict,
        "node_types_display": get_node_display_data_by_class_name(),
        "analysis_variables": analysis_variables,
        "has_write_permission": analysis.can_write(request.user),
        "warnings": analysis.get_toolbar_warnings(request.user),
        "loading_animations": user_settings.grid_loading_animations,
        "tips": [t.tip for t in get_tips(app="analysis")] if user_settings.show_tips else [],
    }
    return render(request, 'analysis/analysis.html', context)


def view_active_node(analysis, active_node=None):
    kwargs = {"analysis_id": analysis.pk}
    url = reverse('analysis', kwargs=kwargs)
    if active_node:
        kwargs["active_node_id"] = active_node.pk
        url = reverse('analysis_node', kwargs=kwargs)
    return redirect(url)


@require_POST
def create_analysis_from_template(request, genome_build_name):
    data = request.POST.dict()
    tag_uuid = data.pop("tag_uuid")
    analysis_template_version_key = f"{tag_uuid}-analysis_template_version"
    analysis_template_version_id = data.pop(analysis_template_version_key)
    template_version = get_object_or_404(AnalysisTemplateVersion, pk=analysis_template_version_id)
    # A draft is only for the people who can change the template - the active version is for everyone
    template_version.template.check_permission(request.user, write=not template_version.active)

    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    template_run = AnalysisTemplateRun.create(template_version.template, genome_build, user=request.user,
                                              template_version=template_version)
    template_run.populate_arguments(data)

    # Block creation against archived sources. populate_arguments has already resolved
    # each variable to its model instance; just check the common `data_archived` interface.
    for arg in template_run.analysistemplaterunargument_set.all():
        if not arg.object_pk:
            continue
        klass = apps.get_model(arg.variable.class_name)
        obj = klass.objects.filter(pk=arg.object_pk).first()
        if obj is not None and getattr(obj, "data_archived", False):
            raise PermissionDenied(
                "Cannot create new analysis: source data is archived. "
                "Restore the data first."
            )

    populate_analysis_from_template_run(template_run)

    return view_active_node(template_run.analysis, None)


def analysis_editor_and_grid(request, analysis_id, stand_alone=False):
    analysis = get_analysis_or_404(request.user, analysis_id)
    context = {"select_grid_column_form": SelectGridColumnForm(),
               "stand_alone": stand_alone,
               "analysis": analysis}
    return render(request, 'analysis/analysis_editor_and_grid.html', context)


def stand_alone_analysis_editor_and_grid(request, analysis_id):
    return analysis_editor_and_grid(request, analysis_id, True)


def _get_form(request, formcls, prefix, **kwargs):
    # From https://stackoverflow.com/a/33028994/295724
    data = request.POST if prefix in request.POST else None
    return formcls(data, prefix=prefix, **kwargs)


def analysis_template_settings(request, pk):
    analysis_template = AnalysisTemplate.get_for_user(request.user, pk)
    has_write_permission = analysis_template.can_write(request.user)

    at_form = _get_form(request, AnalysisTemplateForm, 'at-pre', instance=analysis_template)
    formset_data = request.POST if 'auto-launch' in request.POST else None
    formset = AutoLaunchFormSet(formset_data, prefix='auto-launch', instance=analysis_template)

    # The draft's settings need to be editable before it goes live, so it gets a form as well
    editable_versions = [atv for atv in (analysis_template.active, analysis_template.draft) if atv]
    atv_forms = {atv.pk: _get_form(request, AnalysisTemplateVersionForm, f"atv-{atv.pk}", instance=atv)
                 for atv in editable_versions}

    if not has_write_permission:
        set_form_read_only(at_form)
        for atv_form in atv_forms.values():
            set_form_read_only(atv_form)
        set_form_read_only(formset)
        messages.add_message(request, messages.WARNING, "You can view but not modify this data.")

    if request.method == 'POST':
        if not has_write_permission:
            raise PermissionDenied(f"Don't have permission to modify {analysis_template}")

        if at_form and at_form.is_bound:
            valid = at_form.is_valid()
            if valid:
                at_form.save()
            add_save_message(request, valid, "Analysis Template")

        if formset.is_bound:
            valid = formset.is_valid()
            if valid:
                formset.save()
                # Rebuild unbound, so saved rows are re-displayed along with a new blank one to add another
                formset = AutoLaunchFormSet(prefix='auto-launch', instance=analysis_template)
            add_save_message(request, valid, "Auto Launch Config")

        for atv in editable_versions:
            atv_form = atv_forms[atv.pk]
            if atv_form.is_bound:
                valid = atv_form.is_valid()
                if valid:
                    atv_form.save()
                add_save_message(request, valid, f"Template Version v.{atv.version}")

    versions = analysis_template.analysistemplateversion_set.order_by("-version")
    context = {
        "analysis_template": analysis_template,
        "analysis_template_versions": [(atv, atv_forms.get(atv.pk)) for atv in versions],
        "has_versions": versions.exists(),
        "at_form": at_form,
        "at_formset": formset,
        "has_write_permission": has_write_permission,
    }
    return render(request, 'analysis/analysis_template_settings.html', context)


@require_POST
def analysis_template_version_activate(request, pk):
    atv = get_object_or_404(AnalysisTemplateVersion, pk=pk)
    atv.template.check_can_write(request.user)
    atv.activate()
    messages.add_message(request, messages.SUCCESS, f"v.{atv.version} is now the active version")
    return redirect("analysis_template_settings", pk=atv.template_id)


def analysis_export(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    filename = f"analysis_{analysis.pk}_{slugify(analysis.name)[:40]}.json"
    response = JsonResponse(analysis_export_to_dict(analysis))
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


@user_passes_test(is_superuser)
def view_analysis_issues(request):
    if request.method == "POST":
        if "unpause-jobs" in request.POST:
            JobsControl.resume(by=str(request.user))
            messages.add_message(request, messages.INFO,
                                 "Resumed analysis + annotation job dispatch")
        if "pause-jobs" in request.POST:
            JobsControl.pause(reason=f"Paused from analysis issues page by {request.user}",
                              by=str(request.user))
            messages.add_message(request, messages.WARNING,
                                 "Paused analysis + annotation job dispatch")

    all_nodes = AnalysisNode.objects.all()
    field_counts = get_field_counts(all_nodes, "status")
    summary_data = Counter()
    for status, count in field_counts.items():
        summary = NodeStatus.get_summary_state(status)
        summary_data[summary] += count

    field_counts = {NodeStatus(k).label: v for k, v in field_counts.items()}
    # Don't force-create the singleton just by viewing the page
    jobs_control = JobsControl.objects.filter(pk=JobsControl.SINGLETON_PK).first()
    context = {"nodes_status_summary": summary_data,
               "field_counts": field_counts,
               "jobs_control": jobs_control,
               "jobs_paused": bool(jobs_control and jobs_control.paused)}
    return render(request, 'analysis/view_analysis_issues.html', context)


def view_mutational_signature(request, pk):
    mutational_signature = get_object_or_404(MutationalSignature, pk=pk)
    # Ensure you can see this sample - sample
    mutational_signature.sample.check_can_view(request.user)
    results_qs = mutational_signature.mutationalsignatureminimisationresult_set.all()
    bootstrapped_results_qs = results_qs.filter(result_type=MinimisationResultType.BOOTSTRAPPED)

    sorted_index = np.flip(np.argsort(mutational_signature.mean), axis=0)
    bootstrapped_data = []
    for r in bootstrapped_results_qs.order_by("iteration"):
        bootstrapped_data.append(r.solution_array)
    sorted_bootstrap_transposed = np.array(bootstrapped_data).T[sorted_index]
    sorted_bootstrap_std = sorted_bootstrap_transposed.std(axis=1)
    sorted_mean = np.array(mutational_signature.mean)[sorted_index]
    labels = [f"Sig {i}" for i in sorted_index + 1]  # signatures are 1 based
    signature_info = list(MutationalSignatureInfo.objects.all().order_by("signature_id"))
    sorted_signature_info = [signature_info[i] for i in sorted_index]
    sorted_data = zip(labels, sorted_mean, sorted_bootstrap_std, sorted_signature_info)

    # Data Fit plot
    mut_count_qs = mutational_signature.mutationalsignaturemutationcount_set.all().order_by("pk")
    mutation_types = []
    mutation_type_counts = []

    for mutation_type, count in mut_count_qs.values_list("mutation_type", "count"):
        mutation_types.append(mutation_type)
        mutation_type_counts.append(count)

    context = {"mutational_signature": mutational_signature,
               "labels": labels,
               "sorted_bootstrap_transposed": sorted_bootstrap_transposed.tolist(),
               "mutation_types": mutation_types,
               "mutation_type_counts": mutation_type_counts,
               "sorted_data": sorted_data}
    return render(request, 'analysis/view_mutational_signature.html', context)


class CreateClassificationForVariantTagView(CreateClassificationForVariantView):
    template_name = "analysis/create_classification_for_variant_tag.html"

    def _get_variant(self):
        return self.variant_tag.variant

    def _get_genome_build(self) -> GenomeBuild:
        return self.variant_tag.genome_build

    def _get_form_post_url(self) -> str:
        if self.variant_tag.analysis:
            return reverse("create_classification_for_analysis", kwargs={"analysis_id": self.variant_tag.analysis_id})
        else:
            # Just for variant
            return super()._get_form_post_url()

    def _get_sample_form(self):
        # If we have a node with input samples, use that. Then fall back on all samples in analysis.
        # Otherwise fall back on default (all samples in DB visible to user)
        samples = None
        if self.variant_tag.analysis:
            if self.variant_tag.node:
                samples = self.variant_tag.node.get_subclass().get_samples()

            if not samples:
                samples = self.variant_tag.analysis.get_samples()

        if samples:
            form = SampleChoiceForm()
            form.fields['sample'].required = False
            form.fields['sample'].queryset = Sample.objects.filter(pk__in=[s.pk for s in samples])
        else:
            form = super()._get_sample_form()
        return form

    @cached_property
    def variant_tag(self):
        variant_tag_id = self.kwargs["variant_tag_id"]
        variant_tag = VariantTag.get_for_user(self.request.user, pk=variant_tag_id)
        return variant_tag

    def get_context_data(self, *args, **kwargs):
        try:
            context = super().get_context_data(*args, **kwargs)
            context["variant_tag"] = self.variant_tag

            if not self.variant_tag.can_write(self.request.user):
                if self.variant_tag.analysis:
                    read_only_message = "You have read-only access to this analysis. You can create a " \
                                        "classification but it will not be linked to the analysis and the " \
                                        f"{settings.TAG_REQUIRES_CLASSIFICATION} tag will not be deleted."
                else:
                    read_only_message = "You have read-only access to this tag. You can create a classification " \
                                        f"but the {settings.TAG_REQUIRES_CLASSIFICATION} tag will not be deleted."
                messages.add_message(self.request, messages.WARNING, read_only_message)
        except VariantTag.DoesNotExist:
            variant_tag_id = self.kwargs["variant_tag_id"]
            msg = f"The VariantTag ({variant_tag_id}) does not exist. It may have been deleted or already classified."
            context = {"error_message": msg}
        return context


@require_POST
def create_classification_for_analysis(request, analysis_id):
    classification = create_classification_object(request)
    analysis = Analysis.get_for_user(request.user, pk=analysis_id)

    if analysis.can_write(request.user):
        AnalysisClassification.objects.create(analysis=analysis, classification=classification)
        retire_requires_classification_tags(classification, analysis, request.user)
    return redirect(classification.get_edit_url())
