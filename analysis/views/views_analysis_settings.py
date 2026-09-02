import json
from collections import defaultdict

from auditlog.models import LogEntry
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.core.exceptions import PermissionDenied
from django.http.response import JsonResponse
from django.shortcuts import redirect, render
from django.utils import timezone
from django.views.decorators.http import require_POST

from analysis import forms
from analysis.models import AnalysisLock, AnalysisTemplateRunArgument
from analysis.models.enums import AnalysisTemplateType
from analysis.models.nodes import node_utils
from analysis.models.nodes.analysis_node import NodeCount
from analysis.models.nodes.node_counts import get_node_counts_mine_and_available
from analysis.views.analysis_permissions import get_analysis_or_404
from analysis.views.views import get_analysis_settings
from library.django_utils import add_save_message, set_form_read_only
from library.guardian_utils import is_superuser
from library.utils import defaultdict_to_dict


def view_analysis_settings(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    analysis_settings = get_analysis_settings(request.user, analysis)

    context = {"analysis": analysis,
               "new_analysis_settings": analysis_settings,
               "has_write_permission": analysis.can_write(request.user),
               "can_unlock": analysis.can_unlock(request.user)}
    return render(request, 'analysis/analysis_settings.html', context)


def _get_hidden_nodes(analysis) -> list[dict]:
    """ Nodes a template run hid due to configuration errors, with the reason it stored against them """
    hidden_qs = analysis.analysisnode_set.filter(visible=False).order_by("pk")
    errors_by_node_id = defaultdict(list)
    args_qs = AnalysisTemplateRunArgument.objects.filter(variable__node__in=hidden_qs, error__isnull=False)
    for node_id, error in args_qs.values_list("variable__node_id", "error"):
        errors_by_node_id[node_id].append(error)
    return [{"name": node.name, "errors": errors_by_node_id[node.pk]} for node in hidden_qs]


def analysis_settings_details_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    old_annotation_version = analysis.annotation_version
    old_horizontal_mode = analysis.analysis_horizontal_mode
    form = forms.AnalysisForm(request.POST or None, user=request.user, instance=analysis)
    has_write_permission = analysis.can_write(request.user)
    if not has_write_permission:
        set_form_read_only(form)

    reload_analysis = False
    reload_page = False
    if request.method == "POST":
        analysis.check_can_write(request.user)
        if form.has_changed:
            valid = form.is_valid()
            if valid:
                analysis = form.save()
                analysis.save()
                if reload_analysis := (old_annotation_version != analysis.annotation_version):
                    node_utils.reload_analysis_nodes(analysis.pk)
                # The panel layout and node anchors are rendered per orientation, and the save
                # turned the node positions - the open page needs to come back as the other mode
                reload_page = old_horizontal_mode != analysis.analysis_horizontal_mode

            add_save_message(request, valid, "Analysis Settings")

    for error in analysis.get_errors():
        messages.add_message(request, messages.ERROR, error)

    for warning in analysis.get_warnings():
        messages.add_message(request, messages.WARNING, warning)

    analysis_settings = get_analysis_settings(request.user, analysis)
    context = {"analysis": analysis,
               "form": form,
               "new_analysis_settings": analysis_settings,
               "has_write_permission": has_write_permission,
               "hidden_nodes": _get_hidden_nodes(analysis),
               "reload_analysis": reload_analysis,
               "reload_page": reload_page}
    return render(request, 'analysis/analysis_settings_details_tab.html', context)


def analysis_settings_node_counts_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    return analysis_settings_node_counts_tab_for(request, analysis, has_write_permission=analysis.can_write(request.user))


def analysis_settings_node_counts_tab_for(request, analysis, is_analysis=True, has_write_permission=True):
    """ analysis - an Analysis, or the FakeAnalysis (is_analysis=False) the user/lab/org default
        node counts settings pages use @see snpdb.views.views_user_settings._settings_override_node_counts_tab """
    if request.method == "POST":
        if has_write_permission is False:
            raise PermissionDenied()
        if is_analysis:
            auto_add_tags = bool(request.POST.get("node_count_auto_add_tags"))
            if auto_add_tags != analysis.node_count_auto_add_tags:
                analysis.node_count_auto_add_tags = auto_add_tags
                analysis.save()
        node_counts_str = request.POST.get("node_counts")
        if node_counts_str is not None:
            node_counts_array = node_counts_str.split(',')
            analysis.set_node_count_types(node_counts_array)
            if is_analysis:
                # Unlike tagging (where this is done in a task) the user is already waiting on this
                # form, and filling the counts in here means the tab comes back with them populated
                node_utils.update_analysis_tag_node_counts(analysis)
        add_save_message(request, True, "Node Counts")

    my_list, available_filters_list, available_tags_list = get_node_counts_mine_and_available(analysis)
    context = {"my_node_counts_list": my_list,
               "available_node_counts_list": available_filters_list,
               "available_tag_node_counts_list": available_tags_list,
               "is_analysis": is_analysis,
               "has_write_permission": has_write_permission}

    if is_analysis:
        context["node_count_auto_add_tags"] = analysis.node_count_auto_add_tags
        analysis_settings = get_analysis_settings(request.user, analysis)
        context["new_analysis_settings"] = analysis_settings

    return render(request, 'analysis/analysis_settings_node_counts_tab.html', context)


def analysis_settings_template_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)

    context = {"analysis": analysis}

    # Template run info (if this analysis was generated from a template)
    if hasattr(analysis, 'analysistemplaterun'):
        analysis_template_run = analysis.analysistemplaterun
        node_variables = defaultdict(list)
        for node in analysis.analysisnode_set.filter(analysisvariable__isnull=False).distinct().order_by("y"):
            for av in node.analysisvariable_set.all().order_by("field"):
                node_variables[node].append(av)
        context["analysis_template_run"] = analysis_template_run
        context["node_variables"] = defaultdict_to_dict(node_variables)

    # Create template form (if this analysis is not itself a template)
    if analysis.template_type != AnalysisTemplateType.TEMPLATE:
        form = forms.CreateAnalysisTemplateForm(request.POST or None, user=request.user, analysis=analysis)
        if request.method == "POST":
            if form.is_valid():
                analysis_template = form.save()
                return JsonResponse({"analysis_id": analysis_template.analysis_id})
        context["create_analysis_template_form"] = form

    return render(request, 'analysis/analysis_settings_template_tab.html', context)


class AnalysisLogEntryWrapper:
    def __init__(self, log_entry: LogEntry):
        self.log_entry = log_entry


def analysis_settings_audit_log_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    log_entry_wrappers = [
        AnalysisLogEntryWrapper(le) for le in analysis.log_entry_qs()
    ]

    context = {
        "analysis": analysis,
        "log_entry_wrappers": log_entry_wrappers
    }
    return render(request, 'analysis/analysis_settings_audit_log_tab.html', context)


@user_passes_test(is_superuser)
def analysis_settings_benchmark_tab(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    nodes_qs = analysis.analysisnode_set.all().select_subclasses()

    node_rows = []
    total_load_seconds = 0.0
    nodes_with_load = 0
    for node in nodes_qs:
        if node.load_seconds is not None:
            total_load_seconds += node.load_seconds
            nodes_with_load += 1
        node_rows.append({
            "pk": node.pk,
            "name": node.name or "",
            "node_type": type(node).__name__,
            "count": node.count,
            "load_seconds": node.load_seconds,
            "status": node.get_status_display(),
            "version": node.version,
        })
    node_rows.sort(key=lambda r: (r["load_seconds"] is None, -(r["load_seconds"] or 0.0)))

    nc_times = list(NodeCount.objects.filter(node_version__node__analysis=analysis)
                    .values_list("created", flat=True))
    if nc_times:
        wall_start = min(nc_times)
        wall_end = max(nc_times)
        wall_seconds = (wall_end - wall_start).total_seconds()
    else:
        wall_start = None
        wall_end = None
        wall_seconds = None

    context = {
        "analysis": analysis,
        "node_rows": node_rows,
        "total_load_seconds": total_load_seconds,
        "nodes_with_load": nodes_with_load,
        "node_total": len(node_rows),
        "wall_seconds": wall_seconds,
        "wall_start": wall_start,
        "wall_end": wall_end,
        "node_count_sample_size": len(nc_times),
    }
    return render(request, 'analysis/analysis_settings_benchmark_tab.html', context)


@require_POST
def analysis_settings_lock(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    if not analysis.can_unlock(request.user):  # check_can_write returns false if locked
        raise PermissionDenied(f"You do not have write access to {analysis.pk}")
    lock = json.loads(request.POST["lock"])
    AnalysisLock.objects.create(analysis=analysis, locked=lock, user=request.user, date=timezone.now())
    # Bump version to expire cache
    analysis.version += 1
    analysis.save()
    return redirect(analysis)  # Reload


def analysis_input_samples(request, analysis_id):
    analysis = get_analysis_or_404(request.user, analysis_id)
    input_nodes = analysis.analysisnode_set.filter(analysisnode_parent__isnull=True)
    input_nodes = input_nodes.select_subclasses()

    node_inputs = []
    all_samples = set()
    for node in input_nodes:
        if samples := node.get_samples_from_node_only_not_ancestors():
            all_samples.update(samples)
            node_input_data = {
                "node": node,
                "class_label": node.get_node_class_label(),
                "name": node.name or node.get_name_or_identifier(),
                "samples": samples,
            }
            for field in ["trio", "pedigree", "cohort"]:
                if source := getattr(node, field, None):
                    node_input_data["source"] = source
                    node_input_data["source_label"] = field.title()
                    break

            node_inputs.append(node_input_data)

    context = {"analysis": analysis,
               "node_inputs": node_inputs,
               "num_samples": len(all_samples)}
    return render(request, 'analysis/analysis_input_samples.html', context)
