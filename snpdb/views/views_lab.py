from collections.abc import Iterable

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404, redirect, render
from django.urls.base import reverse
from django.views.decorators.http import require_POST

from classification.classification_stats import get_grouped_classification_counts
from classification.enums import AlleleOriginBucket
from classification.models.clinvar_export_sync import clinvar_export_sync
from classification.views.classification_accumulation_graph import (
    AccumulationReportMode,
    get_accumulation_graph_data,
)
from library.django_utils import (
    add_save_message,
    get_model_fields,
    set_form_read_only,
)
from snpdb.forms import (
    LabForm,
    LabMemberForm,
    LabUserSettingsOverrideForm,
    OrganizationForm,
    OrganizationUserSettingsOverrideForm,
    SettingsInitialGroupPermissionForm,
)
from snpdb.models import (
    ClinVarKey,
    Lab,
    LabHead,
    LabUserSettingsOverride,
    Organization,
    OrganizationUserSettingsOverride,
    State,
    UserSettings,
)
from snpdb.utils import LabNotificationBuilder


def _add_read_only_settings_message(request, lab_list: Iterable[Lab]):
    """ lab_list: labs where lab heads can modify settings """
    lab_heads_qs = LabHead.objects.filter(lab__in=lab_list).distinct()
    lab_head_names = ", ".join([str(lh.user) for lh in lab_heads_qs])
    if lab_head_names:
        lab_head_msg = f" or lab heads: {lab_head_names}"
    else:
        lab_head_msg = ""
    read_only_message = f"Only administrators{lab_head_msg} can modify these settings"
    messages.add_message(request, messages.INFO, read_only_message)


def view_lab(request, lab_id: int):
    lab = get_object_or_404(Lab, pk=lab_id)

    lab_form = LabForm(request.POST or None, instance=lab)
    lab_settings_override = LabUserSettingsOverride.objects.get_or_create(lab=lab)[0]

    override_fields = set(get_model_fields(LabUserSettingsOverride)) - {"id", "settingsoverride_ptr", "lab"}
    parent_overrides = UserSettings.get_settings_overrides(organization=lab.organization)
    override_source, override_values = UserSettings.get_override_source_and_values(override_fields, parent_overrides)
    settings_overrides = parent_overrides + [lab_settings_override]
    read_groups, write_groups = UserSettings.get_initial_perm_read_and_write_groups([lab.group], settings_overrides)

    initial = {"read": lab.group in read_groups, "write": lab.group in write_groups}
    group_initial_perm_form = None
    if settings.USER_SETTINGS_SHOW_GROUPS:
        group_initial_perm_form = SettingsInitialGroupPermissionForm(request.POST or None, initial=initial,
                                                                 settings_override=lab_settings_override,
                                                                 group=lab.group)
    lab_settings_override_form = LabUserSettingsOverrideForm(request.POST or None, instance=lab_settings_override)

    has_write_permission = lab.can_write(request.user)
    all_forms = [form for form in [lab_form, group_initial_perm_form, lab_settings_override_form] if form]

    if request.method == "POST":
        lab.check_can_write(request.user)

        if debug_method := request.POST.get("debug_method"):
            if "Test Slack" == debug_method:
                if not lab.slack_webhook:
                    messages.add_message(request, messages.ERROR, "Slack URL not configured correctly")
                else:
                    #try:
                    notification_builder = LabNotificationBuilder(lab=lab, message="Testing Slack Integration", notification_type=LabNotificationBuilder.NotificationType.SLACK_ONLY)
                    notification_builder.add_header(f"{settings.SITE_NAME} -> Slack Integration Test")
                    notification_builder.add_markdown("If you can see this, then integration has worked! :smile:")
                    notification_builder.send()
                    messages.add_message(request, messages.SUCCESS, "Message sent, check your Slack to confirm")
                    #except:
                    #    report_exc_info()
                    #    messages.add_message(request, messages.ERROR, "Unable to send test notification")
                return redirect(reverse('view_lab', kwargs={"lab_id": lab_id}))
            else:
                raise ValueError(f"Un-supported debug method {debug_method}")
        else:
            all_valid = True
            for form in all_forms:
                if form.is_valid():
                    form.save()
                else:
                    all_valid = False
            add_save_message(request, all_valid, "Lab Settings")

    if has_write_permission is False:
        for form in all_forms:
            set_form_read_only(form)
        # we just hide the form now
        # _add_read_only_settings_message(request, [lab])

    if settings.CLASSIFICATION_STATS_USE_SHARED:
        visibility = "Shared"
    else:
        visibility = "Created"

    context = {
        "lab": lab,
        "visibility": visibility,
        "is_member": lab.is_member(request.user) or request.user.is_superuser,
        "can_manage_members": lab.can_manage_members(request.user),

        "lab_form": lab_form,
        'settings_override_form': lab_settings_override_form,
        'group_initial_perm_form': group_initial_perm_form,
        'override_source': override_source,
        'override_values': override_values,
        'has_write_permission': has_write_permission,
        'clinvar_export_enabled': clinvar_export_sync.is_enabled
    }
    return render(request, 'snpdb/settings/view_lab.html', context)


def lab_members_tab(request, pk):
    lab = get_object_or_404(Lab, pk=pk)
    lab.check_can_manage_members(request.user)

    lab_users = lab.lab_users
    removable_user_pks = {lu.user.pk for lu in lab_users
                          if not lab.remove_member_blocked_reason(lu.user, request.user)}

    context = {
        "lab": lab,
        "lab_users": lab_users,
        "removable_user_pks": removable_user_pks,
        "lab_member_form": LabMemberForm(lab=lab, for_user=request.user),
    }
    return render(request, 'snpdb/settings/lab_members_tab.html', context)


@require_POST
def lab_add_member(request, pk):
    lab = get_object_or_404(Lab, pk=pk)
    lab.check_can_manage_members(request.user)

    form = LabMemberForm(request.POST, lab=lab, for_user=request.user)
    if form.is_valid():
        user = form.cleaned_data["user"]
        lab.add_member(user, added_by=request.user)
        messages.add_message(request, messages.SUCCESS, f"{user} added to {lab}")
    else:
        messages.add_message(request, messages.ERROR, "Could not add that user to the lab")
    return redirect(reverse('view_lab', kwargs={"lab_id": lab.pk}))


@require_POST
def lab_remove_member(request, pk):
    lab = get_object_or_404(Lab, pk=pk)
    lab.check_can_manage_members(request.user)

    members_qs = lab.group.user_set.all() if lab.group else User.objects.none()
    user = get_object_or_404(members_qs, pk=request.POST.get("user") or 0)
    lab.remove_member(user, removed_by=request.user)
    messages.add_message(request, messages.SUCCESS, f"{user} removed from {lab}")
    return redirect(reverse('view_lab', kwargs={"lab_id": lab.pk}))


def view_clinvar_key(request, pk: str):
    clinvar_key = get_object_or_404(ClinVarKey, pk=pk)
    clinvar_key.check_user_can_access(request.user)

    return render(request, 'snpdb/settings/clinvar_key.html', {
        'clinvar_key': clinvar_key,
        'labs': Lab.objects.filter(clinvar_key=clinvar_key).order_by('name')
    })


def view_organization(request, organization_id: int):
    organization = get_object_or_404(Organization, pk=organization_id)
    organization_form = OrganizationForm(request.POST or None, instance=organization)
    org_settings_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=organization)[0]
    override_fields = set(get_model_fields(OrganizationUserSettingsOverride)) - {"id", "settingsoverride_ptr", "organization"}
    parent_overrides = UserSettings.get_settings_overrides()
    override_source, override_values = UserSettings.get_override_source_and_values(override_fields, parent_overrides)
    org_settings_override_form = OrganizationUserSettingsOverrideForm(request.POST or None, instance=org_settings_override)
    all_forms = [organization_form, org_settings_override_form]

    if request.method == "POST":
        organization.check_can_write(request.user)
        all_valid = True
        for form in all_forms:
            if form.is_valid():
                form.save()
            else:
                all_valid = False
        add_save_message(request, all_valid, "Organization Settings")

    has_write_permission = organization.can_write(request.user)
    if has_write_permission is False:
        for form in all_forms:
            set_form_read_only(form)
        # put on individual tabs now
        # _add_read_only_settings_message(request, organization.lab_set.all())

    context = {
        "organization": organization,
        "is_member": organization.is_member(request.user) or request.user.is_superuser,
        "organization_form": organization_form,
        'settings_override_form': org_settings_override_form,
        'override_source': override_source,
        'override_values': override_values,
        'has_write_permission': has_write_permission,
    }
    return render(request, 'snpdb/settings/view_organization.html', context)


def labs(request):
    # Use short names is available
    show_unclassified = not settings.CLASSIFICATION_STATS_USE_SHARED

    """
    vc_state_data_json = get_grouped_classification_counts(
        user=request.user,
        field=state_field,
        max_groups=15,
        show_unclassified=show_unclassified)
    """
    active_organizations = Organization.objects.filter(active=True).order_by('name')
    organization_labs = {}
    for org in active_organizations:
        if settings.CLASSIFICATION_STATS_USE_SHARED:
            org_labs = org.sharing_labs
        else:
            org_labs = org.classifying_labs
        if org_labs:
            organization_labs[org] = list(org_labs)
    lab_list = [l for ll in organization_labs.values() for l in ll]
    context = {
        "organization_labs": organization_labs,
        "labs": lab_list,
        "shared_classifications": settings.CLASSIFICATION_STATS_USE_SHARED,
        # "vc_state_data": vc_state_data_json,
        "show_unclassified": show_unclassified,
    }

    return render(request, "snpdb/labs.html", context)


def labs_graph_detail(request):
    short_names_qs = Organization.objects.filter(short_name__isnull=False)
    name_to_short_name = dict(short_names_qs.values_list("name", "short_name"))
    state_field = "classification__lab__state"
    show_unclassified = not settings.CLASSIFICATION_STATS_USE_SHARED
    org_field = "classification__lab__organization__name"

    germline_buckets = {AlleleOriginBucket.GERMLINE, AlleleOriginBucket.UNKNOWN}

    vc_org_data_json = get_grouped_classification_counts(
        user=request.user,
        field=org_field,
        max_groups=15,
        field_labels=name_to_short_name,
        show_unclassified=show_unclassified,
        allele_level=True,
        allele_origin_buckets=germline_buckets)

    context = {}
    context["vc_org_data"] = vc_org_data_json

    graph_data = get_accumulation_graph_data(mode=AccumulationReportMode.Allele)
    context["accumulation_by_status"] = graph_data["status"]
    if request.user.is_superuser:
        # TODO, do we really need to hide this graph away?
        context["accumulation_by_lab"] = graph_data["lab"]

        state_pop_multiplier = {}
        for state in State.objects.filter(population__gt=0):
            state_pop_multiplier[state.name] = 100_000 / state.population

        vc_normalized_state_data_json = get_grouped_classification_counts(
            user=request.user,
            field=state_field,
            max_groups=15,
            show_unclassified=show_unclassified,
            norm_factor=state_pop_multiplier,
            allele_level=True,
            allele_origin_buckets=germline_buckets)

        context["vc_normalized_state_data_json"] = vc_normalized_state_data_json
    return render(request, "snpdb/labs_graph_detail.html", context)
