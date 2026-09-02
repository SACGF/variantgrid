import itertools
import json

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import Group, User
from django.forms.models import ALL_FIELDS, inlineformset_factory
from django.forms.widgets import TextInput
from django.http.response import (
    HttpResponse,
    HttpResponseRedirect,
)
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse
from django.views.decorators.http import require_POST

from analysis.models import Analysis
from analysis.views.views import _analysis_settings_node_counts_tab
from classification.models import Classification
from library.django_utils import (
    add_save_message,
)
from library.keycloak import Keycloak
from snpdb import forms
from snpdb.forms import (
    SettingsInitialGroupPermissionForm,
    UserContactForm,
    UserForm,
    UserSettingsOverrideForm,
)
from snpdb.models import (
    VCF,
    AbstractNodeCountSettings,
    AllVariantsFilter,
    AvatarDetails,
    CustomColumn,
    CustomColumnsCollection,
    GenomeBuild,
    Lab,
    LabUserSettingsOverride,
    NodeCountSettingsCollection,
    Organization,
    OrganizationUserSettingsOverride,
    Sample,
    UserContact,
    UserDataPrefix,
    UserGridConfig,
    UserSettings,
    UserSettingsOverride,
    VariantGridColumn,
)
from snpdb.models.models_enums import (
    BuiltInFilters,
)
from snpdb.sample_file_path import get_example_replacements
from snpdb.views.views_lab import _add_read_only_settings_message


@require_POST
def set_user_row_config(request):
    """ Rows per page, per user. Posted by DataTables when the page length changes, where the
        config opted in with a grid_name (@see DatatableConfig.grid_name) """

    grid_name = request.POST["grid_name"]
    grid_rows = int(request.POST["grid_rows"])

    UserGridConfig.objects.update_or_create(user=request.user, grid_name=grid_name, defaults={"rows": grid_rows})
    return HttpResponse()


@require_POST
def set_user_data_grid_config(request):
    """ This is set from datatables/user_data_grid_filter.html, should contain either filter_level+checked or filter_name """

    grid_name = request.POST["grid_name"]
    user_grid_config = UserGridConfig.get(request.user, grid_name)
    filter_level = request.POST.get("filter_level")
    if filter_level:
        checked = json.loads(request.POST["checked"])

        if filter_level == 'groups':
            user_grid_config.show_group_data = checked
        elif filter_level == 'incomplete':
            user_grid_config.show_incomplete_data = checked
        elif filter_level == 'hidden':
            user_grid_config.show_hidden_data = checked
        else:
            msg = f"Unknown value for filter_level: '{filter_level}'"
            raise ValueError(msg)
    else:
        user_grid_config.filter_name = request.POST["filter_name"]

    user_grid_config.save()
    return HttpResponse()


@require_POST
def set_show_tips(request):
    """ Turn off feature tips - posted by the tip box's dismiss button, see tips.js """

    UserSettingsOverride.objects.update_or_create(user=request.user, defaults={"show_tips": False})
    return HttpResponse()


@require_POST
def set_all_variants_filter(request, genome_build_name):
    """ Remember the All Variants page filter selections - set from variants.html whenever a filter changes """

    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    filters = json.loads(request.body)
    AllVariantsFilter.objects.update_or_create(user=request.user, genome_build=genome_build,
                                               defaults={"filters": filters})
    return HttpResponse()


def view_user_settings(request):
    user = request.user
    user_contact = UserContact.get_for_user(user)

    action = request.POST.get('action') if request.POST else None
    post = request.POST or None if not action else None
    user_form = UserForm(post, instance=user)
    user_contact_form = UserContactForm(post, instance=user_contact)
    user_settings = UserSettings.get_for_user(user)
    override_source, override_values = user_settings.get_override_source_and_values_before_user()
    user_settings_override = UserSettingsOverride.objects.get(user=user)
    user_settings_override_form = UserSettingsOverrideForm(post, instance=user_settings_override)
    labs_by_group_name = {l.group_name: l for l in Lab.valid_labs_qs(user)}

    group_initial_perm_forms = {}
    if settings.USER_SETTINGS_SHOW_GROUPS:
        read_groups, write_groups = user_settings.initial_perm_read_and_write_groups
        for group in user.groups.all().order_by('name'):
            initial = {"read": group in read_groups, "write": group in write_groups}
            group_initial_perm_forms[group] = SettingsInitialGroupPermissionForm(request.POST or None, initial=initial,
                                                                                 settings_override=user_settings_override,
                                                                                 group=group)
    if request.method == "POST":
        all_valid = True

        action = request.POST.get('action')
        if action == 'password-reset':
            keycloak = Keycloak()
            keycloak.change_password(user)
            messages.add_message(request, level=messages.INFO, message='Password reset email sent',
                                 extra_tags='save-message')
        else:
            if not settings.USE_OIDC:
                if user_form.is_valid():
                    user = user_form.save()
                else:
                    all_valid = False

            for form in itertools.chain([user_contact_form, user_settings_override_form],
                                        group_initial_perm_forms.values()):
                if form.is_valid():
                    form.save()
                else:
                    all_valid = False
            add_save_message(request, all_valid, "User Settings")

    context = {
        'user': user,
        'user_form': user_form,
        'user_contact_form': user_contact_form,
        'user_settings_form': user_settings_override_form,
        'group_initial_perm_forms': group_initial_perm_forms,
        'override_source': override_source,
        'override_values': override_values,
        'labs_by_group_name': labs_by_group_name,
        'avatar_details': AvatarDetails.avatar_for(user)
    }
    return render(request, 'snpdb/settings/view_user_settings.html', context)


def user_settings_node_counts_tab(request):
    user_settings_override = UserSettingsOverride.objects.get_or_create(user=request.user)[0]
    return _settings_override_node_counts_tab(request, user_settings_override)


def lab_settings_node_counts_tab(request, pk):
    lab = get_object_or_404(Lab, pk=pk)
    has_write_permission = lab.can_write(request.user)
    if has_write_permission is False:
        _add_read_only_settings_message(request, [lab])
    lab_settings_override = LabUserSettingsOverride.objects.get_or_create(lab=lab)[0]
    return _settings_override_node_counts_tab(request, lab_settings_override, has_write_permission=has_write_permission)


def organization_settings_node_counts_tab(request, pk):
    organization = get_object_or_404(Organization, pk=pk)
    has_write_permission = organization.can_write(request.user)
    if has_write_permission is False:
        _add_read_only_settings_message(request, organization.lab_set.all())

    org_settings_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=organization)[0]
    return _settings_override_node_counts_tab(request, org_settings_override, has_write_permission=has_write_permission)


def _settings_override_node_counts_tab(request, settings_override, has_write_permission=True):
    # This calls _analysis_settings_node_counts_tab with a FakeAnalysis object that
    # handles loading/saving a global one against User settings objects instead of analysis
    class FakeAnalysis:

        def set_node_count_types(self, node_counts_array):
            collection, _ = NodeCountSettingsCollection.objects.get_or_create(settings=settings_override)
            AbstractNodeCountSettings.save_count_configs_from_array(collection.nodecountsettings_set, node_counts_array)

        def get_node_count_types(self):
            try:
                node_count_config = settings_override.nodecountsettingscollection
                node_count_filters = node_count_config.get_node_count_filters()
            except:
                node_count_filters = BuiltInFilters.DEFAULT_NODE_COUNT_FILTERS

            return AbstractNodeCountSettings.get_types_from_labels(node_count_filters)

    fake_analysis = FakeAnalysis()
    return _analysis_settings_node_counts_tab(request, fake_analysis,
                                              is_analysis=False, has_write_permission=has_write_permission)


def view_user(request, pk):
    user = get_object_or_404(User, pk=pk)
    user_contact = UserContact.get_for_user(user)
    common_groups = Group.objects.filter(user=user.pk).filter(user=request.user.pk).order_by("name")
    all_groups = user.groups.order_by("name") if request.user.is_superuser else None

    context = {
        "other_user": user,
        'user_contact': user_contact,
        'common_groups': common_groups,
        'all_groups': all_groups,
        'genome_build': UserSettings.get_genome_build_or_default(request.user),
        'classifications_count': Classification.filter_for_user(request.user).filter(user=user).count(),
        'vcfs_count': VCF.filter_for_user(request.user).filter(user=user).count(),
        'samples_count': Sample.filter_for_user(request.user).filter(vcf__user=user).count(),
        'analyses_count': Analysis.filter_for_user(request.user).filter(template_type__isnull=True,
                                                                        user=user).count(),
    }
    return render(request, 'snpdb/settings/view_user.html', context)


def view_group(request, pk):
    group = get_object_or_404(Group, pk=pk)
    lab = Lab.objects.filter(group_name=group.name).first()  # group_name is unique so at most only 1
    users_qs = group.user_set.all().order_by("username")
    user_is_in_group = users_qs.filter(pk=request.user.pk).exists()
    public_groups = (settings.PUBLIC_GROUP_NAME, settings.LOGGED_IN_USERS_GROUP_NAME)

    context = {
        "group": group,
        "lab": lab,
        "public_group": group.name in public_groups,
        "user_is_in_group": user_is_in_group,
        "users_qs": users_qs,
    }
    return render(request, 'snpdb/settings/view_group.html', context)


def custom_columns(request):
    context = {}
    form = forms.CustomColumnsCollectionForm(request.POST or None, user=request.user)

    if request.method == "POST":
        if form.is_valid():
            ccc = form.save()
            return HttpResponseRedirect(reverse("view_custom_columns", kwargs={"custom_columns_collection_id": ccc.pk}))

        add_save_message(request, False, "Columns", created=True)

    context["form"] = form
    return render(request, 'snpdb/settings/custom_columns.html', context)


def view_custom_columns(request, custom_columns_collection_id):
    ccc = CustomColumnsCollection.get_for_user(request.user, custom_columns_collection_id)

    # One pass over the catalogue with the composites prefetched - the page draws every column, and
    # each card needs to know whether it is a composite or is drawn inside one
    variant_grid_columns = {vgc.pk: vgc for vgc in
                            VariantGridColumn.objects.prefetch_related("composite_members__column")}
    VariantGridColumn.annotate_composite_membership(variant_grid_columns.values())
    my_column_ids = list(ccc.customcolumn_set.order_by("sort_order").values_list("column_id", flat=True))
    my_columns = [variant_grid_columns[column_id] for column_id in my_column_ids]
    available_columns = [vgc for pk, vgc in variant_grid_columns.items() if pk not in set(my_column_ids)]

    has_write_permission = ccc.can_write(request.user)
    if not has_write_permission:
        msg = "You do not have permission to edit these columns. " \
              "If you wish to customise them, click 'clone' and modify the copy"
        messages.add_message(request, messages.WARNING, msg)

    if request.method == "POST":
        ccc.check_can_write(request.user)
        if name := request.POST.get("name"):
            ccc.name = name
            ccc.save()
        elif my_columns_str := request.POST.get("columns"):
            def update_user_columns(id_list, active):
                for i, col in enumerate(id_list):
                    column = variant_grid_columns[col]
                    CustomColumn.objects.update_or_create(custom_columns_collection=ccc, column=column,
                                                          defaults={"sort_order": i})
                # Delete any not in id_list
                CustomColumn.objects.filter(custom_columns_collection=ccc).exclude(column__in=id_list).delete()

            my_columns_list = my_columns_str.split(',') if my_columns_str else []
            active = 'my_columns' in request.POST
            update_user_columns(my_columns_list, active)
        return HttpResponse()  # Nobody ever looks at this

    context_dict = {
        'available_columns_list': available_columns,
        'my_columns_list': my_columns,
        'custom_columns': ccc,
        'has_write_permission': has_write_permission,
    }
    return render(request, 'snpdb/settings/view_custom_columns.html', context_dict)


def igv_integration(request):
    widgets = {"prefix": TextInput(attrs={'placeholder': 'from...'}),
               "replacement": TextInput(attrs={'placeholder': 'to...'})}

    UserDataPrefixFormSet = inlineformset_factory(User,
                                                  UserDataPrefix,
                                                  can_delete=True,
                                                  fields=ALL_FIELDS,
                                                  widgets=widgets,
                                                  max_num=10,
                                                  extra=3)
    if request.method == "POST":
        formset = UserDataPrefixFormSet(request.POST, instance=request.user)
        valid = formset.is_valid()
        if valid:
            formset.save()
        add_save_message(request, valid, "IGV Integration")

    formset = UserDataPrefixFormSet(instance=request.user)
    context_dict = {'user': request.user,
                    'formset': formset,
                    'example_replacements': get_example_replacements(request.user)}
    return render(request, 'snpdb/settings/igv_integration.html', context_dict)
