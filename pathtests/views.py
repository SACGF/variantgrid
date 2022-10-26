import json
import logging
from collections import defaultdict

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models.query_utils import Q
from django.http.response import JsonResponse, HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.utils import timezone
from django.views.decorators.http import require_POST

from genes.models import GeneListCategory, GeneList, GeneSymbol
from genes.views import add_gene_list_unmatched_genes_message
from library.django_utils import add_save_message
from pathtests.forms import SelectPathologyTestForm, SelectPathologyTestVersionForm, \
    PathologyTestOrderForm, CaseForm, CreatePathologyTestForm, PathologyTestVersionForm
from pathtests.models import PathologyTest, PathologyTestVersion, \
    PathologyTestGeneModificationRequest, PathologyTestGeneModificationOutcome, get_cases_qs, PathologyTestOrder, \
    Case, get_external_order_system_last_checked, ActivePathologyTestVersion
from patients.forms import external_pk_autocomplete_form_factory
from patients.models import Clinician, get_lead_scientist_users_for_user, FollowLeadScientist
from snpdb.models.models_enums import ImportStatus


def user_is_clinician(user):
    try:
        user.clinician
        return True
    except Clinician.DoesNotExist:
        return False


def clinician_login(request):
    try:
        request.user.clinician
    except Clinician.DoesNotExist:
        msg = "Only clinicians can view this page!"
        raise PermissionDenied(msg)
    return redirect('pathology_test_requests')


def clinician_cases(request, clinician_id):
    clinician = get_object_or_404(Clinician, pk=clinician_id)
    q_clinician = Q(caseclinician__clinician=clinician)
    q_requested = Q(pathologytestorder__sapathologyrequestgenelistpathologytestorderlink__sapathology_request_gene_list__sapathology_request__user=clinician.user)
    cases_qs = get_cases_qs().filter(q_clinician | q_requested).order_by("-created")

    last_checked = get_external_order_system_last_checked()
    # No current way to find out cases by clinician - waiting on Helix to come back to us...
    context = {"last_checked": last_checked,
               "cases_qs": cases_qs}
    return render(request, 'pathtests/clinician_cases.html', context)


def my_cases(request, has_menu=True):
    users_list = get_lead_scientist_users_for_user(request.user)
    title = "My cases"
    return view_cases_for_users(request, title, users_list, has_menu=has_menu)


def scientist_cases(request, user_id, has_menu=True):
    user = get_object_or_404(User, pk=user_id)
    title = f"Cases for {user}"
    return view_cases_for_users(request, title, [user], has_menu=has_menu)


@require_POST
def follow_scientist(request, follow_user_id):
    checked = json.loads(request.POST["checked"])
    follow_user = get_object_or_404(User, pk=follow_user_id)

    if checked:
        FollowLeadScientist.objects.get_or_create(user=request.user, follow=follow_user)
    else:
        FollowLeadScientist.objects.filter(user=request.user, follow=follow_user).delete()

    return HttpResponse()


def view_cases_for_users(request, title, users_list, has_menu=True):
    cases_qs = get_cases_qs().filter(lead_scientist__in=users_list).order_by("-created")

    last_checked = get_external_order_system_last_checked()
    if has_menu:
        template_name = 'pathtests/cases/cases_menu.html'
    else:
        template_name = 'pathtests/cases/show_cases.html'

    is_following = False
    if len(users_list) == 1:
        single_lead_scientist = users_list[0]
        is_following = request.user.following.filter(follow=single_lead_scientist).exists()
    else:
        single_lead_scientist = None

    context = {"title": title,
               "single_lead_scientist": single_lead_scientist,
               "is_following": is_following,
               "is_clinician": False,
               "cases_qs": cases_qs,
               "last_checked": last_checked}

    return render(request, template_name, context)


def my_cases_tab(request):
    return my_cases(request, has_menu=False)


def pathology_tests(request):
    context = {}
    if settings.PATHOLOGY_TEST_EXTERNAL_CODE:
        external_test_order_form = external_pk_autocomplete_form_factory(settings.PATHOLOGY_TEST_EXTERNAL_CODE)
        last_checked = get_external_order_system_last_checked()
        context = {"external_test_order_form": external_test_order_form,
                   "last_checked": last_checked}
    return render(request, 'pathtests/pathology_tests.html', context)


def cases(request):
    external_code = settings.PATHOLOGY_TEST_CASE_EXTERNAL_CODE
    external_case_form = external_pk_autocomplete_form_factory(external_code)

    context = {"external_case_form": external_case_form}
    return render(request, 'pathtests/cases.html', context)


def pathology_test_requests(request):
    # This redirects to a custom page

    url = settings.PATHOLOGY_TEST_REQUESTS_REDIRECT_URL
    if not url:
        msg = "Your system does not have 'PATHOLOGY_TEST_REQUESTS_REDIRECT_URL' setting configured"
        raise PermissionDenied(msg)
    return redirect(url)


def manage_pathology_tests(request):
    create_pathology_test_form = CreatePathologyTestForm(request.POST or None)
    if request.method == "POST":
        valid = create_pathology_test_form.is_valid()
        if valid:
            name = create_pathology_test_form.cleaned_data['name']
            pathology_test, created = PathologyTest.objects.get_or_create(name=name, curator=request.user)
            if created:
                pathology_test_version = create_pathology_test_form.cleaned_data['pathology_test_version']
                if pathology_test_version:
                    gene_list = pathology_test_version.gene_list
                else:
                    gene_list = create_pathology_test_form.cleaned_data['gene_list']

                if gene_list:
                    gene_list = gene_list.clone()
                    gene_list.user = request.user
                    gene_list.name = f"Clone of {gene_list.name}"
                    gene_list.category = GeneListCategory.get_pathology_test_gene_category()
                    gene_list.save()
                else:
                    category = GeneListCategory.get_pathology_test_gene_category()
                    gene_list = GeneList.objects.create(name=name, user=request.user,
                                                        category=category,
                                                        import_status=ImportStatus.SUCCESS)

                PathologyTestVersion.objects.create(pathology_test=pathology_test,
                                                    version=1,
                                                    gene_list=gene_list)

                return redirect(pathology_test)

            msg = f"Pathology Test '{name}' already existed"
            messages.add_message(request, messages.ERROR, msg, extra_tags='save-message')

    context = {"pathology_test_form": SelectPathologyTestForm(),
               "pathology_test_version_form": SelectPathologyTestVersionForm(),
               "create_pathology_test_form": create_pathology_test_form}

    return render(request, 'pathtests/manage_pathology_tests.html', context)


def get_gene_modification_request(pathology_test_version):
    """ returns PathologyTestGeneModification separated as (gene_addition_requests, gene_deletion_requests, handled_requests) """
    current_genes = set()
    if pathology_test_version.gene_list:
        current_genes.update(pathology_test_version.gene_list.get_gene_names())
    #logging.info(current_genes)

    gene_addition_requests = defaultdict(list)
    gene_deletion_requests = defaultdict(list)
    handled_requests = defaultdict(list)

    for gmr in pathology_test_version.pathologytestgenemodificationrequest_set.all():
        gene_symbol = gmr.gene_symbol_id
        if gmr.outcome == PathologyTestGeneModificationOutcome.PENDING:
            if gene_symbol in current_genes:
                d = gene_deletion_requests
            else:
                d = gene_addition_requests
        else:
            d = handled_requests

        d[gene_symbol].append(gmr)

    # TODO: Go through and add summary - counts of add/del for gene

    # Cast to dict as Django can't iterate over defaultdict
    return dict(gene_addition_requests), dict(gene_deletion_requests), dict(handled_requests)


def handle_modification_requests(post_dict, pathology_test_version, op_key, modification_requests, genes):
    IGNORE = 'ignore'
    REJECT = 'reject'
    ACCEPT = 'accept'

    current_time = timezone.now()

    gene_modification_info = {}
    for gene_symbol, gene_modification_requests in modification_requests.items():
        # Looks like 'add-ALK' or 'del-ALK'
        key = f'{op_key}-{gene_symbol}'
        outcome = None
        op = post_dict.get(key)
        if op == IGNORE:
            continue

        if op == REJECT:
            outcome = PathologyTestGeneModificationOutcome.REJECTED
        elif op == ACCEPT:
            outcome = PathologyTestGeneModificationOutcome.ACCEPTED
            genes.add(gene_symbol)

            request_info = []
            for gmr in gene_modification_requests:
                req_info = f"{gmr.user} {gmr.get_operation_display()} on {gmr.modified}"
                if gmr.comments:
                    req_info += ": " + gmr.comments
                request_info.append(req_info)
            request_info_str = ", ".join(request_info)
            curator = pathology_test_version.pathology_test.curator
            modification_info = f"request: {request_info_str}, accepted by {curator} on {current_time}"
            gene_modification_info[gene_symbol] = modification_info

        for gmr in gene_modification_requests:
            gmr.pathology_test_version = pathology_test_version  # if new test, move over
            gmr.outcome = outcome
            gmr.save()

    return gene_modification_info


def apply_pathology_test_modifications(request, pathology_test_version, gene_addition_requests, gene_deletion_requests):
    gene_symbol_additions = set()
    gene_symbol_deletions = set()

    gami = handle_modification_requests(request.POST, pathology_test_version, "add", gene_addition_requests, gene_symbol_additions)
    handle_modification_requests(request.POST, pathology_test_version, "del", gene_deletion_requests, gene_symbol_deletions)

    if any([gene_symbol_additions, gene_symbol_deletions]):
        pathology_test_version.gene_list.add_and_remove_gene_symbols(gene_symbol_additions, gene_symbol_deletions, gene_additions_modification_info=gami)
        modifications = []
        if gene_symbol_additions:
            num_genes = len(gene_symbol_additions)
            genes_str = ', '.join(gene_symbol_additions)
            modifications.append(f"Added {num_genes} genes: {genes_str}")
        if gene_symbol_deletions:
            num_genes = len(gene_symbol_deletions)
            genes_str = ', '.join(gene_symbol_deletions)
            modifications.append(f"Removed {num_genes} genes: {genes_str}")

        modification_summary_message = '. '.join(modifications)
    else:
        modification_summary_message = "No changes performed."

    return modification_summary_message


def any_accepted_requests(request, gene_addition_requests, gene_deletion_requests):
    ACCEPT = 'accept'
    op_keys = {"add": gene_addition_requests,
               "del": gene_deletion_requests}

    for op_key, modification_requests in op_keys.items():
        for gene_symbol in modification_requests:
            # Looks like 'add-ALK' or 'del-ALK'
            key = f'{op_key}-{gene_symbol}'
            op = request.POST.get(key)
            if op == ACCEPT:
                return True

    return False


def view_pathology_test_version(request, pk):
    pathology_test_version = get_object_or_404(PathologyTestVersion, pk=pk)
    pathology_test_version_form = PathologyTestVersionForm(request.POST or None,
                                                           instance=pathology_test_version,
                                                           prefix='pathology-test-version')

    if request.method == 'POST':
        if not pathology_test_version.is_curator(request.user):
            msg = f"You are not the curator of pathology test {pathology_test_version.pathology_test}"
            raise PermissionError(msg)

        confirm_test = request.POST.get("confirm_test")
        if confirm_test:
            if not pathology_test_version.can_confirm:
                msg = "Can't set as active test!"
                raise ValueError(msg)

            pathology_test_version.confirmed_date = timezone.now()
            pathology_test_version.save()
            pathology_test_version.set_as_active_test()
            messages.add_message(request, messages.INFO, "Test Confirmed")
        else:
            valid = pathology_test_version_form.is_valid()
            if pathology_test_version_form.cleaned_data:  # Got something in POST
                if valid:
                    pathology_test_version = pathology_test_version_form.save()
                add_save_message(request, valid, "Pathology Test Version")
            else:
                gene_addition_requests, gene_deletion_requests, _ = get_gene_modification_request(pathology_test_version)
                test_confirmed = pathology_test_version.confirmed_date is not None
                new_test = test_confirmed and any_accepted_requests(request, gene_addition_requests, gene_deletion_requests)
                if new_test:
                    logging.info("New Test required")
                    pathology_test_version = pathology_test_version.next_version()
                modification_summary_message = apply_pathology_test_modifications(request, pathology_test_version, gene_addition_requests, gene_deletion_requests)

                if new_test:
                    messages.add_message(request, messages.INFO, f"Created New Test: {modification_summary_message}")
                    return redirect(pathology_test_version)
                messages.add_message(request, messages.INFO, modification_summary_message)

    other_active_test = None
    if not pathology_test_version.is_active_test:
        other_active_test = pathology_test_version.pathology_test.get_active_test_version()

    add_gene_list_unmatched_genes_message(request, pathology_test_version.gene_list)

    # If changed via post, will reload as need to shift from requests -> handled
    gene_addition_requests, gene_deletion_requests, handled_requests = get_gene_modification_request(pathology_test_version)

    context = {"pathology_test_version": pathology_test_version,
               "pathology_test_version_form": pathology_test_version_form,
               "other_active_test": other_active_test,
               "is_curator": pathology_test_version.is_curator(request.user),
               "gene_addition_requests": gene_addition_requests,
               "gene_deletion_requests": gene_deletion_requests,
               "handled_requests": handled_requests,
               "pathology_test_versions": [pathology_test_version]}
    return render(request, 'pathtests/view_pathology_test_version.html', context)


def view_pathology_test(request, name):
    pathology_test = get_object_or_404(PathologyTest, name=name)

    if request.method == 'POST':
        if not pathology_test.is_curator(request.user):
            msg = f"You are not the curator of pathology test {pathology_test}"
            raise PermissionError(msg)

        msg = None
        level = messages.INFO
        restore_test = request.POST.get("restore_test")
        if restore_test:
            pathology_test.deleted = False
            pathology_test.save()
            msg = "Test restored."

        delete_test = request.POST.get("delete_test")
        if delete_test:
            delete_text = request.POST.get("delete_text")
            if delete_text == "delete":
                pathology_test.deleted = True
                pathology_test.save()
                ActivePathologyTestVersion.objects.filter(pathology_test=pathology_test).delete()
                msg = "This test has been deleted."
            else:
                msg = "This test was NOT deleted."
                level = messages.WARNING

        if msg:
            messages.add_message(request, level, msg, extra_tags='save-message')

    gene_grid_columns = []
    for ptv in pathology_test.pathologytestversion_set.all():
        gene_grid_columns.append(f"pathology-test-version-{ptv.pk}")

    pathology_test_versions = pathology_test.pathologytestversion_set.all().order_by("-version")

    context = {"pathology_test": pathology_test,
               "is_curator": pathology_test.is_curator(request.user),
               "gene_grid_columns_from_url": "/".join(gene_grid_columns),
               "pathology_test_versions": pathology_test_versions}
    return render(request, 'pathtests/view_pathology_test.html', context)


@require_POST
def modify_pathology_test_version(request, pk):
    pathology_test_version = get_object_or_404(PathologyTestVersion, pk=pk)

    operation = request.POST['operation']
    symbol = request.POST['gene_symbol']
    comments = request.POST['comments']

    gene_symbol = get_object_or_404(GeneSymbol, pk=symbol)
    mod_request = PathologyTestGeneModificationRequest.objects.create(pathology_test_version=pathology_test_version,
                                                                      operation=operation,
                                                                      gene_symbol=gene_symbol,
                                                                      user=request.user,
                                                                      comments=comments)
    response_data = {"pk": mod_request.pk}
    return JsonResponse(response_data)


def add_external_manager_notice(request, obj):
    if obj.external_manager:
        explaination = obj.external_manager.explaination
        if explaination:
            messages.add_message(request, messages.INFO, explaination)


def view_pathology_test_order(request, pk):
    pathology_test_order = get_object_or_404(PathologyTestOrder, pk=pk)
    form = PathologyTestOrderForm(instance=pathology_test_order)
    add_external_manager_notice(request, pathology_test_order)

    context = {"pathology_test_order": pathology_test_order,
               "form": form}
    return render(request, 'pathtests/view_pathology_test_order.html', context)


def view_external_pathology_test_order(request, external_pk):
    pathology_test_order = get_object_or_404(PathologyTestOrder, external_pk=external_pk)
    return redirect(pathology_test_order)


def view_case(request, pk):
    case = get_object_or_404(Case, pk=pk)
    form = CaseForm(user=request.user, instance=case)

    add_external_manager_notice(request, case)
    context = {"case": case,
               "form": form,
               "can_write": case.can_write(request.user)}
    return render(request, 'pathtests/view_case.html', context)


def view_external_case(request, external_pk):
    case = get_object_or_404(Case, external_pk=external_pk)
    return redirect(case)
