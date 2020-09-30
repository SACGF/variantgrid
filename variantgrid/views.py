from django.conf import settings
from django.contrib import messages
from django.contrib.sites.models import Site
from django.http.response import HttpResponseServerError, JsonResponse, \
    HttpResponseForbidden
from django.shortcuts import redirect, render
from django.template.loader import get_template, render_to_string
from django.urls.base import resolve, reverse
from django.urls.exceptions import Resolver404
from global_login_required import login_not_required

from library.email import Email
from library.keycloak import Keycloak, KeycloakError, KeycloakNewUser
from library.log_utils import report_exc_info
from library.utils import get_git_last_modified_date, get_git_hash
from snpdb.forms import KeycloakUserForm


@login_not_required
def oidc_login(request):
    """
    Go straight to OIDC login page (with next parameter included)
    """
    auth_view = resolve(reverse('oidc_authentication_init')).func
    return auth_view(request)


@login_not_required
def index(request):
    if request.user.is_authenticated:
        return redirect(settings.LOGIN_REDIRECT_URL)
    return redirect('auth_login')


@login_not_required
def page_not_found(request, exception):
    message = "The requested object could not be found"
    if isinstance(exception, Resolver404):  # don't want to show paths
        message = "Could not determine page."
        exception = None
    context = {"message": message, "exception": exception}
    return custom_error_view(request, context)


@login_not_required
def server_error(request):
    message = "There was an internal server error (sorry!). This has been logged and will be investigated by the development team. "
    context = {"message": message}
    return custom_error_view(request, context)


@login_not_required
def csrf_error(request, reason=''):
    message = "Your session token has changed. Please try your operation again."
    context = {"message": message}
    return custom_error_view(request, context)


@login_not_required
def custom_error_view(request, context=None):
    # Django 404 and 500 pages dont pass context to templates.
    # This is a hack to do so, as I want to use static files etc.

    if request.user.is_authenticated:
        template_name = "server_error_logged_in.html"
    else:
        template_name = "server_error_unauth.html"
    t = get_template(template_name)
    response = t.render(request=request, context=context)
    return HttpResponseServerError(response)


@login_not_required
def authenticated(request):
    return JsonResponse({"authenticated": request.user.is_authenticated})


def version(request):
    git_last_modified = get_git_last_modified_date(settings.BASE_DIR)

    context = {"git_hash": get_git_hash(settings.BASE_DIR),
               "git_last_modified": git_last_modified}

    return render(request, 'version.html', context)


#@cache_page(WEEK_SECS)
@login_not_required
def external_help(request):
    """ Made this a template so we can override it for branding """

    context = {}
    return render(request, 'external_help.html', context)


def keycloak_admin(request):
    if not request.user.is_superuser:
        return HttpResponseForbidden()

    response = ''
    if request.method == "POST":
        form = KeycloakUserForm(request.POST)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            first_name = cleaned_data['first_name']
            last_name = cleaned_data['last_name']
            email = cleaned_data['email']
            lab = cleaned_data['lab']
            ku = KeycloakNewUser(email=email, first_name=first_name, last_name=last_name, lab=lab)

            try:
                context = vars(ku)
                site = Site.objects.get_current(request)
                context['site'] = site
                context['accounts_email'] = settings.ACCOUNTS_EMAIL

                welcome_text = render_to_string('keycloak/welcome.txt', context=context)
                welcome_email = Email(
                    from_email=settings.ACCOUNTS_EMAIL,
                    to_email=email,
                    subject=f'Welcome to {site.name}',
                    text=welcome_text
                )
                ku.welcome_email = welcome_email

                Keycloak().add_user(ku)

                messages.success(request, 'User created and notified')
            except KeycloakError as kce:
                messages.error(request, str(kce))
            except Exception as e:
                report_exc_info(request=request)
                messages.error(request, str(e))

    else:
        form = KeycloakUserForm()

    context = {
        "form": form
    }

    return render(request, 'keycloak_admin.html', context)
