from django import forms
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User
from django.contrib.sites.models import Site
from django.http.response import HttpResponseServerError, JsonResponse, \
    HttpResponseForbidden
from django.shortcuts import redirect, render
from django.template.loader import get_template, render_to_string
from django.urls.base import resolve, reverse
from django.urls.exceptions import Resolver404
from django.views.generic import FormView
from global_login_required import login_not_required

from library.email import Email
from library.git import Git
from library.keycloak import Keycloak, KeycloakError, KeycloakNewUser
from library.log_utils import report_exc_info, AdminNotificationBuilder
from manual.models import Deployment
from snpdb.forms import KeycloakUserForm
from snpdb.models import UserSettings


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
    git = Git(settings.BASE_DIR)

    deployments = list()
    is_first = True
    for deployment in Deployment.objects.order_by('-created').all()[0:10]:
        if is_first:
            if deployment.git_hash != git.hash:
                # git hash doesn't match most recent deploy
                # inject one to highlight that
                deployments.append({
                    "git_hash": git.hash,
                    "created": None,
                    "git_link": None,
                })
        is_first = False

        deployment_git_hash = deployment.git_hash
        deployment_git_link = None
        if git.site and git.hash and deployment_git_hash:
            if git.hash != deployment_git_hash:
                deployment_git_link = f"{git.site}/compare/{deployment_git_hash}...{git.hash}"
            else:
                # note: the links are for what's different between the current version and a previous deploy
                # if the hashes are the same, it's the same deploy - so only showing the changes in the most recent
                # commit is misleading
                # deployment_git_link = f"{git.site}/commit/{git.hash}"
                pass

        deployments.append({
            "git_hash": deployment.git_hash,
            "created": deployment.created,
            "git_link": deployment_git_link
        })

    weekly_update_users = list()
    if request.user.is_superuser:
        all_users = User.objects.filter(is_active=True, email__isnull=False).order_by('email')
        for user in all_users:
            # user email could be blank instead of null
            if user.email and UserSettings.get_for_user(user).email_weekly_updates:
                weekly_update_users.append(user)

    context = {
        "git": git,
        "deployment_history": deployments,
        "weekly_update_users": weekly_update_users
    }
    return render(request, 'version.html', context)


def changelog(request):
    return render(request, 'changelog.html')


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


class ContactUsForm(forms.Form):
    name = forms.CharField(max_length=200)
    email = forms.EmailField(max_length=200)
    message = forms.CharField(widget=forms.Textarea)


@login_not_required
class ContactFormView(FormView):
    form_class = ContactUsForm
    template_name = "registration/contact_us.html"

    def post(self, request, *args, **kwargs):
        form = ContactUsForm(data=request.POST)
        context = {"form": form}
        if form.is_valid():
            try:
                self.send_message(form)
                messages.success(request, "Your details have been sent")
                context["sent"] = True
            except:
                report_exc_info(request=request)
                # TODO, maybe grab an email from settings at this point?
                messages.error(request, "There was an error with the contact form. Please try again later")

        return render(request, self.template_name, context)

    def send_message(self, form):
        nb = AdminNotificationBuilder("Contact Us", is_communication=True)
        nb.add_header("Contact Us")
        nb.add_field("From", form.cleaned_data.get('name'))
        nb.add_field("Email", form.cleaned_data.get('email'))
        nb.add_markdown(form.cleaned_data.get('message'))
        nb.send()

"""
def contact_us(request):
    name = ""
    email = ""
    message = ""
    if request.method == "POST":
        name = request.POST.get("name")
        email = request.POST.get("email")
        message = request.POST.get("message")
        if not email:


    context = {
        "name": name,
        "email": email,
        "message": message
    }

    return render(request, 'registration/contact_us.html', {})
"""
