from django.apps import apps
from django.conf import settings
from django.conf.urls import include
from django.conf.urls.static import static
from django.contrib import admin
from django.contrib.auth.decorators import login_not_required, login_required
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.urls import path
from django.views.generic.base import TemplateView
from drf_spectacular.views import SpectacularAPIView, SpectacularRedocView, SpectacularSwaggerView

from library.django_utils.view_utils import login_not_required_include
from variantgrid import views
from variantgrid.views import ContactFormView, OneStepRegistrationView
from variantgrid.views_rest import CapabilitiesView

admin.autodiscover()
# Django marks the admin login login_not_required - send anonymous users to the site login (OIDC on Shariant) instead
admin.site.login = login_required(admin.site.login)

APPS_WITH_URLS = [
    "analysis",
    "annotation",
    "eventlog",
    "email_manager",
    "flags",
    "genes",
    "pathtests",
    "patients",
    "pedigree",
    "ontology",
    "sapath",
    "seqauto",
    "snpdb",
    "upload",
    "classification",
    "variantopedia",
    "manual",
    "review",
    "mme",
    "beacon",
]

urlpatterns = [
    path('', views.index),
    path('loading_animations', views.loading_animations, name='loading_animations'),
    path('admin/', admin.site.urls),
    path('authenticated', views.authenticated, name='authenticated'),
    path('external_help', views.external_help, name='external_help'),
    path('system/version', views.version, name='version'),
    path('system/changelog', views.changelog, name='changelog'),
    path('system/keycloak_admin', views.keycloak_admin, name='keycloak_admin'),
    path('terms/', include('termsandconditions.urls')),
    path('avatar/', include('avatar.urls')),
    path('api/schema', SpectacularAPIView.as_view(), name='openapi-schema'),
    path('api/docs', SpectacularSwaggerView.as_view(url_name='openapi-schema'), name='api-docs'),
    path('api/redoc', SpectacularRedocView.as_view(url_name='openapi-schema'), name='api-redoc'),
    path('api/v1/capabilities', CapabilitiesView.as_view(), name='api_capabilities'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

if settings.INBOX_ENABLED:
    urlpatterns += [path('messages/', include('user_messages.urls'))]

if settings.DEBUG:
    if 'debug_toolbar' in settings.INSTALLED_APPS:
        import debug_toolbar
        urlpatterns += [
            path('__debug__/', include(debug_toolbar.urls)),
        ]

if settings.CONTACT_US_ENABLED:
    urlpatterns += [
        path('contact_us', ContactFormView.as_view(), name='contact_us')
    ]

if getattr(settings, "REGISTRATION_OPEN", False):
    # registration.backends.simple.urls, but with our own view (see OneStepRegistrationView)
    urlpatterns += [
        path('accounts/register/closed/',
             login_not_required(TemplateView.as_view(template_name='registration/registration_closed.html')),
             name='registration_disallowed'),
        path('accounts/register/',
             login_not_required(OneStepRegistrationView.as_view(
                 success_url=getattr(settings, 'SIMPLE_BACKEND_REDIRECT_URL', '/'))),
             name='registration_register'),
        path('accounts/', login_not_required_include('registration.auth_urls')),
    ]
else:
    urlpatterns += [path('accounts/', login_not_required_include('registration.backends.default.urls'))]


if settings.USE_OIDC:
    urlpatterns += [
        path('oidc/', login_not_required_include('mozilla_django_oidc.urls')),
        path('oidc_login/', views.oidc_login),
    ]

handler404 = views.page_not_found
handler500 = views.server_error


for app_name in APPS_WITH_URLS:
    if apps.is_installed(app_name):
        if settings.URLS_APP_REGISTER[app_name]:
            app_urls = f"{app_name}.urls"
            urlpatterns.append(path(f"{app_name}/", include(app_urls)))

# Fix for gunicorn setup
urlpatterns += staticfiles_urlpatterns()
