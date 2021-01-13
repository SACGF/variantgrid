from django.apps import apps
from django.conf import settings
from django.conf.urls import include, url
from django.conf.urls.static import static
from django.contrib import admin
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.urls import path
import debug_toolbar

from variantgrid import views


admin.autodiscover()

APPS_WITH_URLS = ["analysis", "annotation", "eventlog",
                  "expression", "flags", "genes", "pathtests", "patients", "pedigree", "ontology",
                  "sapath", "seqauto", "snpdb", "upload", "classification", "variantopedia"]

urlpatterns = [
    path('', views.index),
    path('admin/', admin.site.urls),
    path('accounts/', include('registration.backends.default.urls')),
    path('authenticated', views.authenticated, name='authenticated'),
    path('messages/', include('django_messages.urls')),
    path('external_help', views.external_help, name='external_help'),
    path('version', views.version, name='version'),
    path('changelog', views.changelog, name='changelog'),
    path('keycloak_admin', views.keycloak_admin, name='keycloak_admin'),
    path('terms/', include('termsandconditions.urls')),
    path('__debug__/', include(debug_toolbar.urls)),
    url('avatar/', include('avatar.urls')),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

if settings.USE_OIDC:
    urlpatterns += [
        path('oidc/', include('mozilla_django_oidc.urls')),
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
