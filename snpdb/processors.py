import socket
from functools import cached_property

from django.apps import apps
from django.conf import settings
from django.contrib.sites.models import Site

from snpdb.models import SiteMessage
from snpdb.user_settings_manager import UserSettingsManager
from uicore.utils.form_helpers import FORM_HELPER_HELPER
from variantgrid.perm_path import get_visible_url_names
from variantopedia.forms import SearchForm


class LazyUserProperties:
    """
    A class that we can send to every page that can provide access to common user values but doesn't need to
    pre-process anything unless the details are asked for
    """

    def __init__(self, request):
        self.request = request

    @cached_property
    def avatar_details(self):
        return UserSettingsManager.get_avatar_details()


def settings_context_processor(request):
    context = {
        'form_helper': FORM_HELPER_HELPER,
        'discordance_enabled': settings.DISCORDANCE_ENABLED,
        'classifications_new_grouping': settings.CLASSIFICATION_NEW_GROUPING,
        'help_url': settings.HELP_URL,
        'inbox_enabled': settings.INBOX_ENABLED,
        'pathology_tests_enabled': settings.PATHOLOGY_TESTS_ENABLED,
        'rollbar_access_token': settings.ROLLBAR.get('client_access_token'),
        'rollbar_enabled': settings.ROLLBAR.get('enabled', False),
        'rollbar_environment': settings.ROLLBAR.get('environment'),
        'sapath_enabled': apps.is_installed("sapath"),
        'seqauto_enabled': settings.SEQAUTO_ENABLED,
        'site': Site.objects.get_current(),
        'site_messages': SiteMessage.get_site_messages(),
        'site_name': settings.SITE_NAME,
        'timezone': settings.TIME_ZONE,
        'top_right_search_form': SearchForm(search_allow_blank=True),
        'url_name_visible': get_visible_url_names(),
        'use_oidc': settings.USE_OIDC,  # whether user is managed by django or externally by open connect
        'user_feedback_enabled': settings.ROLLBAR.get('enabled', False) and settings.USER_FEEDBACK_ENABLED,
        "contact_us_enabled": settings.CONTACT_US_ENABLED,
        "user_properties": LazyUserProperties(request)
    }

    # This can fail on bad urls
    if rm := request.resolver_match:
        context['url_name'] = rm.url_name

    if settings.DEBUG:
        context['debug'] = True
        context['hostname'] = socket.gethostname()

    if settings.SOMALIER.get("enabled"):
        context['somalier_enabled'] = request.user.is_superuser or not settings.SOMALIER.get("admin_only")

    if r_match := request.resolver_match:
        context['url_name'] = r_match.url_name

    # We extend templates to provide the menus
    # For clinicians, set them all to a restricted view with less menus
    MENU_BASE_TEMPLATES = [
        "menu_analysis_base",
        "menu_annotation_base",
        "menu_classifications_base",
        "menu_data_base",
        "menu_genes_base",
        "menu_help_base",
        "menu_pathtests_base",
        "menu_patients_base",
        "menu_settings_base",
        "menu_variants_base",
    ]

    DEFAULT_TEMPLATE_PATTERN = "snpdb/menu/%s.html"
    for base_template in MENU_BASE_TEMPLATES:
        context[base_template] = DEFAULT_TEMPLATE_PATTERN % base_template

    return context
