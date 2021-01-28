from django.conf import settings
from django.template.library import Library
from classification.models.classification_utils import UserClassificationStats

register = Library()

"""
todo, put in support for

    if Clinician.user_is_clinician(request.user):
        for base_template in MENU_BASE_TEMPLATES:
            context[base_template] = "snpdb/clinician_view_base.html"
            
also, put in support for

    def get_sequencing_software_versions_template():
        if settings.SEQAUTO_ENABLED:
            base_template = "seqauto/menu_sequencing_data_base.html"
        else:
            base_template = "snpdb/menu/menu_settings_base.html"

"""


@register.inclusion_tag("uicore/menus/menu_bar_annotations.html", takes_context=True)
def menu_bar_annotations(context):
    return {}


@register.inclusion_tag("uicore/menus/menu_bar_analysis.html", takes_context=True)
def menu_bar_analysis(context):
    return context


@register.inclusion_tag("uicore/menus/menu_bar_main.html", takes_context=True)
def menu_bar_main(context):
    return context


@register.inclusion_tag("uicore/menus/menu_bar_data.html", takes_context=True)
def menu_bar_data(context):
    return {}


@register.inclusion_tag("uicore/menus/menu_bar_genes.html", takes_context=True)
def menu_bar_genes(context):
    return {}


@register.inclusion_tag("uicore/menus/menu_bar_patients.html", takes_context=True)
def menu_bar_patients(context):
    return {}


@register.inclusion_tag("uicore/menus/menu_bar_tests.html", takes_context=True)
def menu_bar_tests(context):
    return {}


@register.inclusion_tag("uicore/menus/menu_bar_sequencing.html", takes_context=True)
def menu_bar_sequencing(context):
    return {}


@register.inclusion_tag("uicore/menus/menu_bar_settings.html", takes_context=True)
def menu_bar_settings(context):
    return {
        'seqauto_enabled': settings.SEQAUTO_ENABLED
    }


@register.inclusion_tag("uicore/menus/menu_bar_classifications.html", takes_context=True)
def menu_bar_classifications(context):
    return {
        'classification_issue_count': UserClassificationStats(context.request.user).issue_count
    }


@register.inclusion_tag("uicore/menus/menu_bar_variants.html", takes_context=True)
def menu_bar_variants(context):
    return {}


@register.inclusion_tag("uicore/site_messages/site_messages.html")
def site_messages(site_messages):
    return {"site_messages": site_messages}
