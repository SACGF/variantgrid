from dataclasses import dataclass
from typing import List, Any, Optional

from django import template
from django.http import HttpRequest
from django.template.base import FilterExpression
from django.utils.text import slugify

from uicore.templatetags.ui_utils import parse_tag, TagUtils

register = template.Library()


@dataclass
class TabBuilderTab:
    tab_builder: 'TabBuilder'
    tab_number: int
    tab_id: str
    label: str
    badge: Optional[int] = None
    badge_status: Optional[str] = None
    admin_only: Optional[bool] = None
    url: Optional[str] = None
    resolved_url: Optional[str] = None
    param: Any = None
    content: str = ''

    @property
    def active(self):
        return self.tab_builder.active_tab == self.tab_number


class TabBuilder:

    def __init__(self, tab_set: str):
        self.tab_set = tab_set
        self.rendered = False
        self.active_tab = 0
        self.tabs: List[TabBuilderTab] = list()

    def __del__(self):
        if not self.rendered:
            print("Generated TabBuilder but didn't render it with ui_render_tabs")

    @property
    def tabs_required(self) -> bool:
        """
        Returns if we need tabs, or if it can be rendered directly on the page.
        Currently if there is just one tab but it's an ajax tab, this will return True
        :return:
        """
        if len(self.tabs) > 1:
            return True
        if len(self.tabs) == 0:
            return False
        return self.tabs[0].url is not None


def check_active_tab(tab_set: str, tab_id: str, request: HttpRequest) -> bool:
    active = False
    if activeTab := request.GET.get("activeTab"):
        parts = activeTab.split(":")
        if len(parts) == 2 and parts[0] == tab_set:
            active = parts[1] == tab_id
        elif len(parts) == 1:
            active = parts[0] == tab_id
    return active


@register.simple_tag(takes_context=True)
def ui_register_tab(
        context,
        tab_set: str,
        label: str,
        url: str = None,
        param: Any = None,
        url_check=False,
        badge: Optional[int] = None,
        badge_status: Optional[str] = None,
        active=False):

    if url_check:
        if not context["url_name_visible"].get(url):
            return ""

    tab_key = f"ui-tab-{tab_set}"
    builder: TabBuilder = context.get(tab_key)
    if not builder:
        builder = TabBuilder(tab_set)
        context[tab_key] = builder

    tab_number = len(builder.tabs)
    param_id = "_" + str(param) if param else ""
    tab_id = url + param_id

    if active or check_active_tab(tab_set, tab_id, context.request):
        builder.active_tab = tab_number

    builder.tabs.append(TabBuilderTab(tab_builder=builder, tab_number=tab_number,
                                      tab_id=url + param_id, label=label, badge=badge, badge_status=badge_status, url=url, param=param))
    return ""


@register.simple_tag(takes_context=True)
def ui_register_tabs(context, tab_set: str):
    """
    Used because the creation of he tabs in subcontext e.g. for loops, will dissapear when poping to a higher context
    """
    tab_key = f"ui-tab-{tab_set}"
    context[tab_key] = TabBuilder(tab_set)
    return ""


@register.inclusion_tag("uicore/tags/tabs.html", takes_context=True)
def ui_render_tabs(context, tab_set: str, css: str = None, mode='tabs'):
    tab_key = f"ui-tab-{tab_set}"
    tab_builder: TabBuilder = context.get(tab_key)
    if tab_builder is None:
        # TODO, maybe include this text just as text if in debug mode
        raise ValueError('No TabBuilder found - if defining tabs within a loop or other sub-context, include ui_register_tabs beforehand at the desired context level')
    tab_builder.rendered = True
    return {
        "tab_builder": tab_builder,
        "id": tab_set,
        "css": css or "",
        "mode": mode
    }


@register.tag(name='ui_register_tab_embedded')
def ui_register_tab_embedded(parser, token):
    tag_name, args, kwargs = parse_tag(token, parser)
    nodelist = parser.parse(('end_ui_register_tab_embedded',))
    parser.delete_first_token()
    return LocalTabContent(
        nodelist,
        tab_set=kwargs.get('tab_set'),
        label=kwargs.get('label'),
        admin_only=kwargs.get('admin_only'),
        badge=kwargs.get('badge'),
        badge_status=kwargs.get('badge_status')
    )


class LocalTabContent(template.Node):
    def __init__(self,
                 nodelist,
                 tab_set: FilterExpression,
                 label: FilterExpression,
                 admin_only: FilterExpression,
                 badge: FilterExpression,
                 badge_status: FilterExpression):
        self.nodelist = nodelist
        self.tab_set = tab_set
        self.label = label
        self.admin_only = admin_only
        self.badge = badge
        self.badge_status = badge_status

    def render(self, context):
        admin_only = TagUtils.value_bool(context, self.admin_only)
        tab_set = TagUtils.value_str(context, self.tab_set)
        label = TagUtils.value_str(context, self.label)
        badge = TagUtils.value_int(context, self.badge)
        badge_status = TagUtils.value_str(context, self.badge_status)
        if not tab_set:
            raise ValueError("UI Tab requires a value for 'tab_set'")
        if not label:
            raise ValueError("UI Tab requires a value for 'label'")

        # Needs to not start with a number
        tab_id = 't' + slugify(label)
        tab_key = f"ui-tab-{tab_set}"
        builder: TabBuilder = context.get(tab_key)
        if not builder:
            builder = TabBuilder(tab_set)
            context[tab_key] = builder

        tab_number = len(builder.tabs)
        content: str = self.nodelist.render(context)

        if admin_only and not context.request.user.is_superuser:
            return

        if check_active_tab(tab_set, tab_id, context.request):
            builder.active_tab = tab_number

        if content.startswith('/'):
            builder.tabs.append(TabBuilderTab(tab_builder=builder, tab_number=tab_number, admin_only=admin_only,
                                              tab_id=tab_id, label=label, badge=badge, badge_status=badge_status, resolved_url=content))
        else:
            builder.tabs.append(TabBuilderTab(tab_builder=builder, tab_number=tab_number, admin_only=admin_only,
                                              tab_id=tab_id, label=label, badge=badge, badge_status=badge_status, content=content))
        return ""
