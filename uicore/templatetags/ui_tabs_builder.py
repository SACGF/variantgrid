from dataclasses import dataclass
from typing import List, Any, Optional
from django import template
from django.template.base import FilterExpression

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
    url: Optional[str] = None
    resolved_url: Optional[str] = None
    param: Any = None
    content: str = ''

    @property
    def active(self):
        return self.tab_builder.active_tab == self.tab_number


class TabBuilder:

    def __init__(self):
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
        builder = TabBuilder()
        context[tab_key] = builder

    tab_number = len(builder.tabs)
    if active:
        builder.active_tab = tab_number
    builder.tabs.append(TabBuilderTab(tab_builder=builder, tab_number=tab_number,
                                      tab_id=url, label=label, badge=badge, badge_status=badge_status, url=url, param=param))
    return ""


@register.simple_tag(takes_context=True)
def ui_register_tabs(context, tab_set: str):
    """
    Used because the creation of he tabs in subcontext e.g. for loops, will dissapear when poping to a higher context
    """
    tab_key = f"ui-tab-{tab_set}"
    context[tab_key] = TabBuilder()
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
        badge=kwargs.get('badge'),
        badge_status=kwargs.get('badge_status')
    )


class LocalTabContent(template.Node):
    def __init__(self,
                 nodelist,
                 tab_set: FilterExpression,
                 label: FilterExpression,
                 badge: FilterExpression,
                 badge_status: FilterExpression):
        self.nodelist = nodelist
        self.tab_set = tab_set
        self.label = label
        self.badge = badge
        self.badge_status = badge_status

    def render(self, context):
        tab_set = TagUtils.value_str(context, self.tab_set)
        label = TagUtils.value_str(context, self.label)
        badge = TagUtils.value_int(context, self.badge)
        badge_status = TagUtils.value_str(context, self.badge_status)
        if not tab_set:
            raise ValueError("UI Tab requires a value for 'tab_set'")
        if not label:
            raise ValueError("UI Tab requires a value for 'label'")

        # FIXME just replace everything that's not A-Z 0-9 and make sure dont start with a number
        tab_id = 't' + label.replace(' ', '-').replace('/', '_').replace('!', 'not').replace('=', '_')

        tab_key = f"ui-tab-{tab_set}"
        builder: TabBuilder = context.get(tab_key)
        if not builder:
            builder = TabBuilder()
            context[tab_key] = builder

        tab_number = len(builder.tabs)
        content: str = self.nodelist.render(context)
        if content.startswith('/'):
            builder.tabs.append(TabBuilderTab(tab_builder=builder, tab_number=tab_number,
                                              tab_id=tab_id, label=label, badge=badge, badge_status=badge_status, resolved_url=content))
        else:
            builder.tabs.append(TabBuilderTab(tab_builder=builder, tab_number=tab_number,
                                              tab_id=tab_id, label=label, badge=badge, badge_status=badge_status, content=content))
        return ""
