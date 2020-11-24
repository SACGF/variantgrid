from dataclasses import dataclass
from typing import List, Any, Optional
from django import template
from django.template.base import FilterExpression

from uicore.templatetags.ui_utils import parse_tag, TagUtils

register = template.Library()

@dataclass
class TabBuilderTab:
    tab_id: str
    label: str
    url: Optional[str] = None
    resolved_url: Optional[str] = None
    param: Any = None
    content: str = ''


class TabBuilder:

    def __init__(self):
        self.tabs: List[TabBuilderTab] = []

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
def ui_register_tab(context, tab_set: str, label: str, url: str = None, param: Any = None, url_check=False):
    if url_check:
        if url not in context["url_name_visible"]:
            return ""

    tab_key = f"ui-tab-{tab_set}"
    builder: TabBuilder = context.get(tab_key)
    if not builder:
        builder = TabBuilder()
        context[tab_key] = builder

    builder.tabs.append(TabBuilderTab(tab_id=url, label=label, url=url, param=param))
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
    return LocalTabContent(nodelist,
                            tab_set=kwargs.get('tab_set'),
                            label=kwargs.get('label')
    )


class LocalTabContent(template.Node):
    def __init__(self,
                 nodelist,
                 tab_set: FilterExpression,
                 label: FilterExpression):
        self.nodelist = nodelist
        self.tab_set = tab_set
        self.label = label

    def render(self, context):
        tab_set = TagUtils.value_str(context, self.tab_set)
        label = TagUtils.value_str(context, self.label)
        if not tab_set:
            raise ValueError("UI Tab requires a value for 'tab_set'")
        if not label:
            raise ValueError("UI Tab requires a value for 'label'")

        tab_id = label.replace(' ', '-').replace('/', '_')

        tab_key = f"ui-tab-{tab_set}"
        builder: TabBuilder = context.get(tab_key)
        if not builder:
            builder = TabBuilder()
            context[tab_key] = builder

        content: str = self.nodelist.render(context)
        if content.startswith('/'):
            builder.tabs.append(TabBuilderTab(tab_id=tab_id, label=label, resolved_url=content))
        else:
            builder.tabs.append(TabBuilderTab(tab_id=tab_id, label=label, content=content))
        return ""
