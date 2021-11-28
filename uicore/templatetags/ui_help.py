import os

from django import template
from django.conf import settings
from django.contrib.staticfiles import finders
from django.template import Library, loader
from django.template.base import FilterExpression
from django.utils.safestring import SafeText

from library.log_utils import report_message
from library.utils import html_id_safe
from uicore.templatetags.ui_utils import parse_tag, TagUtils

register = Library()


@register.inclusion_tag("uicore/tags/help.html")
def page_help(page_id: str = None, title=None, show_title=True, header_tag="h3"):
    """
    Best to use page_help_embedded now so help can be dynamic.

    Displays a heading and page help loaded from a file.
    :param page_id: Used to determine the help page
    :param title: The title for the heading and the help drop down
    :param show_title: Should we show a heading at all - deprecated, generally always show title (use page_help instead
    of header if page help is available).
    :param header: Default header if there's no attached help
    """
    suffix_help = True

    if not title or title == SafeText(' '):
        title = 'Page Help'

    help_url = settings.HELP_URL
    page_help_html = None
    page_help_path = f"page_help/{page_id}.html"
    page_help_filename = finders.find(page_help_path)
    file_exists = False
    if page_help_filename and os.path.exists(page_help_filename):
        with open(page_help_filename) as ph:
            file_exists = True
            page_help_html = ph.read()
    if not file_exists:
        report_message(f"Could not find help for {page_help_path}")

    return {
        'page_id': html_id_safe(page_id),
        'page_title': title if show_title else None,
        'help_page_title': title,
        "page_help_html": page_help_html,
        "header_tag": header_tag,
        "help_url": help_url}


@register.tag(name='page_help_embedded')
def page_help_embedded(parser, token):
    tag_name, args, kwargs = parse_tag(token, parser)
    nodelist = parser.parse(('end_page_help_embedded',))
    parser.delete_first_token()
    if title := kwargs.get("title") or args[0]:
        return PageHelpContent(nodelist, title=title)
    else:
        raise ValueError("page_help_embedded must have attribute 'title'")


class PageHelpContent(template.Node):

    def __init__(self,
                 nodelist,
                 title: FilterExpression):
        self.nodelist = nodelist
        self.title = title

    def render(self, context):

        title = TagUtils.value_str(context, self.title)
        page_id = html_id_safe(title)
        content = self.nodelist.render(context).strip()

        return loader.render_to_string("uicore/tags/help.html", context={
            'page_id': page_id,
            'page_title': title,
            'page_help_html': content,
            'header_tag': 'h4',  # can revert back to just a tag if there's no help content
            'help_url': settings.HELP_URL
        })
