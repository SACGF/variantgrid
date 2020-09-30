from django.conf import settings
from django.contrib.staticfiles import finders
from django import template
from django.template import Library, loader
from django.template.base import FilterExpression
from django.utils.safestring import SafeText
import os

from uicore.templatetags.ui_utils import parse_tag, TagUtils

register = Library()


@register.inclusion_tag("uicore/tags/help.html")
def page_help(page_id: str = None, title=None, show_title=True):
    """
    Displays a heading and page help
    :param page_id: Used to determine the help page
    :param title: The title for the heading and the help drop down
    :param show_title: Should we show a heading at all
    """
    suffix_help = True

    if not title or title == SafeText(' '):
        help_page_title = 'Page Help'
    else:
        help_page_title = title
        if suffix_help:
            help_page_title += ' Help'

    help_url = settings.HELP_URL
    page_help_html = None
    page_help_path = f"page_help/{page_id}.html"
    page_help_filename = finders.find(page_help_path)
    if page_help_filename and os.path.exists(page_help_filename):
        with open(page_help_filename) as ph:
            page_help_html = ph.read()

    return {
        'page_id': page_id.replace('/', '-').replace(' ', '-'),
        'page_title': title if show_title else None,
        'help_page_title': help_page_title,
        "page_help_html": page_help_html,
        "help_url": help_url}


@register.inclusion_tag("uicore/tags/help.html")
def page_help_manual(page_label: str, content: str):

    page_id = page_label.replace(' ', '-').replace('/', '_')
    help_url = settings.HELP_URL

    return {
        'page_id': page_id,
        'page_title': None,
        'help_page_title': page_label + ' Help',
        'page_help_html': content,
        'help_url': help_url
    }


@register.tag(name='page_help_embedded')
def page_help_embedded(parser, token):
    tag_name, args, kwargs = parse_tag(token, parser)
    nodelist = parser.parse(('end_page_help_embedded',))
    parser.delete_first_token()
    return PageHelpContent(nodelist,
                            title=kwargs.get('title')
    )


class PageHelpContent(template.Node):

    def __init__(self,
                 nodelist,
                 title: FilterExpression):
        self.nodelist = nodelist
        self.title = title

    def render(self, context):
        title = TagUtils.value_str(context, self.title)
        page_id = title.replace(' ', '-').replace('/', '_')
        content = self.nodelist.render(context)

        return loader.render_to_string("uicore/tags/help.html", context={
            'page_id': page_id,
            'page_title': None,
            'help_page_title': title + ' Help',
            'page_help_html': content,
            'help_url': settings.HELP_URL
        })
