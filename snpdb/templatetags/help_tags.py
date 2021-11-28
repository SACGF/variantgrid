import os

from django.conf import settings
from django.contrib.staticfiles import finders
from django.template import Library
from django.utils.safestring import SafeText

register = Library()


@register.inclusion_tag("help_tags/help_tags.html", takes_context=True)
def page_help(context, user=None, page_id: str = None, title=None, show_title=True):
    """
    Displays a heading and page help
    :param context:
    :param user: deprecated, if provided old behaviour is used (no h3, don't affix Help to the end of the help)
    :param page_id: Used to determine the help page
    :param title: The title for the heading and the help drop down
    :param show_title: Should we show a heading at all
    """
    suffix_help = True
    if not user:
        user = context.request.user
    else:
        show_title = False
        suffix_help = False

    if not title or title == SafeText(' '):
        help_page_title = 'Page Help'
    else:
        help_page_title = title
        if suffix_help:
            help_page_title += ' Help'

    help_url = settings.HELP_URL
    page_help_html = None
    show_page_help = False
    page_help_path = f"page_help/{page_id}.html"
    page_help_filename = finders.find(page_help_path)
    if page_help_filename and os.path.exists(page_help_filename):
        with open(page_help_filename) as ph:
            page_help_html = ph.read()

    return {
        'show': show_page_help,
        'page_id': page_id.replace('/', '-').replace(' ', '-'),
        'page_title': title if show_title else None,
        'help_page_title': help_page_title,
        "page_help_html": page_help_html,
        "help_url": help_url}
