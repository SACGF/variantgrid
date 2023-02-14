import html
import re
import uuid
from html import escape
from typing import Optional, Set, Dict

from bs4 import BeautifulSoup
from django.utils.safestring import mark_safe, SafeString


def html_id_safe(text: str) -> str:
    """
    Makes a string that can be made into an HTML id and referenced easily
    """
    if not text:
        return str(uuid.uuid4())
    text = re.sub("[^0-9A-Za-z]", "-", text)
    text = re.sub("_{2,}", "-", text)
    if not text[0].isalpha():
        text = "x" + text
    return text


def html_link(url: str, title: str) -> SafeString:
    if not url:
        return mark_safe(title)
    return mark_safe(f"<a href='{url}'>{html.escape(title)}</a>")


# note this tags expected in a single line of text
# don't catch too many tags in case you get some false positives
EXPECTED_HTML_TAGS_SINGLE_LINE = {'div', 'b', 'i', 'u', 'strong', 'em', 'sup', 'sub'}


def cautious_attempt_html_to_text(text: str, whitelist: Set[str] = None) -> str:
    """
    Given some text, and an expected whitelist of tags, will convert the text from possible HTML content to plain text.
    Converts things like &#039; to ' and will strip out expected tags, if a tag is found that's not in the whitelist
    the text will be returned untouched. This is to avoid treating text like "Patient had ouchies <painful>" as HTML
    :param text: Text that may contain HTML elements
    :param whitelist: HTML tags that we expect and want stripped out
    :return: Text stripped of HTML (if only tags present were whitelisted)
    """
    if whitelist is None:
        whitelist = EXPECTED_HTML_TAGS_SINGLE_LINE

    if not text:
        return text
    bs = BeautifulSoup(text, features="html.parser")
    for tag in bs.find_all():
        if tag.name.lower() not in whitelist:
            return text
    return bs.get_text()


def html_to_text(html: str, preserve_lines: bool = False) -> Optional[str]:
    if not html:
        return None
    bs = BeautifulSoup(f'<body>{html}</body>', features="html.parser")

    if not preserve_lines:
        return bs.get_text()
    else:

        def replace_with_newlines(element):
            text = ''
            for elem in element.recursiveChildGenerator():
                if isinstance(elem, str):
                    text += elem.strip() + " "
                elif elem.name == 'br':
                    text += '\n'
            return text.strip()

        def get_plain_text(soup):
            lines = soup.find("body")
            return replace_with_newlines(lines)

        return get_plain_text(bs).strip()


class IconWithTooltip:

    ERROR_ICON = 'fas fa-exclamation-circle text-danger'
    WARNING_ICON = 'fas fa-exclamation-triangle text-warning'
    HOURGLASS_START = 'fa-solid fa-hourglass-start'
    HOURGLASS_MID = 'fa-regular fa-hourglass-half'
    HOURGLASS_END = 'fa-solid fa-hourglass-end'

    def __init__(self, icon: str, tooltip: Optional[str] = None):
        self.icon = icon
        self.tooltip = tooltip

    def __str__(self):
        title = ""
        if tooltip := self.tooltip:
            title = f'title="{escape(tooltip)}"'
        return SafeString(f'<i class="{escape(self.icon)}" {title}></i>')

    def as_json(self) -> Dict:
        return {
            "icon": self.icon,
            "tooltip": self.tooltip
        }
