import html
import re
import uuid
from typing import Optional

import nh3
from bs4 import BeautifulSoup
from django.utils.html import format_html
from django.utils.safestring import SafeString


def html_id_safe(text: str) -> str:
    """
    Makes a string that can be made into an HTML id and referenced easily
    """
    if not text:
        return str(uuid.uuid4())
    text = str(text)
    text = re.sub("[^0-9A-Za-z]", "-", text)
    text = re.sub("_{2,}", "-", text)
    if not text[0].isalpha():
        text = "x" + text
    return text


# note this tags expected in a single line of text
# don't catch too many tags in case you get some false positives
EXPECTED_HTML_TAGS_SINGLE_LINE = {'div', 'b', 'i', 'u', 'strong', 'em', 'sup', 'sub'}


def cautious_attempt_html_to_text(text: str, whitelist: Optional[set[str]] = None) -> str:
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


# Tags/attributes/protocols allowed when displaying untrusted text that may contain HTML
# (e.g. user-to-user inbox messages). Keeps basic formatting and links, escapes everything else.
SANITIZE_HTML_ALLOWED_TAGS = {'a', 'b', 'i', 'u', 'strong', 'em', 'sup', 'sub', 'br', 'p',
                              'ul', 'ol', 'li', 'span', 'div'}
SANITIZE_HTML_ALLOWED_ATTRIBUTES = {'a': {'href', 'title'}}
SANITIZE_HTML_ALLOWED_PROTOCOLS = {'http', 'https', 'mailto'}

# A comment, a start/end tag (quoted attribute values may contain '<' or '>'), or a stray '<' / '>'
_ANGLE_BRACKET_RE = re.compile(r"<!--.*?-->|</?([a-zA-Z][a-zA-Z0-9-]*)(?:[^<>\"']|\"[^\"]*\"|'[^']*')*>|[<>]", re.DOTALL)


def _escape_disallowed_tags(text: str) -> str:
    """ nh3 drops a disallowed tag, so text like "c.123<DEL>" would lose "<DEL>". Escape anything that isn't
        an allowed tag (or a comment, which nh3 strips) so it displays as text. """
    def escape_unless_allowed(m: re.Match) -> str:
        token = m.group(0)
        if token.startswith("<!--") or (m.group(1) and m.group(1).lower() in SANITIZE_HTML_ALLOWED_TAGS):
            return token
        return html.escape(token, quote=False)

    return _ANGLE_BRACKET_RE.sub(escape_unless_allowed, text)


def sanitize_html(text: str) -> SafeString:
    """
    Sanitize untrusted text that may contain HTML so it is safe to render unescaped.
    A small whitelist of formatting tags and links is preserved; any other tags/attributes
    (e.g. <script>, onerror=, javascript: hrefs) are escaped or stripped.
    :param text: Untrusted text that may contain HTML
    :return: A SafeString containing only whitelisted HTML
    """
    if not text:
        return SafeString("")
    cleaned = nh3.clean(_escape_disallowed_tags(text),
                        tags=SANITIZE_HTML_ALLOWED_TAGS,
                        attributes=SANITIZE_HTML_ALLOWED_ATTRIBUTES,
                        url_schemes=SANITIZE_HTML_ALLOWED_PROTOCOLS,
                        link_rel=None)
    return SafeString(cleaned)


def html_to_text(html_str: str, preserve_lines: bool = False) -> Optional[str]:
    if not html_str:
        return None
    bs = BeautifulSoup(f'<body>{html_str}</body>', features="html.parser")

    if not preserve_lines:
        return bs.get_text()
    else:

        def replace_with_newlines(element):
            text = ''
            for elem in element.recursiveChildGenerator():
                if isinstance(elem, str):
                    # why did we used to add a space after a text snippet?
                    text += elem  # .strip() + " "
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
        if self.tooltip:
            return format_html('<i class="{}" title="{}"></i>', self.icon, self.tooltip)
        return format_html('<i class="{}"></i>', self.icon)

    def as_json(self) -> dict:
        return {
            "icon": self.icon,
            "tooltip": self.tooltip
        }
