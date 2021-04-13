import json
import uuid
from typing import Optional, Any

from django import template
from django.template.base import FilterExpression, kwarg_re
import re

from django.utils.safestring import SafeString

register = template.Library()


@register.simple_tag(takes_context=True)
def update_django_messages(context):
    """
    Use when you've loaded messages in an ajax tab, and need to update the page messages with save etc messages
    """
    message_json = []
    if messages := context.get("messages"):
        for message in messages:
            tags = message.tags
            text = str(message)
            message_json.append({"tags": tags, "text": text})
    message_json_string = json.dumps(message_json)
    return SafeString(f"update_django_messages({message_json_string});")


# Taken from https://www.caktusgroup.com/blog/2017/05/01/building-custom-block-template-tag/
# how is this not built in?!
def parse_tag(token, parser):
    """
    Generic template tag parser.
    Returns a three-tuple: (tag_name, args, kwargs)
    tag_name is a string, the name of the tag.
    args is a list of FilterExpressions, from all the arguments that didn't look like kwargs,
    in the order they occurred, including any that were mingled amongst kwargs.
    kwargs is a dictionary mapping kwarg names to FilterExpressions, for all the arguments that
    looked like kwargs, including any that were mingled amongst args.
    (At rendering time, a FilterExpression f can be evaluated by calling f.resolve(context).)
    """
    # Split the tag content into words, respecting quoted strings.
    bits = token.split_contents()

    # Pull out the tag name.
    tag_name = bits.pop(0)

    # Parse the rest of the args, and build FilterExpressions from them so that
    # we can evaluate them later.
    args = []
    kwargs = {}
    for bit in bits:
        # Is this a kwarg or an arg?
        match = kwarg_re.match(bit)
        kwarg_format = match and match.group(1)
        if kwarg_format:
            key, value = match.groups()
            kwargs[key] = FilterExpression(value, parser)
        else:
            args.append(FilterExpression(bit, parser))

    return (tag_name, args, kwargs)


@register.tag(name='install-instructions')
def render_install_instructions(parser, token):
    tag_name, args, kwargs = parse_tag(token, parser)
    nodelist = parser.parse(('endinstall-instructions',))
    parser.delete_first_token()
    label = kwargs.get('label')
    installed = kwargs.get('installed')
    if not label:
        if len(args) > 0:
            label = args[0]
        else:
            raise ValueError("install-instructions required requires 'label' argument")

    return InstallInstructionsTag(nodelist, label=label, installed=installed)


class InstallInstructionsTag(template.Node):

    def __init__(self, nodelist, installed: bool, label: FilterExpression = None):
        self.nodelist = nodelist
        self.installed = installed
        self.label = label

    def render(self, context):
        if not context.request.user.is_superuser:
            return ""

        label_str = TagUtils.value_str(context, self.label)
        if not label_str:
            label_str = ""
            id_safe = str(uuid.uuid4()).replace("-", "_") + "_instructions"
        else:
            id_safe = re.sub(r"\W", "_", label_str).lower()
        css_classes = ["install-instructions"]
        if self.installed.resolve(context):
            css_classes.append("collapse")
        else:
            css_classes.append("not-installed")

        return f"""
        <div>
        <a class='toggle-link install-instructions-toggle' data-toggle='collapse' href='#{id_safe}'>{label_str} Install/Update Instructions</a>
        <div class='{' '.join(css_classes)}' id='{id_safe}'>
        {self.nodelist.render(context)}
        </div>
        </div>
        """


@register.tag(name='labelled')
def render_labelled(parser, token):
    tag_name, args, kwargs = parse_tag(token, parser)
    nodelist = parser.parse(('endlabelled',))
    parser.delete_first_token()
    return LabelledValueTag(nodelist,
                            id_prefix=kwargs.get('id_prefix'),
                            value_id=kwargs.get('id'),
                            label=kwargs.get('label'),
                            hint=kwargs.get('hint'),
                            label_css=kwargs.get('label_css'),
                            value_css=kwargs.get('value_css'),
                            row_css=kwargs.get('row_css'),
                            shorten_label=kwargs.get('shorten_label')
                            )


class LabelledValueTag(template.Node):
    def __init__(self, nodelist,
                 id_prefix: FilterExpression,
                 value_id: FilterExpression,
                 label: FilterExpression = None,
                 hint: FilterExpression = None,
                 label_css: FilterExpression = None,
                 value_css: FilterExpression = None,
                 row_css: FilterExpression = None,
                 shorten_label: FilterExpression = None):
        self.id_prefix = id_prefix
        self.nodelist = nodelist
        self.value_id = value_id
        self.label = label
        self.hint = hint
        self.label_css = label_css
        self.value_css = value_css
        self.row_css = row_css
        self.shorten_label = shorten_label

    id_regex = re.compile(r"id=[\"|'](.*?)[\"|']")
    big_zero = re.compile(r"^0([.]0+)?$")

    def render(self, context):
        prefix_id = TagUtils.value_str(context, self.id_prefix)
        value_id = TagUtils.value_str(context, self.value_id)
        complete_id = value_id
        if prefix_id and value_id:
            complete_id = f"{prefix_id}-{value_id}"
        label = TagUtils.value_str(context, self.label, (value_id.replace('_', ' ') if value_id else ""))

        popover = None

        if TagUtils.value_bool(context, self.shorten_label):
            first_fullstop = label.find('. ')
            if first_fullstop != -1:
                first_fullstop += 1

            if first_fullstop == -1:
                first_fullstop = label.find('</a>')
                if first_fullstop != -1:
                    first_fullstop += 4

            if first_fullstop != -1:
                popover = label
                label = label[0:first_fullstop]

        hint = TagUtils.value_str(context, self.hint)
        label_css = ""
        value_css = ""
        row_css = "form-group row mb-4 mb-md-3"

        if hint == "tiny":
            label_css = "col-12 col-md-6 text-md-right"
            value_css = "col-12 col-md-6 text-left align-self-center text-break"
        elif hint == "chunky":
            label_css = "col-12"
            value_css = "col-12 text-break"
        elif hint == "inline":
            label_css = "m-2 align-self-center"
            value_css = "m-2"
        else:
            label_css = "col-12 col-md-3 text-md-right align-self-center"
            value_css = "col-12 col-md-9 text-left text-break"

        label_css_extra = TagUtils.value_str(context, self.label_css, '')
        label_css = f"{label_css} {label_css_extra}".strip()

        value_css_extra = TagUtils.value_str(context, self.value_css, '')
        value_css = f"{value_css} {value_css_extra}".strip()

        row_css_extra = TagUtils.value_str(context, self.row_css, "")
        row_css = f"{row_css} {row_css_extra}".strip()

        output = self.nodelist.render(context)
        give_div_id = complete_id
        if not complete_id and not '<label' in output:
            if found_id := LabelledValueTag.id_regex.search(output):
                complete_id = found_id.group(1)
                give_div_id = False

        div_id = ""
        for_id = ""
        if give_div_id:
            div_id = f"id=\"{complete_id}\""
            for_id = f"for=\"{complete_id}\""

        if output in ("", "None"):
            output = "<span class=\"no-value\">-</span>"
        elif LabelledValueTag.big_zero.match(output):
            output = f"<span class=\"zero-value\">{output}</span>"

        label_tag = f'<label {for_id} class="{label_css}">{label}</label>'
        if popover:
            popover = popover.replace('"', '&quot;')
            label_tag = f'<label {for_id} class="{label_css} hover-detail" data-toggle="popover" data-content="{popover}">{label}</label>'

        content = f"""{label_tag}<div {div_id} class="{value_css}">{output}</div>"""

        if hint == "inline":
            return content
        return f'<div class="{row_css}">{content}</div>'


@register.filter()
def severity_icon(severity: str) -> str:
    if not severity:
        severity = 'I'
    severity = severity.upper()
    if severity.startswith('C'):  # critical
        return SafeString('<i class="fas fa-bomb text-danger"></i>')
    if severity.startswith('E'):  # error
        return SafeString('<i class="fas fa-exclamation-circle text-danger"></i>')
    if severity.startswith('W'):  # warning
        return SafeString('<i class="fas fa-exclamation-triangle text-warning"></i>')
    if severity.startswith('I'):  # info
        return SafeString('<i class="fas fa-info-circle text-info"></i>')
    # debug
    return SafeString('<i class="fas fa-question-circle text-secondary"></i>')


@register.filter()
def danger_badge(count: Optional[int]) -> str:
    """
    If count is None, does nothing
    If count is 0, shows a success badge with 0 in it
    If count is anything else, shown in a danger badge
    """
    if count is None:
        return ""
    if count == 0:
        return SafeString(' <span class="d-inline-block ml-1 badge badge-success">0</span>')
    return SafeString(f' <span class="d-inline-block ml-1 badge badge-danger">{count}</span>')


@register.filter()
def checked(test: bool) -> str:
    if test:
        return SafeString('checked="checked"')
    return ''


class TagUtils:

    @staticmethod
    def value(context, value: Any, default: Any = None) -> Optional[Any]:
        if isinstance(value, FilterExpression):
            resolved = value.resolve(context)
        else:
            resolved = value
        if resolved is None:
            resolved = default
        return resolved

    @staticmethod
    def value_str(context, value: Any, default: Any = None) -> Optional[str]:
        val = TagUtils.value(context, value, default)
        if val is not None:
            return str(val)
        return None

    @staticmethod
    def value_int(context, value: Any, default: Any = None) -> Optional[int]:
        val = TagUtils.value(context, value, default)
        if val is not None and val != '':
            try:
                return int(val)
            except ValueError:
                pass
        return None

    @staticmethod
    def value_bool(context, value: Any, default: Any = None) -> Optional[bool]:
        val = TagUtils.value(context, value, default)
        if val is not None:
            return bool(val)
        return None
