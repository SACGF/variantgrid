from typing import Optional

from django import template
from django.utils.safestring import SafeString

from upload.file_type_icons import file_type_icon_html

register = template.Library()


@register.simple_tag
def file_type_icon(file_type: Optional[str]) -> SafeString:
    """ The icon for an UploadedFileTypes code - see upload/file_type_icons.py """
    return file_type_icon_html(file_type)
