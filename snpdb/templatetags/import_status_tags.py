"""
`{% import_status_icon status %}` - an ImportStatus as its icon and label, the same icon the data
grids draw through TableFormat.importStatus (@see ImportStatus.icon_css)
"""
from django.template import Library

from snpdb.enums import ImportStatus

register = Library()


@register.inclusion_tag("snpdb/tags/import_status_icon.html")
def import_status_icon(status: str, show_label: bool = True):
    import_status = ImportStatus(status)
    return {"icon_css": import_status.icon_css, "label": import_status.label, "show_label": show_label}
