from dataclasses import dataclass
from typing import Optional, Dict

from django.core.exceptions import ValidationError
from django.forms import Field, Widget

from django.db.models import TextChoices


class DifferenceResolution(TextChoices):
    RESOLVED = "R", "Resolved"
    UNRESOLVED = "U", "Still Outstanding"


@dataclass
class DescribeDifference:
    selected: bool = False
    details: Optional[str] = None
    resolution: Optional[DifferenceResolution] = None

    @staticmethod
    def from_json(data: Dict) -> 'DescribeDifference':
        if not data:
            return DescribeDifference()
        return DescribeDifference(
            selected=bool(data.get("selected")),
            details=data.get("details") or "",
            resolution=data.get("resolution")
        )

    def __bool__(self):
        return self.selected

    def as_json(self):
        return {
            "selected": self.selected,
            "details": self.details,
            "resolution": self.resolution.value
        }


class DescribeDifferenceWidget(Widget):
    template_name = 'uicore/widgets/describe_difference_widget.html'

    def get_context(self, name, value, attrs):
        context = super().get_context(name, self.format_value(value), attrs)
        return context

    def value_from_datadict(self, data, files, name):
        is_selected = data.get(name) == "on"
        if is_selected:
            details = data.get(name + "-details")
            resolution = DifferenceResolution(data.get(name + "-resolution"))
            return DescribeDifference(
                selected=True,
                details=details,
                resolution=resolution
            )
        else:
            return DescribeDifference()

    def format_value(self, value):
        if value is None:
            value = DescribeDifference()
        return value


class DescribeDifferenceField(Field):
    widget = DescribeDifferenceWidget

    def __init__(self,  category: str = None, **kwargs):
        super().__init__(**kwargs)
        self.category = category

    def widget_attrs(self, widget):
        # hack to get this information to the widget
        # normally this is rendered elsewhere but the label is built into the widget rendering here
        return {
            "label": self.label,
            "help_text": self.help_text
        }

    def validate(self, value):
        super().validate(value)
        if isinstance(value, DescribeDifference):
            if value.selected and (not value.details or not value.resolution):
                raise ValidationError("If selected, details and a resolution must be populated")