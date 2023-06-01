from dataclasses import dataclass
from typing import List, Tuple

from django.db.models import Model
from django.forms import ChoiceField
import copy

from django.forms.widgets import ChoiceWidget
from more_itertools import first


@dataclass(frozen=True)
class OptionData:
    name: str
    label: str
    value: str
    selected: bool
    index: int


class RadioOtherWidget(ChoiceWidget):
    template_name = 'uicore/widgets/radio_other_widget.html'
    allow_multiple_selected = False

    def __init__(self, attrs=None, choices=()):
        super().__init__(attrs)
        self.choices = list(choices)

    def get_context(self, name, value, attrs):
        context = super().get_context(name, value, attrs)

        value = self.format_value(value)
        if value is None:
            value = []
        options, other_option = self.other_options(name=name, value=value, attrs=attrs)

        context["widget"]["options"] = options
        context["widget"]["other"] = other_option
        return context

    def other_options(self, name: str, value: List[str], attrs) -> Tuple[List[OptionData], OptionData]:
        option_datas: List[OptionData] = []
        value_set = set(value)
        if '' in value_set:
            value_set.remove('')

        has_selected = False
        index = 0
        for index, (option_value, option_label) in enumerate(self.choices):
            if option_value is None:
                option_value = ""

            matches_value = False
            if option_value in value_set:
                value_set.remove(option_value)
                matches_value = self.allow_multiple_selected or not has_selected
                has_selected = True

            option_datas.append(OptionData(
                name=name,
                label=option_label,
                value=option_value,
                selected=matches_value,
                index=index
            ))

        other_option = OptionData(
            name=name + "-other",
            label="Other",
            value=first(value_set, ""),
            selected=bool(value_set),
            index=index+1
        )
        return option_datas, other_option

    def __deepcopy__(self, memo):
        obj = copy.copy(self)
        obj.attrs = self.attrs.copy()
        obj.choices = copy.copy(self.choices)
        memo[id(self)] = obj
        return obj

    def value_from_datadict(self, data, files, name):
        if self.allow_multiple_selected:
            values = set()
            for choice in self.choices:
                choice_key = choice[0]
                if data.get(f"{name}-{choice_key}") == "on":
                    values.add(choice[0])
            if data.get(f"{name}-other") == "on":
                if other_value := data.get(f"{name}-other-value"):
                    values.add(other_value)
            return list(sorted(values))
        else:
            selected_value = super().value_from_datadict(data, files, name)
            if selected_value == "other":
                selected_value = data.get(f"{name}-other")
            if not selected_value:
                return None
            return selected_value

    def format_value(self, value):
        """Return selected values as a list."""
        if value is None and self.allow_multiple_selected:
            return []
        elif not isinstance(value, (tuple, list, set)):
            value = [value]

        def format_single_value(v):
            if isinstance(v, Model):
                return str(v.pk)
            else:
                return str(v)

        return [format_single_value(v) if v is not None else "" for v in value]


class CheckboxOtherWidget(RadioOtherWidget):
    template_name = 'uicore/widgets/checkbox_other_widget.html'
    allow_multiple_selected = True


class ChoiceFieldWithOther(ChoiceField):
    widget = RadioOtherWidget
    other_enabled = True

    def widget_attrs(self, widget):
        return {"other_enabled": self.other_enabled}

    def valid_value(self, value):
        return True

    def to_python(self, value):
        # don't want to turn an array into a string representation of an array
        return value


class MultiChoiceFieldWithOther(ChoiceField):
    widget = CheckboxOtherWidget
    other_enabled = True

    def widget_attrs(self, widget):
        return {"other_enabled": self.other_enabled}

    def valid_value(self, value):
        return True

    def to_python(self, value):
        # don't want to turn an array into a string representation of an array
        return value
