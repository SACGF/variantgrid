from typing import Iterable, Any

from django.core.exceptions import ValidationError

from snpdb.models import Lab
from uicore.widgets.radio_other_widget import MultiChoiceFieldWithOther, CheckboxOtherWidget


class MultiChoiceLabWidget(CheckboxOtherWidget):
    pass


class MultiChoiceLabField(MultiChoiceFieldWithOther):
    other_enabled = False

    def __init__(self, labs: Iterable[Lab], initial: Any = None):
        self.labs = {
            lab.pk: lab for lab in sorted(labs)
        }
        choices = [(str(lab.pk), str(lab)) for lab in labs]
        super().__init__(choices=choices, initial=initial)

    def to_python(self, value):
        if isinstance(value, (set, list, tuple)):
            result = []
            for v in value:
                try:
                    result.append(self.labs[int(v)])
                except (KeyError, ValueError):
                    raise ValidationError("Invalid lab selection.")
            return result
        return []
