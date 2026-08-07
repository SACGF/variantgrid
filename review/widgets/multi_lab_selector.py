from collections.abc import Iterable
from typing import Any

from snpdb.models import Lab
from uicore.widgets.radio_other_widget import CheckboxOtherWidget, MultiChoiceFieldWithOther


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
            return [self.labs[int(v)] for v in value]
        return []
