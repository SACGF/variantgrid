"""
OptionalRangeWidget: a slider for a nullable number, behind a checkbox. Unticked hides the slider and the
field saves as null - eg a settings override left blank so the lab / organisation value applies. Behaviour
is uicore/optional_range_widget.js, which the widget's Media loads.
"""
from django.forms.widgets import Input


class OptionalRangeWidget(Input):
    input_type = "range"
    template_name = "uicore/widgets/optional_range_widget.html"

    class Media:
        js = ["js/uicore/optional_range_widget.js"]

    def __init__(self, attrs=None, initial=None, value_suffix=""):
        """ initial: where the slider sits when first ticked, value_suffix: shown after the value (eg "%") """
        super().__init__(attrs)
        self.initial = initial
        self.value_suffix = value_suffix

    def get_context(self, name, value, attrs):
        context = super().get_context(name, value, attrs)
        widget = context["widget"]
        is_set = widget["value"] is not None
        widget["is_set"] = is_set
        widget["read_only"] = bool(widget["attrs"].get("disabled"))
        widget["value_suffix"] = self.value_suffix
        widget["attrs"]["class"] = "custom-range w-auto optional-range-input"
        if not is_set:
            widget["attrs"].update({"disabled": True, "hidden": True})
            if self.initial is not None:
                widget["value"] = str(self.initial)
        return context

    def value_from_datadict(self, data, files, name):
        if data.get(f"{name}-set"):
            return data.get(name)
        return None

    def value_omitted_from_data(self, data, files, name):
        # An unticked box is a value (null), so a ModelForm must not keep the old one
        return False
