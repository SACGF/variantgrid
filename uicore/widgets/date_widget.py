from datetime import date

from django.forms import DateInput


class NativeDateInput(DateInput):
    """ <input type="date"> - the browser supplies the calendar.

        The value rendered and the value submitted are both ISO yyyy-mm-dd, which heads both
        DATE_INPUT_FORMATS and (via Django's fallback) DATETIME_INPUT_FORMATS, so date and
        datetime fields alike round-trip through form parsing. """

    input_type = "date"
    EARLIEST_YEARS_AGO = 120

    def __init__(self, attrs=None, allow_future: bool = False):
        self.allow_future = allow_future
        super().__init__(attrs=attrs, format="%Y-%m-%d")

    def get_context(self, name, value, attrs):
        context = super().get_context(name, value, attrs)
        widget_attrs = context["widget"]["attrs"]
        today = date.today()
        widget_attrs.setdefault("min", f"{today.year - self.EARLIEST_YEARS_AGO}-01-01")
        if not self.allow_future:
            widget_attrs.setdefault("max", today.isoformat())
        return context
