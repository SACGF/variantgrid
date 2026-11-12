from django import forms
from django.test import SimpleTestCase

from uicore.widgets.optional_range_widget import OptionalRangeWidget


class ZoomForm(forms.Form):
    zoom = forms.IntegerField(required=False, widget=OptionalRangeWidget())


class OptionalRangeWidgetTest(SimpleTestCase):
    """ The slider only counts when its checkbox is ticked - unticked is null whatever the slider posts """

    def _clean_zoom(self, data):
        form = ZoomForm(data)
        self.assertTrue(form.is_valid(), form.errors)
        return form.cleaned_data["zoom"]

    def test_unticked_is_null(self):
        self.assertIsNone(self._clean_zoom({"zoom": "120"}))

    def test_ticked_is_slider_value(self):
        self.assertEqual(self._clean_zoom({"zoom-set": "on", "zoom": "120"}), 120)
