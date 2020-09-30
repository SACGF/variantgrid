from crispy_forms.helper import FormHelper


def form_helper_horizontal() -> FormHelper:
    helper = FormHelper()
    helper.form_class = 'form-horizontal'
    helper.label_class = 'col-12 col-md-3 text-md-right'
    helper.field_class = 'col-12 col-md-9 text-left'
    return helper


class FormHelperHelper:

    def __init__(self):
        pass

    instance: 'FormHelperHelper'

    @property
    def horizontal_nested(self) -> FormHelper:
        helper = self.horizontal
        helper.form_tag = False
        return helper

    @property
    def horizontal(self) -> FormHelper:
        helper = FormHelper()
        helper.form_tag = False
        helper.form_class = 'form-horizontal'
        helper.label_class = 'col-12 col-md-3 text-md-right'
        helper.field_class = 'col-12 col-md-9 text-left'
        return helper

    @property
    def horizontal_no_labels(self) -> FormHelper:
        helper = FormHelper()
        helper.form_tag = False
        helper.form_class = 'form-horizontal'
        helper.field_class = 'col-12'
        helper.form_show_labels = False
        return helper

FORM_HELPER_HELPER = FormHelperHelper()
FormHelperHelper.instance = FORM_HELPER_HELPER