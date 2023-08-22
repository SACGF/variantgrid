from crispy_forms.helper import FormHelper


def form_helper_horizontal() -> FormHelper:
    helper = FormHelper()
    helper.form_class = 'form-horizontal'
    helper.label_class = 'col-12 col-md-3 text-md-right'
    helper.field_class = 'col-12 col-md-9 text-left'
    return helper


class FormHelperHelper:

    """
    Form helper used for CrispyForms
    I'd like to rework this so the methods are cummalitive, and at the end you call a method to generate the actual form helper
    e.g. form_hh.horizontal.no_labels.disable_csrf.helper()
    """

    def __init__(self, disable_csrf: bool = False):
        self.disable_csrf = disable_csrf
        pass

    instance: 'FormHelperHelper'

    def no_csrf(self) -> 'FormHelperHelper':
        return FormHelperHelper(disable_csrf=True)

    def _base_form_helper(self):
        fh = FormHelper()
        fh.disable_csrf = self.disable_csrf
        return fh

    @property
    def horizontal_nested(self) -> FormHelper:
        helper = self.horizontal
        helper.form_tag = False
        return helper

    @property
    def horizontal(self) -> FormHelper:
        helper = self._base_form_helper()
        helper.form_tag = False
        helper.form_class = 'form-horizontal'
        helper.label_class = 'col-12 col-md-3 text-md-right'
        helper.field_class = 'col-12 col-md-9 text-left'
        return helper

    @property
    def horizontal_no_labels(self) -> FormHelper:
        helper = self._base_form_helper()
        helper.form_tag = False
        helper.form_class = 'form-horizontal'
        helper.field_class = 'col-12'
        helper.form_show_labels = False
        return helper

    @property
    def fields_only(self) -> FormHelper:
        helper = self._base_form_helper()
        helper.form_tag = False
        helper.form_show_labels = False
        helper.label_class = "d-none"
        helper.field_class = "col-12"
        return helper


FORM_HELPER_HELPER = FormHelperHelper()
FormHelperHelper.instance = FORM_HELPER_HELPER
