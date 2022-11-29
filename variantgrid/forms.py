from captcha.fields import ReCaptchaField
from captcha.widgets import ReCaptchaV3
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from registration.forms import RegistrationForm


class ReCaptchaSignupForm(RegistrationForm):
    captcha = ReCaptchaField(widget=ReCaptchaV3)

    helper = FormHelper()
    helper.add_input(Submit('submit', 'Submit', css_class='btn-primary w-100'))
    helper.form_method = 'POST'
