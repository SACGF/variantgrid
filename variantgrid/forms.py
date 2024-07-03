from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django_recaptcha.fields import ReCaptchaField
from django_recaptcha.widgets import ReCaptchaV3
from registration.forms import RegistrationForm


class ReCaptchaSignupForm(RegistrationForm):
    captcha = ReCaptchaField(widget=ReCaptchaV3)

    helper = FormHelper()
    # name can't be 'submit' see https://github.com/django-recaptcha/django-recaptcha/issues/203#issuecomment-700997270
    helper.add_input(Submit('my_submit', 'Submit', css_class='btn-primary w-100'))
    helper.form_method = 'POST'
