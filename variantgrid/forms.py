from captcha.fields import ReCaptchaField
from captcha.widgets import ReCaptchaV3
from registration.forms import RegistrationForm


class ReCaptchaSignupForm(RegistrationForm):
    captcha = ReCaptchaField(widget=ReCaptchaV3)
