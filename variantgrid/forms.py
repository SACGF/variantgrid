from captcha.fields import ReCaptchaField
from registration.forms import RegistrationForm


class ReCaptchaSignupForm(RegistrationForm):
    captcha = ReCaptchaField()
