from captcha.fields import ReCaptchaField
from captcha.widgets import ReCaptchaV3
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver
from registration.forms import RegistrationForm

from library.log_utils import AdminNotificationBuilder


class ReCaptchaSignupForm(RegistrationForm):
    captcha = ReCaptchaField(widget=ReCaptchaV3)

    helper = FormHelper()
    helper.add_input(Submit('submit', 'Submit', css_class='btn-primary w-100'))
    helper.form_method = 'POST'


# have put the post_save handler here so we're only notified if registration form is being used

@receiver(post_save, sender=User)
def user_saved(sender, instance: User, **kwargs):
    if not instance.last_login and not instance.is_active:
        nb = AdminNotificationBuilder(message="User Registered")
        nb.add_markdown("A new user has registered")
        nb.add_field("Username", instance.username)
        nb.add_field("Email", instance.email)
        nb.send()
