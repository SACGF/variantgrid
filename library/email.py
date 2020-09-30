from django.core.mail import send_mail


class Email:

    def __init__(self,
                 from_email: str,
                 to_email: str,
                 subject: str,
                 text: str,
                 html: str = None):
        self.from_email = from_email
        self.to_email = to_email
        self.subject = subject
        self.text = text
        self.html = html

    def send(self):
        send_mail(
            subject=self.subject,
            message=self.text,
            from_email=self.from_email,
            recipient_list=[self.to_email],
            fail_silently=False,
            html_message=self.html)
