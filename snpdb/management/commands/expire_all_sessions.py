# From https://gist.github.com/tclancy/10269504

import datetime
from django.conf import settings
from django.contrib.auth import logout
from django.contrib.auth.models import User
from django.contrib.sessions.models import Session
from django.core.management.base import BaseCommand
from django.http import HttpRequest
from importlib import import_module


def init_session(session_key):
    """
    Initialize same session as done for ``SessionMiddleware``.
    """
    engine = import_module(settings.SESSION_ENGINE)
    return engine.SessionStore(session_key)


class Command(BaseCommand):
    help = "Kill all active sessions"

    def handle(self, *args, **options):
        """
        Read all available users and all available not expired sessions. Then
        logout from each session. Start with a day ago to hack around the
        timezone issue instead of doing something smart.
        """
        start = datetime.datetime.now() - datetime.timedelta(days=1)
        request = HttpRequest()

        sessions = Session.objects.filter(expire_date__gt=start)

        print('Found %d not-expired session(s).' % len(sessions))

        for session in sessions:
            user_id = session.get_decoded().get('_auth_user_id')
            request.session = init_session(session.session_key)

            logout(request)

            user = User.objects.get(pk=user_id)
            print('Successfully logged out user: %r' % user)

        print('All OK!')
