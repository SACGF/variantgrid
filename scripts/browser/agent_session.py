"""
Log the headless browser in as `claude_agent`: creates a session in the configured session engine (the
cache backend here, so a db SessionStore would never be found) and writes its key to
$VG_BROWSER_HOME/sessionid, where shot.py and jsq.py read it.

    python3 manage.py shell -c "$(cat scripts/browser/agent_session.py)"
"""
import os
from importlib import import_module
from pathlib import Path

from django.conf import settings
from django.contrib.auth import BACKEND_SESSION_KEY, HASH_SESSION_KEY, SESSION_KEY

from library.vg.page import create_agent_user

user, _ = create_agent_user()
session = import_module(settings.SESSION_ENGINE).SessionStore()
session[SESSION_KEY] = str(user.pk)
session[BACKEND_SESSION_KEY] = "django.contrib.auth.backends.ModelBackend"
session[HASH_SESSION_KEY] = user.get_session_auth_hash()
session.set_expiry(7 * 24 * 3600)
session.create()

home = Path(os.environ.get("VG_BROWSER_HOME", Path.home() / ".cache" / "claude-browser"))
home.mkdir(parents=True, exist_ok=True)
(home / "sessionid").write_text(session.session_key)
print(f"session for {user.username} written to {home / 'sessionid'}")
