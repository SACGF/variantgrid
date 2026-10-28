"""
Log the headless browser in as `claude_agent`: creates a session in the configured session engine and writes
its key to $VG_BROWSER_HOME/sessionid, where shot.py and jsq.py read it.

    python3 manage.py shell -c "$(cat scripts/browser/agent_session.py)"
"""
import os
from pathlib import Path

from library.vg.page import agent_session_key, create_agent_user

user, _ = create_agent_user()
key = agent_session_key(user, 7 * 24 * 3600)

home = Path(os.environ.get("VG_BROWSER_HOME", Path.home() / ".cache" / "claude-browser"))
home.mkdir(parents=True, exist_ok=True)
(home / "sessionid").write_text(key)
print(f"session for {user.username} written to {home / 'sessionid'}")
