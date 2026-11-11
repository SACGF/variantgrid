"""
PreToolUse hook for Edit/Write:
  CLAUDE.local.md -> refused. While one exists anywhere from the working directory up, Claude Code skips every
      AGENTS.md in the repo (AGENTS.md, This box), so the file is never wanted here.
  a migration or a one_off_* command, in any app -> once per session, points at manual/AGENTS.md. Those rules live
      there rather than in the root AGENTS.md, and an app's AGENTS.md only loads for work under its own directory.
"""
import json
import os
import re
import sys
import tempfile
from pathlib import Path

MIGRATION_RULES_RE = re.compile(r"/migrations/[^/]+\.py$|/management/commands/one_off_[^/]+\.py$")
REMINDED_DIR = Path(tempfile.gettempdir()) / "vg_migration_rules_reminders"


def _migration_rules_reminder(path: str, session_id: str) -> str:
    if not MIGRATION_RULES_RE.search(path):
        return ""
    marker = REMINDED_DIR / re.sub(r"[^A-Za-z0-9_-]+", "_", session_id or "no_session")
    if marker.exists():
        return ""
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.touch()
    return ("Migrations, ManualOperations and one-off fix/backfill commands follow manual/AGENTS.md (editing a "
            "pushed migration, registering a deploy step, how a one-off is shaped). Read it before this edit.")


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except ValueError:
        return 0
    path = (payload.get("tool_input") or {}).get("file_path") or ""
    if os.path.basename(path) == "CLAUDE.local.md":
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": "AGENTS.md: a CLAUDE.local.md makes Claude Code skip every AGENTS.md in the repo. "
                                        "Per-host facts go in the This box table.",
        }}))
    elif reminder := _migration_rules_reminder(path, payload.get("session_id") or ""):
        print(json.dumps({"hookSpecificOutput": {"hookEventName": "PreToolUse", "additionalContext": reminder}}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
