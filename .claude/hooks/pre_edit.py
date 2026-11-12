"""
PreToolUse hook for Edit/Write:
  CLAUDE.local.md -> refused. While one exists anywhere from the working directory up, Claude Code skips every
      AGENTS.md in the repo (AGENTS.md, This box), so the file is never wanted here.
  a file matching one of REMINDERS, in any app -> once per session per reminder, points at the doc holding the rules
      for that kind of file. Those rules live there rather than in the root AGENTS.md, and an app's AGENTS.md only
      loads for work under its own directory, so without the reminder they would never be read.
"""
import json
import os
import re
import sys
import tempfile
from pathlib import Path

# (name, path pattern, reminder): the name keys the once-per-session marker
REMINDERS = [
    ("migrations", re.compile(r"/migrations/[^/]+\.py$|/management/commands/one_off_[^/]+\.py$"),
     "Migrations, ManualOperations and one-off fix/backfill commands follow manual/AGENTS.md (editing a pushed "
     "migration, registering a deploy step, how a one-off is shaped). Read it before this edit."),
    ("tests", re.compile(r"/tests/.+\.py$"),
     "Tests follow claude/guides/testing.md: which tests to keep (#which-tests-to-keep - no framework or "
     "change-detector tests), the base classes and the fixture builders. Read it before writing a test."),
    ("imports", re.compile(r"/(models|enums)/[^/]+\.py$"),
     "Imports in a models/ or enums/ package follow library/AGENTS.md#imports: import from the defining module, not "
     "a package __init__; an enums module never imports an app's models; the app layering lint-imports enforces."),
    ("templates", re.compile(r"/templates/.+\.html$"),
     "Templates follow variantgrid/static_files/AGENTS.md: Bootstrap 4 data-toggle/data-target (not data-bs-*), and "
     "page JS in a static file (a JS block in a template you edit moves out with the change)."),
    ("pages", re.compile(r"/templates/.+\.html$|/views(/[^/]+)?\.py$|/views_[^/]+\.py$"),
     "A template or view change: run `vg page <url> --queries` before this edit and after it - a rising query count "
     "is an N+1 (uicore/AGENTS.md#verifying-pages)."),
    ("celery", re.compile(r"/tasks(/[^/]+)?\.py$|/tasks_[^/]+\.py$"),
     "Celery tasks: one with no queue= or CELERY_TASK_ROUTES entry runs on db_workers, which suits most; the queues "
     "and when to route elsewhere are in claude/guides/operations.md#services-queues-logs."),
]
REMINDED_DIR = Path(tempfile.gettempdir()) / "vg_agent_rule_reminders"


def _reminders(path: str, session_id: str) -> str:
    session = re.sub(r"[^A-Za-z0-9_-]+", "_", session_id or "no_session")
    messages = []
    for name, pattern, message in REMINDERS:
        if not pattern.search(path):
            continue
        marker = REMINDED_DIR / f"{session}_{name}"
        if marker.exists():
            continue
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.touch()
        messages.append(message)
    return "\n".join(messages)


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
    elif reminder := _reminders(path, payload.get("session_id") or ""):
        print(json.dumps({"hookSpecificOutput": {"hookEventName": "PreToolUse", "additionalContext": reminder}}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
