"""
PostToolUse hook for Bash: after a `gh issue view N`, or any command naming an issue by URL or by its
`repos/<owner>/<repo>/issues/N` API path, add a line naming the work that already exists for that issue - PRs with
their branch and state, `<issue>_<slug>` branches, a plan, commits on master - so an open PR is seen before any work
is proposed or written. `scripts/vg issue N` gives
the whole picture; this only makes sure the plain `gh` view can't hide it.

After a `git worktree add`, once per session, it points at claude/guides/github.md#worktrees: the gitignored files a
new worktree lacks (settings module, reverse.js, its own test database, collectstatic) are listed there.

The lookup runs `gh api` (a few seconds). Anything that fails - no network, no gh auth, an unparseable command -
is silent: the hook never alters the command's own result. A private-repo issue (the user handed it over) gets the
same line, listing linked items by number and state only.
"""
import json
import os
import re
import sys
import tempfile
from pathlib import Path

REPO = os.environ.get("CLAUDE_PROJECT_DIR") or os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from library.vg.issue import existing_work_summary, origin_repo, parse_issue_ref  # noqa: E402

ISSUE_VIEW_RE = re.compile(r"\bgh\s+issue\s+view\s+(?P<ref>\S+)")
REPO_FLAG_RE = re.compile(r"(?:--repo|-R)[\s=](?P<repo>[\w.-]+/[\w.-]+)")
ISSUE_URL_RE = re.compile(r"(?:github\.com|\brepos)/[\w.-]+/[\w.-]+/issues/\d+")
WORKTREE_ADD_RE = re.compile(r"\bgit\s+worktree\s+add\b")
WORKTREE_REMINDED_DIR = Path(tempfile.gettempdir()) / "vg_agent_rule_reminders"
WORKTREE_REMINDER = ("A new worktree shares the virtualenv but not gitignored files: set it up as "
                     "claude/guides/github.md#worktrees says (settings module, collectstatic_js_reverse, its own test "
                     "database, collectstatic) and `git worktree remove` it when the PR is up.")


def worktree_reminder(command: str, session_id: str) -> str:
    """ Once per session, after the first `git worktree add` """
    if not WORKTREE_ADD_RE.search(command):
        return ""
    marker = WORKTREE_REMINDED_DIR / (re.sub(r"[^A-Za-z0-9_-]+", "_", session_id or "no_session") + "_worktree")
    if marker.exists():
        return ""
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.touch()
    return WORKTREE_REMINDER


def issue_ref_in(command: str):
    """ The issue a `gh issue view` or an issue URL in the command refers to, or None when it names none """
    if re.search(r"\bvg\s+issue\b", command):
        return None
    if match := ISSUE_VIEW_RE.search(command):
        ref = match["ref"]
        if ref.startswith("-"):
            return None
        repo_flag = REPO_FLAG_RE.search(command)
        return parse_issue_ref(ref, default_repo=repo_flag["repo"] if repo_flag else origin_repo())
    if match := ISSUE_URL_RE.search(command):
        return parse_issue_ref(match.group(0))
    return None


def main() -> int:
    try:
        payload = json.load(sys.stdin)
        command = (payload.get("tool_input") or {}).get("command") or ""
        ref = issue_ref_in(command)
        if ref is None:
            summary = worktree_reminder(command, payload.get("session_id") or "")
            if not summary:
                return 0
        else:
            summary = existing_work_summary(ref)
    except Exception:  # pylint: disable=broad-exception-caught
        return 0
    print(json.dumps({"hookSpecificOutput": {"hookEventName": "PostToolUse", "additionalContext": summary}}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
