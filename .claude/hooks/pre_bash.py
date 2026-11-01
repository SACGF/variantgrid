"""
PreToolUse hook for Bash: the rules in AGENTS.md "This box" and "Git commits" that a prompt or a refusal can
enforce, whatever the permission mode.

- GUARDED: state-changing commands on a box other people are testing on get a confirmation prompt with a reason.
  ALLOWED lists the exceptions AGENTS.md grants a host (vg-test2's upgrade restarts), checked first.
- DENIED: commands no box wants - attribution trailers on a commit, `pkill -f` (it matches the shell running it),
  issue searches across the whole org, a CLAUDE.local.md - are refused with the rule and what to do instead.

Everything else passes through untouched. A box set aside for agents (a developer box) opts out of GUARDED only by
setting VG_AGENT_UNGUARDED=1 in the "env" block of that user's ~/.claude/settings.json; DENIED applies everywhere.
"""
import json
import os
import re
import socket
import sys

# (host pattern, command pattern): sanctioned on that host, so no prompt even though GUARDED would match
ALLOWED = [
    (r"^vg-test2$", r"^\s*sudo\s+-n\s+systemctl\s+restart\s+(gunicorn|celeryd_\w+)\s*$"),
]

# (pattern, why it needs a human)
GUARDED = [
    (r"\brestart_services\.sh|\bstop_services\.sh|\bsystemctl\s+(restart|stop|start)\b|\bservice\s+\S+\s+(restart|stop|start)\b",
     "restarts or stops services on this box: anyone using the site or a running import/annotation loses their work"),
    (r"manage\.py\s+migrate\b",
     "migrates the live database this box's testers share; pushed migrations are frozen (AGENTS.md), so check this is what is wanted"),
    (r"manage\.py\s+(vep_run|create_new_variant_annotation_version|annotate_variants|liftover_alleles|delete_variant_annotation_version)\b",
     "changes annotation state for every variant on the box and runs for hours; annotation runs are shared"),
    (r"\b(DROP|TRUNCATE)\s+(TABLE|DATABASE)\b|\bdropdb\b",
     "drops data on the shared database"),
    (r"\b(UPDATE|DELETE\s+FROM|INSERT\s+INTO)\s+(snpdb_variant|snpdb_allele|annotation_variantannotation)\b",
     "writes to a core table every user of the box shares (AGENTS.md, This box)"),
    (r"--js-writes-to-database\b",
     "loads the page in a browser as claude_agent and its ajax writes to the shared database, not rolled back (AGENTS.md, Verifying pages)"),
    (r"\bgit\s+push\s+(-f|--force)|\bgit\s+reset\s+--hard\b|\bgit\s+checkout\s+--\s|\bgit\s+clean\s+-[a-z]*f",
     "discards history or working-tree changes"),
    (r"(^|[;&|]\s*)git\s+pull(?=\s|$)",
     "is a bare git pull: a deployment is upgraded with scripts/upgrade.sh --quick, which also migrates and runs collectstatic (AGENTS.md, This box)"),
    (r"manage\.py\s+test\b(?!.*--keepdb)",
     "runs tests without --keepdb, so the test database is rebuilt from every migration (minutes) and the box's testers share the database server"),
]

# (pattern, the rule and what to do instead)
DENIED = [
    (r"\bgit\s+commit\b.*(Co-Authored-By|Claude-Session)",
     "AGENTS.md (Git commits): no attribution trailers. Commit again with the message ending at its last line of content."),
    (r"\bpkill\s+(-\S+\s+)*-f\b|\bpkill\s+-[a-zA-Z]*f",
     "pkill -f matches the shell running it (its command line holds the pattern) and kills this hook's own shell. "
     "Find the PID (ps -eo pid,args | grep <pattern> | grep -v grep) and kill that."),
    (r"\bgh\s+search\s+\S+.*--owner\s+SACGF\b",
     "AGENTS.md (Which repo an issue goes in): private repos hold patient and lab details, so search only the public repo: "
     "gh search issues --repo SACGF/variantgrid."),
    (r"CLAUDE\.local\.md",
     "AGENTS.md: while a CLAUDE.local.md exists anywhere from the working directory up, Claude Code skips every AGENTS.md "
     "in the repo. Per-host facts go in the This box table."),
]


def decision(command: str, host: str, unguarded: bool) -> tuple[str, str] | None:
    """ ("ask" | "deny", reason) for a command, or None to pass it through """
    for pattern, reason in DENIED:
        if re.search(pattern, command, re.DOTALL):
            return "deny", reason
    if unguarded:
        return None
    for host_pattern, command_pattern in ALLOWED:
        if re.search(host_pattern, host) and all(re.search(command_pattern, line) for line in command.strip().splitlines()):
            return None
    for pattern, reason in GUARDED:
        if re.search(pattern, command, re.MULTILINE):
            return "ask", f"{host} - this {reason}."
    return None


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except ValueError:
        return 0
    command = (payload.get("tool_input") or {}).get("command") or ""
    verdict = decision(command, socket.gethostname(), os.environ.get("VG_AGENT_UNGUARDED") == "1")
    if verdict:
        permission, reason = verdict
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": permission,
            "permissionDecisionReason": reason,
        }}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
