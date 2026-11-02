"""
PreToolUse hook for Edit/Write: refuses a CLAUDE.local.md. While one exists anywhere from the working directory
up, Claude Code skips every AGENTS.md in the repo (AGENTS.md, This box), so the file is never wanted here.
"""
import json
import os
import sys


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
    return 0


if __name__ == "__main__":
    sys.exit(main())
