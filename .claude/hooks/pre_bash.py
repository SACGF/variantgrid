"""
PreToolUse hook for Bash: state-changing commands on a box other people are testing on get a
confirmation prompt with a reason, whatever the permission mode (agent_system.md §4.5). Everything
else passes through untouched.
"""
import json
import re
import sys

# (pattern, why it needs a human)
GUARDED = [
    (r"\brestart_services\.sh|\bstop_services\.sh|\bsystemctl\s+(restart|stop|start)\b|\bservice\s+\S+\s+(restart|stop|start)\b",
     "restarts or stops services on this box: anyone using the site or a running import/annotation loses their work"),
    (r"manage\.py\s+migrate\b",
     "migrates the live database this box's testers share; pushed migrations are frozen (CLAUDE.md), so check this is what is wanted"),
    (r"manage\.py\s+(vep_run|create_new_variant_annotation_version|annotate_variants|liftover_alleles|delete_variant_annotation_version)\b",
     "changes annotation state for every variant on the box and runs for hours; annotation runs are shared"),
    (r"\b(DROP|TRUNCATE)\s+(TABLE|DATABASE)\b|\bdropdb\b",
     "drops data on the shared database"),
    (r"\bgit\s+push\s+(-f|--force)|\bgit\s+reset\s+--hard\b|\bgit\s+checkout\s+--\s|\bgit\s+clean\s+-[a-z]*f",
     "discards history or working-tree changes"),
]


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except ValueError:
        return 0
    command = (payload.get("tool_input") or {}).get("command") or ""
    for pattern, reason in GUARDED:
        if re.search(pattern, command):
            print(json.dumps({"hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "permissionDecision": "ask",
                "permissionDecisionReason": f"vg-test2 is shared with human testers - this {reason}.",
            }}))
            return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
