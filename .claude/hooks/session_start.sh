#!/bin/bash
# SessionStart hook: orientation for free, no Django boot. Printed into the session's context.
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}" || exit 0
echo "host: $(hostname)  branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null)@$(git rev-parse --short HEAD 2>/dev/null)"
dirty=$(git status --short 2>/dev/null | head -10)
if [ -n "$dirty" ]; then
  echo "working tree (first 10):"
  echo "$dirty" | sed 's/^/  /'
else
  echo "working tree: clean"
fi
newest_plan=$(ls -t claude/plans/*.md 2>/dev/null | head -1)
[ -n "$newest_plan" ] && echo "newest plan: $newest_plan ($(sed -n 's/^Status: *//p' "$newest_plan" | head -1))"
echo "orient with: python3 manage.py vg status · scripts/vg tests --explain · vg page <url> --queries (CLAUDE.md 'Start here')"
