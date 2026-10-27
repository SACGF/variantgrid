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
# claude/runbooks/upstream_watch.md: a row per fork / patched copy we carry, with its "Next check" date in column 5
due=$(awk -F'|' -v today="$(date +%F)" '/^\|/ && $5 ~ /[0-9]{4}-[0-9]{2}-[0-9]{2}/ { gsub(/ /, "", $5); if ($5 <= today) n++ } END { print n + 0 }' claude/runbooks/upstream_watch.md 2>/dev/null)
[ "${due:-0}" -gt 0 ] && echo "upstream watch: $due check(s) due - claude/runbooks/upstream_watch.md"
newest_plan=$(ls -t claude/plans/*.md 2>/dev/null | head -1)
[ -n "$newest_plan" ] && echo "newest plan: $newest_plan ($(sed -n 's/^Status: *//p' "$newest_plan" | head -1))"
# claude/maps/*.md are gitignored: regenerate in the background so they are fresh and never in a diff.
(scripts/vg map >/dev/null 2>&1 &)
echo "orient with: python3 manage.py vg status · scripts/vg tests --explain · vg page <url> --queries (AGENTS.md 'Start here')"
