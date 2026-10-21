# claude/plans — agent notes

## Plans
Plans live in `claude/plans/<issue>_<slug>_plan.md`. Directly under the title, record which Claude model wrote it, e.g.
`Written by Claude Fable 5 (claude-fable-5), 2026-08-31` - so when a plan is picked up later it is clear which model's
judgement it reflects. Update the line if a different model revises the plan. Add a `Status:` line
(`draft | approved | in progress | landed <sha> | superseded by <plan>`) and keep it current.

A plan states the chosen design: what the code and data will be, and why. Leave out rejected alternatives,
"considered and dropped" sections and revision history (what changed since the last draft, review notes, struck-through
text) - when the design changes, rewrite the text in place so it reads as the current decision. Once part of a plan
lands, describe that part as what exists now. Where the obvious approach is tempting and still wrong, one line giving
the present reason is enough.

Put the data front and centre. Code can be changed later; data stays in the database for years and limits what can be
built on it, so the database models are what the reviewer most wants to see. When a plan adds or changes a Django model,
show the model as a code block with just its fields, relations, constraints and `Meta` - near the top of the plan, before
the code that uses it. Same for a dataclass or other data holder: show the member variables only. Leave methods and
properties out of the plan; they belong in the implementation.

`scripts/vg docs check` verifies a plan's citations while its `Status:` is draft, approved or in progress. A landed plan
whose knowledge has moved into docs is deleted. A backticked path is a claim the file exists, so write a file the plan
will create as plain text (analysis/tasks/new_tasks.py (new file)).

## Plans and their issues
`claude/plans/` is the record of work still to do, so a plan is deleted once it is done - landed, superseded or
abandoned - with any lasting knowledge moved into docs first. When a plan is written for a GitHub issue, comment on the
issue with the plan's path (prefaced 🤖 Written by Claude, as for any GitHub comment) and add the `Has Plan` label:
`gh issue comment N --body "🤖 Written by Claude: plan at claude/plans/<file>"` and
`gh issue edit N --add-label "Has Plan"`. When the plan is deleted, remove the label
(`gh issue edit N --remove-label "Has Plan"`).

## Implementation prompts
When asked to draft a prompt for an agent to implement a plan in another conversation:
- The plan file is the spec. Reference it; don't restate it.
- Phrase everything positively. Do not include "do not", "don't", "no X", or any "Constraints" section listing things to
  avoid - even for defaults the agent would otherwise do, and even for ideas that came up and were rejected during
  planning. Naming the unwanted thing plants it ("don't think of an elephant"). If a default needs to be overridden,
  either fix the plan to carry the positive instruction, or state the positive behaviour you want ("update all callers to
  use the new kwarg" rather than "don't add a backwards-compat shim").
- The plan reflects the final decision; the agent reading it won't see the alternatives. Mentioning rejected options only
  confuses or implies the plan is incomplete.
- Keep prompts short: read-list, "follow plan §X-§Y", any positive overrides, report-back format. No "pre-resolved
  decisions" section.
