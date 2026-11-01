---
name: handoff
description: Write the text the user pastes into a fresh session after clearing context - "give me text to copy/paste after clearing", "brief handoff", "I'm going to clear context". Puts everything durable in the repo, the issue or the PR first, then writes a half page that points at them.
---

# Handoff to a new session

The chat is gone when the context clears, so the handoff is two steps: make the state durable, then describe it.

1. **Make it durable before writing a word of prose.**
   - Uncommitted work the user wants kept: commit and push it on the task's branch (only if the user asked for
     commits in this session; otherwise say exactly which files are uncommitted and what they hold).
   - Results, measurements and decisions: a comment on the issue (prefaced 🤖 Written by Claude) or the PR body.
   - Helper scripts worth keeping: in the repo (`scripts/`, or the app) and committed; a scratch directory does not
     survive.
   - A fact the next session could not work out from the code: the doc it belongs in (`AGENTS.md`, `claude/guides/`,
     a module docstring), in the same change.
   - Background work still running: list each shell or agent, its PID or task and what it is doing; stop anything
     that must not outlive this session (a runserver, a benchmark) unless the user wants it kept.
2. **Write the handoff**, in a fenced block the user can copy, at most half a page:
   - **Where**: host, repo path, branch and commit, the issue and PR numbers, the worktree if any.
   - **Done**: what landed, in one line each, each with where it is recorded (commit, PR, issue comment).
   - **Next**: the remaining steps in order, each as the command or the edit to make, plus what "done" looks like.
   - **Decided**: choices the user made in this session that the next session must not reopen, one line each.
   - **Open**: anything waiting on the user, and anything still running on the box.
   - **Read first**: the two or three files or comments the next session should open before doing anything.

Write it as instructions to a competent colleague, not a diary: present tense, no history of what was tried, no
mention of this session. If a durable step in part 1 was impossible (the user hasn't asked for a commit), the handoff
says so in its Open section rather than hoping the working tree survives.
