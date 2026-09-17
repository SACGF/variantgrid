"""
`manage.py vg` — introspection commands for driving VariantGrid from an agent (or a terminal).

Owns the implementation behind snpdb/management/commands/vg.py: what is running on this box (`vg status`),
generated maps of the codebase (`vg map`), test selection from changed files (`vg tests --changed`), page
rendering through the Django test client (`vg page`), module outlines (`vg outline`) and where a setting's
value came from (`vg settings`). Everything here is read-only against the database.

Modules that need no Django (repo, import_graph, test_selection, markdown, outline, settings_chain,
maps.signals, maps.settings, maps.tasks) are kept import-light so scripts/vg can run them in well under a
second without booting the project. See claude/plans/agent_system.md §4.2 for the command contract.
"""
