"""
`manage.py vg` — introspection commands for driving VariantGrid from an agent (or a terminal).

Owns the implementation behind snpdb/management/commands/vg.py: generated maps of the codebase
(`vg map`), test selection from changed files (`vg tests --changed`) and page rendering through the
Django test client (`vg page`). Everything here is read-only against the database.

Modules that need no Django (repo, import_graph, test_selection, markdown, maps.signals, maps.settings)
are kept import-light so scripts/vg can run them in well under a second without booting the project.
See claude/plans/agent_system.md §4.2 for the command contract.
"""
