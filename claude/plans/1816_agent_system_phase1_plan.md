# Agent-legible VariantGrid — Phase 1

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-02

Issue: https://github.com/SACGF/variantgrid/issues/1816
Parent plan: [agent_system.md](agent_system.md) §4.1–4.3, §6 "Phase 1".
Status: implemented 2026-09-02 (uncommitted)

Phase 0 (CLAUDE.md restructure, `.claude/` hooks, `vg status`) has not been started, so this phase also
lays down the `vg` command scaffold it would have created. The root CLAUDE.md restructure stays out of
scope here; the nested app notes are written so the root file can later route to them.

## Deliverables

1. `manage.py vg` command family scaffold and three subcommands: `map`, `page`, `tests`.
2. `claude/maps/` — six generated maps, committed; `vg map --check` fails when stale; a DB-free CI job
   runs the check.
3. `vg page` renders a URL through the Django test client against the live DB as the `claude_agent`
   user, with `--queries` for production query count and N+1 detection.
4. `vg tests --changed` selects test modules from changed files through the first-party import graph.
5. `<app>/CLAUDE.md` for snpdb, genes, annotation, analysis, classification, upload, uicore, library.
6. `claude/guides/testing.md` fixture index.
7. `--parallel` for the suite: tried locally, result recorded below.

## Data

No Django models. The generators produce plain dicts that the Markdown and JSON renderers share:

```python
@dataclass
class MapTable:            # one map = one or more tables
    title: str             # e.g. "snpdb"
    columns: list[str]
    rows: list[list[str]]
    anchor: str            # heading slug, so docs can link claude/maps/models.md#snpdb
```

```python
@dataclass
class PageResult:
    status: int
    templates: list[str]
    production_queries: int | None
    repeated: list[tuple[str, int]]   # normalised SQL, count (only > 1)
    outline: list[str]                # "h1 …", "table 25×12 (…)", "form name: field, field", "link text → href"
```

## Layout

```
library/vg/                     pure helpers; Django-free unless noted
    __init__.py                 contract docstring
    repo.py                     repo root, first-party top-level packages, git changed files
    markdown.py                 MapTable → Markdown / JSON
    import_graph.py             AST import graph over first-party modules (no Django)
    test_selection.py           changed files → test labels (no Django)
    page.py                     test-client render + outline + query stats (Django)
    maps/
        __init__.py             registry: name → generator, and check()
        models.py               apps registry (Django)
        urls.py                 resolver walk (Django)
        commands.py             get_commands() + AST for models touched (Django)
        tasks.py                AST for decorators/Task subclasses; enqueuers; queue from CELERY_TASK_ROUTES (settings)
        signals.py              AST: Signal() definitions, @receiver / .connect receivers
        settings.py             AST over components/ and env/: setting → default → env overrides
snpdb/management/commands/vg.py subparsers: map, page, tests
scripts/vg                      thin dispatcher: `tests` and `map signals|settings` run without Django
claude/maps/{models,urls,commands,tasks,signals,settings}.md
.github/workflows/agent-maps.yml   DB-free job: vg map --check (runs on every push, including md/claude changes)
```

## Decisions

- **Where the command lives:** `snpdb/management/commands/vg.py`, next to `create_fake_data`.
  `library/` is not an app and `variantgrid/` is not in INSTALLED_APPS. All logic is in `library/vg/`
  so the command file is only argument parsing.
- **Maps are DB-free.** `models` and `urls` need Django set up but never query. The CI job boots
  Django with the `github_actions` settings and no Postgres service; anything that queries at import
  time is a bug the job will surface.
- **Templates per URL** come from `template_name` on class-based views, or a regex over the view's
  source for `*.html` literals for function views. Approximate by design — the map says so.
- **Task queue** comes from the decorator's `queue=` kwarg, else `CELERY_TASK_ROUTES`, else the
  default queue. "Enqueued by" is any module that calls `<task>.delay/apply_async/s/si/signature`.
- **`claude_agent` user** is created on demand by `vg page --create-user` (member of `all_users`, no
  lab, unusable password). `vg page` otherwise never writes; the render runs inside a transaction that
  is rolled back so a view with a write side-effect cannot dirty the lab.
- **Test selection rules:** a test module is selected if it transitively imports a changed module;
  a changed test module selects itself; a changed template / JS / `urls.py` selects the app's
  `tests.test_urls`; a changed migration selects the app's whole `tests` package; anything under
  `variantgrid/settings` selects nothing on its own (the graph carries settings imports anyway);
  `library/` changes fan out naturally through the graph.
- **`--parallel`:** Django's parallel runner with `--keepdb` clones `test_snpdb` per worker. Result of
  the local trial goes in the "Outcome" section below and in `claude/guides/testing.md`.

## Outcome

- `manage.py vg map|page|tests` and `scripts/vg` landed; six maps committed under `claude/maps/`;
  `.github/workflows/agent-maps.yml` runs `vg map --check` without a database on every push.
- `claude_agent` user created on vg-test2; `vg page /variantopedia/dashboard --queries` renders in ~8s
  including the Django boot and reports 16 production queries.
- `scripts/vg tests --explain` resolves the import graph in 0.6s warm (3s cold, cached in `.vg_cache/`).
- `--parallel 4 --keepdb`: 2,741 tests in 103s (1m59s wall) on vg-test2, all passing, no serial-only tests
  found. CI now runs `--parallel 4`.
- Eight `<app>/CLAUDE.md` files and `claude/guides/testing.md` written; each cited symbol was grep-verified.
  The drafting surfaced ~40 stale claims in `claude/research/*.md` (field names, removed models, split
  views) - the research docs are due the Phase 2 "narrative + Verified against" rewrite.
- Not done: `seqauto_api_client2` fails to import (`dataclasses_json` missing from requirements), so the
  commands map lists it with category `?`. Root CLAUDE.md only gained pointers; its Phase 0 restructure
  (Start here table, This box, memory policy, definition of done) is still open.
