# paper/ — VariantGrid paper working folder

Working materials for the VariantGrid platform paper. Last updated 2026-06-12.

## Start here
- **`todo.md`** — the sequenced action plan (VG4 upgrade → cdot → kiosk → data-mining → write).
- **`strategy.md`** — the strategy guide: the angle, scope (in/out), differentiators, the
  data-mining results plan, Talos positioning, timing, and the "where to spend effort" decision.
  **Read this first.**
- **`journal_targets.md`** — ranked candidate journals + the public-access / kiosk analysis.
- **`data_mining_opportunities.md`** — the real-world usage-data mining ideas (the empirical
  Results section); the unique evidence no competitor can produce.
- **`data_collection_plan.md`** — concrete plan: datasets → paper claim → source model/field,
  governance/safety rules, and the VG3→VG4 approach.
- **`scripts/`** — read-only aggregate mining scripts (`paper_data_mining.py`) for **you** to run
  (Claude does not run on the clinical environment).

## Literature (`literature/`)
- `README.md` — index of the literature folder.
- `comparison.md` — the full landscape analysis: comparator-by-comparator, what's novel, what to
  soften, scope decisions, Talos.
- `feature_comparison.csv` — feature matrix (~67 features × 9 tools); VG column source-verified.
- `feature_comparison_notes.md` — how to read the CSV; STRONG/PARITY/BEHIND groupings.
- `paywalled_articles.md` — references to retrieve manually (ACMG 2015, ClinGen 2015, …).
- `pdfs/` + `text/` — downloaded papers (full text) and web extracts. `convert.py` re-runs the
  PDF→text conversion.

## Reference
- `2018_paper.txt` — the old, out-of-date draft (kept for the abstract/framing seeds).

## The angle in one line
A configurable, clinically-validated, self-hostable platform whose data model enables continuous
ingestion, versioned annotation, and automated re-analysis — with live node-graph analysis — proven
by years of real clinical/research usage. Sharing/discordance → cite Shariant; HGVS → cite cdot
(publish first). Don't chase feature parity; spend effort on the usage-data figures + a cached
public demo (which unlocks NAR Web Server).
