# #1669 — ReleaseGeneMatcher: single-hop alias matching

Written by Claude Fable 5 (claude-fable-5), 2026-08-31

[#1669](https://github.com/SACGF/variantgrid/issues/1669): `ReleaseGeneMatcher` walks the gene symbol
alias table transitively, so an alias string shared by two unrelated genes bridges them. Reported from
prod (SACGF/variantgrid_sapath#426) as a hypogonadotropic hypogonadism panel returning a **PDCD2**
variant via `MT-TS2 → RP8 → PDCD2`.

The fix is one idea: **a symbol matches a release gene directly, or through exactly one
`GeneSymbolAlias` row that touches a symbol the release knows.** Everything else in this plan is
cleaning up the rows the old traversal wrote.

---

## 1. Diagnosis (confirmed against the code and the local DB)

`genes/gene_matching.py`, `ReleaseGeneMatcher`:

* `aliases_dict` (line 165) builds `alias_graph` with every `GeneSymbolAlias` row hung under **both**
  its `alias` and its `gene_symbol_id` (lines 173–176), then starts a walk from both ends of every row.
* `_aliases()` (line 186) recurses through that graph until it reaches any symbol in `self.genes`
  (symbols the release actually has), then assigns those gene ids to **every symbol on the path**
  (lines 196–211). Path length is unbounded; only loops are guarded.

Local DB reproduces the prod rows exactly — `MT-TS2 → 5134 (PDCD2)` in releases 3, 4, 7, 8 and
`MT-TS2 → ENSG05220024333` in release 5.

Two further findings that shape the cleanup:

1. **`match_info` cannot identify chained rows.** `symbol_match_path` is a dict keyed by symbol, so
   when the walk re-enters a symbol (`RP8 → MT-TS2 → RP8 → PDCD2`, which the loop guard allows because
   the root symbol is never added to `visited_symbols`) the earlier hop is overwritten and the stored
   text collapses to a single hop. Release 3's `MT-TS2` row reads `RP8 is an alias for PDCD2 (HGNC)`,
   so the issue's `LIKE '%), %'` audit query undercounts.
2. **Rows also go stale when the alias table changes.** `genes/cached_web_resource/refseq.py` deletes
   and re-imports NCBI aliases, but `ReleaseGeneSymbolGene` rows derived from the deleted aliases
   survive (`match_symbols_to_genes` only inserts). Release 1 has `ADAM12-OT1 → ADAM12` citing
   `ADAM12-OT1 is an alias for ADAM12 (NCBI)`; the only alias rows for `ADAM12-OT1` today are
   `CAR10` and `FLJ31066` (HGNC).

So the cleanup recomputes what the fixed matcher would produce and diffs against the table, rather
than pattern-matching `match_info`.

### What the prototype single-hop matcher says about existing rows

Alias/gene-version derived rows (`match_info IS NOT NULL`) per release on the local DB, classified by
the shortest explanation the single-hop rules give them:

| release | rows | gene-version | forward hop only | backward hop only | multi-hop / stale |
|---|---|---|---|---|---|
| 1 GRCh38 RefSeq 110 | 1110 | 0 | 129 | 111 | 870 |
| 2 GRCh37 RefSeq 105 | 2012 | 82 | 177 | 113 | 1640 |
| 3 GRCh37 RefSeq 105 | 1502 | 91 | 180 | 113 | 1118 |
| 4 GRCh38 RefSeq 2023 | 1288 | 0 | 299 | 110 | 879 |
| 5 T2T Ensembl 2022-06 | 2249 | 0 | 171 | **1110** | 968 |
| 7 GRCh38 RefSeq 2024-08 | 1275 | 0 | 315 | 112 | 848 |
| 8 GRCh38 RefSeq 2025-08 | 606 | 0 | 128 | 125 | 353 |

The "backward hop only" column is why the hop stays bidirectional (§2). The last column is what the
resync command (§4) removes.

---

## 2. Design: one hop, either direction

`GeneSymbolAlias(alias=A, gene_symbol=S)` reads "A is an alias for S". For a queried symbol `Q` that
is **absent** from the release:

* **Forward** — `Q == A`, `S` in release → `Q` matches `S`'s genes. Gene list says `KAL1`, release
  has `ANOS1`. This is the common "list uses an old name" case.
* **Backward** — `Q == S`, `A` in release → `Q` matches `A`'s genes. Gene list says `ANOS1`, release
  still calls the gene `KAL1`. HGNC is refreshed far more often than gene annotation releases, so a
  current-symbol gene list against an older release hits this constantly — 1,110 rows in the T2T
  Ensembl 2022 release alone (`ABTB3`, `ACE2-DT`, `ACTMAP`, `ADISSP`, …), and ~110 per RefSeq release.

One hop is enough for renames: HGNC's `prev_symbol` lists *every* previous symbol of a gene, so a
gene renamed twice has direct alias rows from both old names. Chaining is what lets an alias string
shared by two genes (`RP8`, `ALP`, `CD`, `ADMR`) act as a bridge, and dropping it stops
`MT-TS2 → RP8 → PDCD2`, `ATHS → ALP → {SLPI, ATRNL1, CCL27, PDLIM3, ASRGL1}`,
`CELIAC2 → CD → {NOD2, CTLA4}`, `ACKR5 → ADMR → GPR182`.

Precedence, unchanged from today's `_get_gene_id_and_match_info_for_symbol`:

1. Symbol is in the release → its genes, `match_info=None`, aliases never consulted.
2. Otherwise the union of gene-version matches (`Gene v0/GRCh38`, from `_get_genes_dict`) and
   single-hop alias matches. Where both name the same gene the gene-version text wins (it is set
   first; alias entries use `setdefault`).

A symbol can still legitimately map to several genes (e.g. `RP8` forward-matches both `MT-TS2` and
`PDCD2`; only `PDCD2` is in the release, so `RP8` → `PDCD2` is correct and kept).

The residual risk of the backward hop — HGNC's informal `alias_symbol` occasionally being another
gene's approved symbol (`AURKAIP1` lists `AIP`) — only bites when the queried symbol is absent from
the release *and* the alias is present as a different gene. Distinguishing `prev_symbol` from
`alias_symbol` is [#1668](https://github.com/SACGF/variantgrid/issues/1668) item 2; once that
provenance exists, the backward hop should accept `prev_symbol` rows only. `ReleaseGeneMatcher`
becomes a consumer of #1668's shared resolver at that point.

---

## 3. `genes/gene_matching.py`

Replace `aliases_dict` and delete `_aliases`:

```python
@cached_property
def aliases_dict(self) -> dict[str, dict]:
    """ Upper-case symbol -> {gene_id: match_info} for symbols the release lacks, reached either
        from an older GeneVersion's symbol or by a single GeneSymbolAlias hop (in either direction)
        to a symbol the release has. One hop is deliberate: HGNC lists every previous symbol of a
        gene, so renames never need chaining, and chaining lets an alias string shared by two
        unrelated genes bridge them (#1669). """
    genes_dict = self._get_genes_dict()
    alias_qs = GeneSymbolAlias.objects.all()
    for gsa in alias_qs:
        alias = gsa.alias.upper()
        symbol = gsa.gene_symbol_id.upper()
        if alias == symbol:
            continue
        for query, target in ((alias, symbol), (symbol, alias)):
            if query in self.genes:
                continue
            for gene_id in self.genes.get(target, []):
                genes_dict[query].setdefault(gene_id, gsa.match_info)
    return genes_dict
```

Notes for the implementer:

* `self.genes` keys are `Upper("gene_symbol")`; `_get_genes_dict` keys are `gene_symbol_id.upper()`.
  Keep everything upper-cased on the Python side — `alias` has `case_insensitive` collation in
  Postgres but that does nothing for dict lookups.
* Loading all alias rows replaces the old `gene_symbol__releasegenesymbol__release=self.release`
  filter. With one hop both ends are checked against `self.genes` directly, so the pre-filter buys
  nothing, and removing it means `aliases_dict` (a `cached_property`) no longer depends on which
  `ReleaseGeneSymbol` rows existed at the moment it was first computed.
* `_get_genes_dict`, `_get_gene_id_and_match_info_for_symbol`, `match_symbols_to_genes` and
  everything below are unchanged. `ReleaseGeneSymbolGene.match_info` keeps its existing format
  (`"X is an alias for Y (HGNC)"` / `"Gene v0/GRCh38"`), so `gene_grid.js` keeps rendering it.

---

## 4. Resync command: `fix_rematch_release_symbols_to_genes`

Extend the existing command (`genes/management/commands/fix_rematch_release_symbols_to_genes.py`)
from add-only to a full resync of the derived table, per release:

1. `release_gene_symbols = list(gar.releasegenesymbol_set.all())`
2. `expected = gm._get_gene_id_and_match_info_for_symbol(rgs.gene_symbol_id for rgs in ...)` —
   a dict `gene_symbol_id -> [(gene_id, match_info)]`.
3. Load existing rows: `ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=gar)
   .values_list("pk", "release_gene_symbol__gene_symbol_id", "gene_id", "match_info")`.
4. Diff on `(gene_symbol_id, gene_id)`:
   * existing pair absent from `expected` → **delete** (`filter(pk__in=...)` in batches of 2000)
   * pair present in both with different `match_info` → **update** (`bulk_update`)
   * pair only in `expected` → **insert** (reuse `gm.match_symbols_to_genes(release_gene_symbols)`
     after the deletes; it already does `ignore_conflicts=True`)
5. Print per-release counts: deleted / updated / inserted, and symbols left with no gene (as now).
   At `--verbosity 2` also print each deleted `(symbol, gene_id, match_info)` so an operator can
   eyeball what went.

Add `--dry-run` (report the diff, make no writes) so it can be run on prod for inspection first.

Cache: `GeneAnnotationRelease.genes_for_symbols` is a plain queryset with no Redis caching, and the
analysis gene-list nodes query `ReleaseGeneSymbolGene` at run time, so the change takes effect on
the next node load. Existing `NodeCount` figures for already-loaded gene-list nodes reflect the old
rows until the node re-runs; that is the normal behaviour after any gene-data change and needs no
extra invalidation.

---

## 5. Migration: `genes/migrations/0089_one_off_resync_release_gene_symbol_genes.py`

Follow `genes/migrations/0065_one_off_fix_make_panel_app_gene_lists_public.py`:

```python
def _has_alias_derived_release_gene_matches(apps):
    ReleaseGeneSymbolGene = apps.get_model("genes", "ReleaseGeneSymbolGene")
    return ReleaseGeneSymbolGene.objects.filter(match_info__contains="is an alias for").exists()

operations = [
    ManualOperation(task_id=ManualOperation.task_id_manage(["fix_rematch_release_symbols_to_genes"]),
                    note="Remove chained-alias gene matches (e.g. MT-TS2 -> RP8 -> PDCD2) and resync "
                         "ReleaseGeneSymbolGene to single-hop matching (#1669)",
                    test=_has_alias_derived_release_gene_matches),
]
```

Dependency: `("genes", "0088_one_off_stamp_existing_pfam_domains_imported")`.

---

## 6. Tests — `genes/tests/test_gene_matching.py`, new `TestReleaseGeneMatcher`

Build one small release in `setUpTestData`. `annotation/tests/test_data_fake_genes.py` has
`_create_fake_gene_version(genome_build, gene_id, symbol, consortium)`; the release needs a
`GeneAnnotationImport`, a `GeneAnnotationRelease(version, annotation_consortium, genome_build,
gene_annotation_import)` and a `ReleaseGeneVersion(release, gene_version)` per gene — a local helper
`_release_gene(symbol, gene_id)` wrapping those keeps the fixture readable.

Release genes: `PDCD2`=`5134`, `ANOS1`=`3730`, `OLDNAME`=`111` (a gene the release still carries
under its previous symbol). Aliases (all HGNC): `RP8→MT-TS2`, `RP8→PDCD2`, `KAL1→ANOS1`,
`OLDNAME→NEWNAME`. Also create `GeneSymbol` rows for `MT-TS2`, `KAL1`, `NEWNAME`, `RP8`.

Call `gm = ReleaseGeneMatcher(release)` and assert on
`gm._get_gene_id_and_match_info_for_symbol([...])`:

| test | query | expected |
|---|---|---|
| `test_direct_symbol_wins` | `PDCD2` | `[("5134", None)]` |
| `test_forward_alias_hop` | `KAL1` | `[("3730", "KAL1 is an alias for ANOS1 (HGNC)")]` |
| `test_backward_alias_hop` | `NEWNAME` | `[("111", "OLDNAME is an alias for NEWNAME (HGNC)")]` |
| `test_shared_alias_does_not_bridge` | `MT-TS2` | nothing (`RP8` is absent from the release, so the walk stops) |
| `test_alias_matches_only_release_target` | `RP8` | `[("5134", "RP8 is an alias for PDCD2 (HGNC)")]` |

Plus one command test: create `ReleaseGeneSymbol`s for `PDCD2`, `KAL1`, `MT-TS2`, plant a
`ReleaseGeneSymbolGene(MT-TS2 → 5134, "RP8 is an alias for PDCD2 (HGNC)")` and a direct
`PDCD2 → 5134`, run `call_command("fix_rematch_release_symbols_to_genes")`, and assert the planted
chained row is gone, the direct row remains, and `KAL1 → 3730` was inserted. Run the same with
`--dry-run` and assert nothing changed.

Run with `python3 manage.py test --keepdb genes.tests.test_gene_matching`, then the wider
`genes.tests` and `annotation.tests.test_gene_level_annotation` (it builds `ReleaseGeneSymbolGene`
rows by hand and exercises `genes_for_symbol`).

---

## 7. Changelog

`variantgrid/templates/default_templates/changelog.html`, current release block:

```html
<li>#1669 - Gene list matching follows a single alias hop; symbols sharing an alias string (e.g. MT-TS2 / PDCD2 via "RP8") are no longer linked</li>
```

---

## 8. Surfacing the substitution to the user

GeneGrid already shows `Matched <gene>: <match_info>` for alias-matched symbols
(`variantgrid/static_files/default_static/js/gene_grid.js:384`). Showing the same on the gene list
page is #1668 item 4 and stays with that issue.

---

## 9. Immediate prod workaround (before deploy)

```sql
DELETE FROM genes_releasegenesymbolgene rgsg
USING genes_releasegenesymbol rgs
WHERE rgs.id = rgsg.release_gene_symbol_id
  AND rgs.gene_symbol_id = 'MT-TS2'
  AND rgsg.match_info LIKE '%PDCD2%';
```

The full cleanup is §4 run via the §5 migration task.

---

## Files

| file | change |
|---|---|
| `genes/gene_matching.py` | `aliases_dict` rewritten to single hop; `_aliases` removed |
| `genes/management/commands/fix_rematch_release_symbols_to_genes.py` | add-only → full resync, `--dry-run` |
| `genes/migrations/0089_one_off_resync_release_gene_symbol_genes.py` | `ManualOperation` for the resync |
| `genes/tests/test_gene_matching.py` | `TestReleaseGeneMatcher` |
| `variantgrid/templates/default_templates/changelog.html` | entry |
