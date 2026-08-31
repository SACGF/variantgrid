# #1361 — Duplicate VariantAlleles from the ClinGen skip path

Written by Claude Fable 5 (claude-fable-5), 2026-08-31

[#1361](https://github.com/SACGF/variantgrid/issues/1361): tagging `NC_000007.14:g.140753336_140763336del`
(a 10,001 bp deletion) leaves the variant with 4 `Allele`s and the tag never gets an allele, so it
stops showing in other builds.

Three changes, in dependency order:

1. `populate_clingen_alleles_for_variants` creates a skip-path `VariantAllele` **only when the variant has none**.
2. `VariantAllele.unique_together` becomes `("variant", "genome_build")`, with a data migration that removes the
   dupes first.
3. Tags whose `allele` was left `None` by the crashed task get backfilled.

---

## 1. Root cause (verified against current code)

`ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE = 10_000`, so this variant has
`can_have_clingen_allele == False` and takes the "skip ClinGen" branch of
`populate_clingen_alleles_for_variants` (`snpdb/clingen_allele.py:49-185`).

The variant loop at `snpdb/clingen_allele.py:79-89`:

```python
if variant_id not in variant_ids_with_allele or variant_id in allele_missing_clingen_by_variant_id:
    if v.can_have_clingen_allele:
        variant_ids_without_alleles.append(variant_id)   # ClinGen API path - uses existing_va
    else:
        skip_variant_ids_without_alleles.append(variant_id)  # skip path - ALWAYS creates new Allele + VariantAllele
```

- First call: no row exists → skip path creates `Allele()` + `VariantAllele(clingen_error=None)`. Correct.
- Every later call: the row matches the `allele_missing_clingen_by_variant_id` filter at line 63
  (`clingen_error__isnull=True, allele__clingen_allele__isnull=True` — permanently true for a variant that can
  never get a ClinGen allele), so the variant re-enters the `if`. The API branches (lines 120-139) consult
  `existing_va`; the skip branch ignores it and creates another empty `Allele` + `VariantAllele`.
- `bulk_create(ignore_conflicts=True)` at line 185 lets it through because the constraint is
  `("variant", "genome_build", "allele")` and the allele is new each time.

`populate_clingen_alleles_for_variants` is called on every tag creation (`analysis/tasks/variant_tag_tasks.py:45`),
every TagNode update (`analysis/tasks/analysis_update_tasks.py:307`), classification import, liftover, etc. — one
new dupe per call for every variant over 10 kb (or gene-level / non-DEL-DUP-INV symbolic; see
`Variant.clingen_allele_skip_reason`, `snpdb/models/models_variant.py:857`).

The same root shape has a second, smaller leak: passing the same variant twice in `variants` puts it in
`skip_variant_ids_without_alleles` twice (docstring says "more efficient to make variants distinct", but
correctness currently depends on it too).

### Why tags stop showing

- `_liftover_variant_tag` (`analysis/tasks/variant_tag_tasks.py:46`) does `VariantAllele.objects.get(variant, genome_build)`
  → `MultipleObjectsReturned` → task dies before `variant_tag.allele` is assigned and before liftover is queued.
- `VariantTag.get_for_build` (`analysis/models/models_variant_tag.py:95`) finds cross-build tags through
  `allele__in=tags_qs.values_list("allele")`; a tag with `allele=None` only shows in its native build.

---

## 2. Stop the leak — `snpdb/clingen_allele.py`

Restructure the loop so an existing `VariantAllele` is only "work to do" when the variant can actually get a
ClinGen allele, and make the loop idempotent on repeated variants:

```python
seen_variant_ids = set()
for v in variants:
    variant_id = v.pk
    if variant_id in seen_variant_ids:
        continue
    seen_variant_ids.add(variant_id)

    has_allele = variant_id in variant_ids_with_allele
    if has_allele and variant_id not in allele_missing_clingen_by_variant_id:
        num_existing_records += 1
    elif v.can_have_clingen_allele:
        variant_ids_without_alleles.append(variant_id)
        variant_hgvs.append(hgvs_matcher.variant_to_g_hgvs(v))
    elif has_allele:
        num_existing_records += 1   # linked, and can never have ClinGen - nothing to do
    else:
        skip_variant_ids_without_alleles.append(variant_id)
```

Keep the rest of the function as is. Log line at 92 stays meaningful (`num_existing_records` now counts
skip-path variants that are already linked).

### Race hardening in the single-variant path

With the tighter constraint from §3, the `VariantAllele.objects.create(...)` calls in
`get_variant_allele_for_variant` (`snpdb/clingen_allele.py:211`) and the error branch of `variant_allele_clingen`
(line 263) raise `IntegrityError` when two celery tasks race through the "no row yet" window. Wrap each in
`transaction.atomic()` and on `IntegrityError` re-fetch the existing row with
`VariantAllele.objects.get(variant=variant, genome_build=genome_build)` (the just-created empty `Allele` becomes an
orphan — harmless and invisible, nothing joins through it). The bulk paths already use `ignore_conflicts=True` and
need nothing.

---

## 3. Tighten the constraint — `snpdb/models/models_variant.py:996`

```python
class Meta:
    unique_together = ("variant", "genome_build")
```

Rewrite the class docstring: one `Allele` per variant per build; the same `Variant` may legitimately link through
several `VariantAllele` rows across builds that share a contig (MT, some scaffolds — the comment at line 985-986
already says this); many variants may map to one allele (normalisation differences) and that stays allowed.
Delete the paragraph about the 3-column constraint.

### Why this is safe

Every writer already assumes one row per (variant, build):

| Site | Behaviour today | After |
|---|---|---|
| `link_allele_to_existing_variants` (`snpdb/clingen_allele.py:433`) | `get_or_create(variant, genome_build)` then `merge` | unchanged |
| `get_variant_allele_for_variant` (line 195) | `.order_by("pk").first()` with a comment apologising for the race | `.get()` is now correct; keep `first()` or switch — either works |
| `BulkAlleleLinkingVCFProcessor` (`upload/vcf/bulk_allele_linking_vcf_processor.py:50-85`) | pre-queries existing rows, merges alleles, marks liftover `SKIPPED` | unchanged; `ignore_conflicts` now a true no-op guard |
| `_run_liftover_using_same_contig` (`snpdb/liftover.py:263-286`) | `bulk_create(ignore_conflicts=True)` would add a *second* link if a different allele already owns the variant | see below |
| `Allele.merge` (`snpdb/models/models_variant.py:182-191`) | `except IntegrityError: delete` | fires when `self` already has a row for that variant/build — exactly the case the comment describes |
| `fake_variant_tags._create_alleles` | fresh variants, fresh alleles | unchanged |

`_run_liftover_using_same_contig`: with the new constraint a conflicting row is silently skipped but its
`AlleleLiftover` is still written as `SUCCESS`. Bring it in line with the bulk processor: pre-query
`VariantAllele.objects.filter(variant__in=[v for _, v in av_tuples], genome_build=liftover.genome_build)` into a
`{variant_id: allele_id}` dict; for a tuple whose variant is already linked to a *different* allele, merge into the
lower-pk allele (same rule as `BulkAlleleLinkingVCFProcessor.merge_alleles`) and record the `AlleleLiftover` as
`SKIPPED` with an error message; otherwise proceed as now.

### Migration `snpdb/migrations/0224_variantallele_unique_variant_genome_build.py`

Two operations, in order:

**(a) `RunPython(_dedupe_variant_alleles, noop)`** — runs before the constraint so `AlterUniqueTogether` succeeds.

```
groups = VariantAllele.objects.values("variant", "genome_build").annotate(n=Count("pk")).filter(n__gt=1)
```

For each group, fetch its rows ordered by pk and pick the **keeper**:

1. the row whose allele has `clingen_allele_id`; else
2. the row whose allele is referenced by anything (list below); else
3. the lowest pk.

Every other row in the group is a **loser**. A loser is *orphan* when its allele has no `VariantAllele` rows
outside this group and no references. Orphan losers: delete the `VariantAllele` row, then the `Allele`. That is
the entire population produced by this bug (empty allele, nothing attached, the tag task crashed before it could
point at one), and the fix that runs automatically.

Any loser that is *not* orphan gets collected, and after the loop the migration raises with one line per group —
`variant_id / build / keeper allele / loser alleles and what references them` — so it can be resolved by hand
(`Allele.merge` in a shell, then repoint the PROTECT FKs) and `migrate` re-run. Silent auto-merging of alleles
that carry classifications is worth stopping the upgrade for; the dev-DB check from the issue's second comment
found none.

Reference tables to check on a loser allele (all via `apps.get_model`):

- `snpdb.VariantAllele` (other builds), `snpdb.AlleleLiftover`, `snpdb.AlleleMergeLog` (either side)
- `classification.Classification.allele`, `classification.ClinicalContext.allele`,
  `classification.ImportedAlleleInfo.allele`, `classification.ClinVarExport.allele`,
  `classification.AlleleGrouping.allele`
- `analysis.VariantTag.allele`
- `annotation.ClinVarRecordCollection.allele`
- `Allele.flag_collection_id` non-null

**(b) `AlterUniqueTogether(name="variantallele", unique_together={("variant", "genome_build")})`**

Test the dedupe against the dev DB first by running the issue's diagnostic script before and after
`migrate snpdb 0224`.

---

## 4. Backfill tags left without an allele

`_liftover_variant_tag` needs no structural change once §2/§3 land — `.get()` is now the correct call.

Existing tags with `allele IS NULL` on deployments that hit the bug never got liftover. Add
`analysis/management/commands/one_off_backfill_variant_tag_alleles.py`:

```
for variant_tag in VariantTag.objects.filter(allele__isnull=True).select_related("variant", "genome_build"):
    _liftover_variant_tag(variant_tag)
```

(`_liftover_variant_tag` already does populate → assign allele → `create_liftover_pipelines`; the null-allele
population is small so per-tag calls are fine. Make `_liftover_variant_tag` public — `liftover_variant_tag` — since
it now has a second caller.)

Register it with `analysis/migrations/0125_one_off_backfill_variant_tag_alleles.py` containing a
`ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_backfill_variant_tag_alleles"]), test=...)`
where `test` returns `VariantTag.objects.filter(allele__isnull=True).exists()`. Depend on `snpdb 0224` so the
command runs against a deduped table.

---

## 5. Tests

`snpdb/tests/test_clingen_allele.py` (fixture already builds GRCh37 annotation for `HGVSMatcher`):

- **Skip path is idempotent** — create a 10,001 bp deletion with
  `slowly_create_test_variant("7", 140753336, "A" * 10_002, "A", grch37)` (the helper runs
  `as_internal_symbolic`, which converts any ref longer than `VARIANT_SYMBOLIC_ALT_SIZE` to `<DEL>` with an
  `svlen`, so `can_have_clingen_allele` is `False`; assert that first). Call
  `populate_clingen_alleles_for_variants` three times, assert exactly one `VariantAllele` and one `Allele` for the
  variant. This is the regression test for the issue.
- **Repeated input variants** — pass `[variant, variant]` once, assert one row.
- **Skip path then API path share the row** — a variant with a row from the skip path is left alone; a
  `can_have_clingen_allele` variant with an existing empty row gets its allele's `clingen_allele` set and gains
  no second row (covers the `elif has_allele` branch ordering).

`snpdb/tests/test_liftover.py`:

- **Same-contig liftover onto an already-linked variant** — variant already linked to allele A, lift allele B
  onto it via `_run_liftover_using_same_contig`; assert A and B merged into the lower pk, one `VariantAllele`,
  `AlleleLiftover.status == SKIPPED`.

Skip a test that the constraint itself raises `IntegrityError` — that is Django's behaviour.

---

## 6. Changelog

Add a line to `variantgrid/templates/default_templates/changelog.html` under the next release: "Fixed variants
over 10 kb (and other variants ClinGen cannot register) accumulating a duplicate Allele every time they were
tagged, lifted over or imported; tags on those variants now appear in every build (#1361)".

---

## 7. Out of scope, noted for later

- `Allele.merge` moves flags, clinical contexts and `VariantAllele` rows but leaves `Classification.allele`,
  `VariantTag.allele`, `ImportedAlleleInfo.allele`, `ClinVarExport.allele` and `AlleleGrouping` pointing at the
  old allele. Separate issue.
- `populate_clingen_alleles_for_variants` filters on `clingen_error__isnull=True`, so a row that recorded a
  ClinGen *server* error is never retried in bulk, although `VariantAllele.needs_clingen_call` says it should
  be. Separate issue.
