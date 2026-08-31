# #1361 — Duplicate `VariantAllele`s per (variant, build); tags not showing

[#1361](https://github.com/SACGF/variantgrid/issues/1361): tagging `NC_000007.14:g.140753336_140763336del`
left the variant linked to four `Allele`s in GRCh38, and the tag stopped showing.

The root-cause comment on the issue (June 2026) is correct. This plan re-verifies it against the current
code, adds a second duplicate-writer found in the dev DB, and lays out the fix. Reproduced locally (rolled
back) on `vg-test2` with a 36,005 bp GRCh38 deletion:

```
after call 1: VariantAllele rows = 1, alleles = [53356]
after call 2: VariantAllele rows = 2, alleles = [53356, 53357]
after call 3: VariantAllele rows = 3, alleles = [53356, 53357, 53358]
```

The dev DB already holds **66** duplicate `(variant, genome_build)` groups (all pairs) — see §2.4.

*Second pass (same day) added:* the three models that FK onto `VariantAllele` (§2.4), the missing savepoint
in `Allele.merge()` (§3.5), the fact that moving `Classification.allele` strands clinical contexts and
groupings (§3.5, §3.9), safer `ClinicalContext` handling in the migration (§3.6), exact cross-app migration
dependencies, and a test-fixture audit (§4).

---

## 1. Root cause

### 1.1 The leak — `populate_clingen_alleles_for_variants` (`snpdb/clingen_allele.py:49`)

```python
va_qs = VariantAllele.objects.filter(genome_build=genome_build, variant__in=variants)
variant_ids_with_allele = set(va_qs.values_list("variant_id", flat=True))
allele_missing_clingen_by_variant_id = {}
for va in va_qs.filter(clingen_error__isnull=True, allele__clingen_allele__isnull=True):   # line 63
    allele_missing_clingen_by_variant_id[va.variant_id] = va

for v in variants:                                                                          # line 79
    variant_id = v.pk
    if variant_id not in variant_ids_with_allele or variant_id in allele_missing_clingen_by_variant_id:
        if v.can_have_clingen_allele:
            variant_ids_without_alleles.append(variant_id)          # -> ClinGen API, handles existing_va
        else:
            skip_variant_ids_without_alleles.append(variant_id)     # -> NEW empty Allele + NEW VariantAllele
    else:
        num_existing_records += 1
```

A variant that can *never* get a ClinGen allele (`Variant.clingen_allele_skip_reason()`,
`snpdb/models/models_variant.py:843`: gene-level, un-HGVS-able symbolic, or
`> ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE = 10_000` bp) gets an empty `Allele` on the first call and
then matches the line-63 filter forever (`clingen_error` null, no `clingen_allele`). Every later call routes it
back to the skip branch, which — unlike the API-error and API-success branches, which check `existing_va` —
ignores the existing row and bulk-creates a fresh `Allele` plus `VariantAllele` (lines 96–101, 168–184).

The issue's variant is a 10,001 bp deletion: one over the limit.

### 1.2 Why the constraint doesn't catch it

`VariantAllele.Meta.unique_together = ("variant", "genome_build", "allele")` (`models_variant.py:978`). The
`allele` is brand new each time, so `bulk_create(..., ignore_conflicts=True)` never conflicts.

### 1.3 Why it accumulates

`populate_clingen_alleles_for_variants` is re-run for the same variants from many places, several of them
per-event: `variant_tag_created_task` (`analysis/tasks/variant_tag_tasks.py:45`),
`populate_clingen_alleles_from_analysis_node` (`analysis/tasks/analysis_update_tasks.py:307` — every TagNode
update), `liftover_run_complete_handler` (`classification/signals/classification_liftover.py:20`),
classification imports, `liftover_alleles`, `BulkClinGenAlleleVCFProcessor`, and two management commands.
One duplicate per call.

### 1.4 Why the tag disappears

- `_liftover_variant_tag` (`variant_tag_tasks.py:46`) does
  `VariantAllele.objects.get(variant=..., genome_build=...)`. With duplicates this raises
  `MultipleObjectsReturned`; the celery task dies and `VariantTag.allele` stays `None` for any tag created
  after the first duplicate appeared.
- `VariantTag.get_for_build` (`analysis/models/models_variant_tag.py:95`) matches tags via
  `allele__in=tags_qs.values_list("allele")`. A tag with `allele=None` is only visible in its own build via
  the `q_native_build` fallback, and never in the other build.
- `Variant.allele` (`models_variant.py:887`) is `variantallele_set.first()` with no ordering, so
  `ImportedAlleleInfo.update_and_save` (`classification_variant_info_models.py:1062`) and anything else
  going variant → allele may pick a different duplicate than the tag did.

---

## 2. Design

### 2.1 Tighten the constraint to `("variant", "genome_build")`

Nothing in the codebase wants one variant in one build linked to two alleles:

- `link_allele_to_existing_variants` (`clingen_allele.py:433`) already does
  `get_or_create(variant=, genome_build=)` and merges when the allele differs.
- `BulkAlleleLinkingVCFProcessor.batch_handle_variant_ids` (`upload/vcf/bulk_allele_linking_vcf_processor.py:50`)
  looks up `variants_with_existing_allele` keyed on `(variant, build)`, merges, and comments
  *"Only one VariantAllele allowed"*.
- `get_variant_allele_for_variant` (`clingen_allele.py:194`) apologises for the race the wide constraint permits.

The docstring's "multiple variants → same allele" case is the *other* direction and stays legal. The
"same variant, several builds" case (MT / shared contigs, `models_variant.py:966`) is distinct
`genome_build` values and stays legal.

Once the narrow constraint exists, every `bulk_create(ignore_conflicts=True)` becomes genuine dedupe and the
concurrency race collapses to "one writer wins" (see §2.3 for the loser's orphaned `Allele`).

### 2.2 Every `VariantAllele` writer, and what changes

| Writer | Today | After the constraint | Action |
|---|---|---|---|
| `populate_clingen_alleles_for_variants` skip branch (`clingen_allele.py:96`) | leaks (§1.1) | still leaks an orphan `Allele` per call; the VA insert is now silently ignored | **fix the branch** (§3.1) |
| same function, API-error / API-success new-VA branches (`:127`, `:138`) | dupe only in a race | race loser's VA ignored; error-branch loser leaves an orphan empty `Allele` | delete unlinked new alleles after the bulk_create (§3.1) |
| `get_variant_allele_for_variant` no-ClinGen branch (`:211`) | dupe only in a race | race loser raises `IntegrityError` | catch, drop the just-made `Allele`, re-fetch (§3.2) |
| `variant_allele_clingen` error branch (`:263`) | dupe only in a race | as above | same helper (§3.2) |
| `link_allele_to_existing_variants` (`:433`) | `get_or_create` + merge | unchanged | none |
| `_run_liftover_using_same_contig` (`snpdb/liftover.py:263`) | if the dest-build variant already has its own allele, adds a second link | insert silently dropped, `AlleleLiftover` still marked `SUCCESS`, alleles never merged | **merge like the bulk processor does** (§3.3) |
| `BulkAlleleLinkingVCFProcessor` (`bulk_allele_linking_vcf_processor.py:97`) | detects existing, merges | unchanged | none |
| `analysis/fake_variant_tags.py:_create_alleles` (`:243`, dev only) | blindly makes a new `Allele`+VA per picked variant — this is 57 of the 66 dev dupes | `bulk_create` (no `ignore_conflicts`) raises | **reuse the existing allele** (§3.4) |
| `upload/tasks/import_variant_tags_task.py:184` | reads VAs into `dict(variant_id → allele_id)` — with dupes an arbitrary one wins | one row per key, deterministic | none |

### 2.3 Orphan `Allele`s

`Allele.merge()` (`models_variant.py:146`) already leaves the merged-away allele in place with no
`VariantAllele` (dev DB: 2 such orphans). The race losers in §2.2 add a few more empty ones. They are
harmless — nothing navigates to an allele except through the FKs — so the plan tolerates them rather than
adding cleanup machinery. `populate_clingen_alleles_for_variants` can cheaply delete its own losers (§3.1)
because it knows which alleles it just created.

### 2.4 Existing duplicates — what the dedupe must handle

Profile of the 66 dev-DB groups (all size 2, all `origin='D'`):

| Flavour | Groups | Shape |
|---|---|---|
| Skip-path leak (§1.1) | 9 | variant can't have ClinGen; both alleles empty; surplus allele referenced only by `AlleleLiftover` rows |
| Fake tag generator | 57 | old allele has `clingen_allele` + VAs in the other build (`tool=CA`); new allele (`tool=SC`) carries the `VariantTag`s and failed `AlleleLiftover`s |
| Either | 1 | surplus allele referenced by a `Classification` + `ImportedAlleleInfo` |

So the dedupe has to (a) pick a keeper, (b) re-point *every* FK to `Allele`, not just the ones
`Allele.merge()` moves (it moves `VariantAllele`, `clingen_allele`, flags and `ClinicalContext` only —
`VariantTag.allele`, `Classification.allele`, `ImportedAlleleInfo.allele`, `AlleleLiftover.allele`,
`ClinVarAllele`, `ClinVarRecordCollection.allele`, `AlleleGrouping` are left behind). `VariantTag.allele`
and `Classification.allele` are `on_delete=PROTECT`, so nothing can be deleted until they're moved.

FKs to `Allele` outside snpdb (from `grep ForeignKey(Allele`):

| Model | Field | on_delete | Uniqueness on allele |
|---|---|---|---|
| `analysis.VariantTag` | `allele` | PROTECT | – |
| `classification.Classification` | `allele` | PROTECT | – |
| `classification.ImportedAlleleInfo` | `allele` | SET_NULL | – |
| `classification.ClinicalContext` | `allele` | CASCADE | `(allele, allele_origin_bucket, name)` |
| `classification.ClinVarAllele` | `allele` | CASCADE | none in the DB; derived from classifications |
| `classification.AlleleGrouping` | `allele` | CASCADE | OneToOne |
| `annotation.ClinVarRecordCollection` | `allele` | SET_NULL | – |
| `snpdb.AlleleLiftover` | `allele` | CASCADE | – |
| `snpdb.AlleleMergeLog` | `old_allele`/`new_allele` | CASCADE | – (history; leave alone) |
| `snpdb.VariantAllele` | `allele` | CASCADE | the thing being fixed |

And the models that FK onto `VariantAllele` itself — these cascade when the migration deletes a surplus row:

| Model | Field | on_delete | Notes |
|---|---|---|---|
| `snpdb.AlleleLiftover` | `variant_allele` | OneToOne, CASCADE | never assigned anywhere in code (0 of 204,260 rows set in dev); null it out before deleting a surplus VA anyway so a liftover record can't vanish |
| `snpdb.VariantAlleleSource` | `variant_allele` | CASCADE | no creators anywhere in the codebase — dead model |
| `snpdb.VariantAlleleCollectionRecord` | `variant_allele` | CASCADE | as above |

No signal receivers listen on `VariantAllele` or `Allele` save/delete (only `preview_extra_signal`), so
neither the migration (historical models never fire signals anyway) nor `Allele.merge()` has hidden
side-effects to worry about.

`Classification.clinical_context` is `SET_NULL` and `Classification.allele` is derived from
`ImportedAlleleInfo.allele` on save (`classification.py:943`). Clinical contexts are re-homed lazily by
`_update_clinical_context` (`clinical_context_utils.py:57`) and groupings by
`ClassificationGrouping.assign_grouping_for_classification`, both driven by classification publish/withdraw/
allele-info signals — so a bulk `.update(allele=keeper)` on classifications leaves their
`clinical_context` and `ClassificationGroupingEntry` pointing at the *old* allele until something re-saves
each one. §3.9 handles that.

### 2.5 Where the dedupe runs

The upgrader (`scripts/migrator/migrator.py:307`) runs `manage.py migrate` in the standard pass and only
auto-runs `ManualOperation` tasks afterwards. A `ManualOperation` therefore can't precede an
`AlterUniqueTogether` in the same release, and a migration that imports live models breaks fresh-DB builds
as soon as a later migration adds a column. So the dedupe is a `RunPython` on **historical models** in the
same migration as the constraint change, immediately before it. Cross-app dependencies in an snpdb migration
have precedent (`snpdb/migrations/0208_one_off_variant_tag_cleanup.py` depends on `analysis`).

Volume is small: the leak only affects variants that are both never-ClinGen-able *and* repeatedly pushed
through `populate_clingen_alleles_for_variants` (tagged large SVs, gene-level fusion variants). The
grouping query is a single index scan of `snpdb_variantallele`.

---

## 3. Implementation

### 3.1 Stop the leak — `snpdb/clingen_allele.py`

Rewrite the classification loop so a variant that can never have a ClinGen allele counts as *existing* when it
already has a `VariantAllele`:

```python
for v in variants:
    variant_id = v.pk
    has_allele = variant_id in variant_ids_with_allele
    if has_allele and variant_id not in allele_missing_clingen_by_variant_id:
        num_existing_records += 1
    elif v.can_have_clingen_allele:
        variant_ids_without_alleles.append(variant_id)
        variant_hgvs.append(hgvs_matcher.variant_to_g_hgvs(v))
    elif has_allele:
        num_existing_records += 1  # Can never get ClinGen - the VariantAllele it has is all it will ever have
    else:
        skip_variant_ids_without_alleles.append(variant_id)
```

Then, after the final `VariantAllele.objects.bulk_create(variant_allele_list, ignore_conflicts=True)`
(line 184), drop any empty allele this call created that didn't end up linked — this is the race loser
under the new constraint:

```python
new_empty_allele_ids = [a.pk for a in itertools.chain(reference_alleles, allele_no_clingen_list)]
Allele.objects.filter(pk__in=new_empty_allele_ids, variantallele__isnull=True).delete()
```

(`allele_no_clingen_list` is consumed by `.pop()` at line 161 — capture the ids before that loop, or keep
a separate list of created alleles.) Alleles with a `clingen_allele` are deduped by the OneToOne already and
need no cleanup.

The docstring says callers *should* pass distinct variants but don't have to. A repeated variant currently
lands in the skip branch twice (two alleles, two VAs). Track a `seen` set of variant ids in the loop so a
repeat is ignored — that also avoids sending the same HGVS to ClinGen twice.

### 3.2 Race-safe single creates — `snpdb/clingen_allele.py`

`get_variant_allele_for_variant` (`:211`) and `variant_allele_clingen` (`:263`) both do
`Allele.objects.create()` then `VariantAllele.objects.create(...)`. Under the new constraint a concurrent
creator makes the second call raise `IntegrityError`. Extract one helper used by both:

```python
def _create_variant_allele_with_new_allele(variant, genome_build, **va_kwargs) -> VariantAllele:
    """ Another task may link the variant between our Allele and VariantAllele inserts - on a clash use
        theirs and discard the Allele we just made """
    allele = Allele.objects.create()
    try:
        with transaction.atomic():
            return VariantAllele.objects.create(variant_id=variant.pk, genome_build=genome_build,
                                                allele=allele,
                                                origin=AlleleOrigin.variant_origin(variant, allele, genome_build),
                                                **va_kwargs)
    except IntegrityError:
        allele.delete()
        return VariantAllele.objects.get(variant=variant, genome_build=genome_build)
```

`variant_allele_clingen` passes `clingen_error=api_response`. With uniqueness guaranteed,
`get_variant_allele_for_variant`'s `.order_by("pk").first()` and its "very rare race condition" comment
(`:194`) can become a plain `.filter(...).first()` — leave the ordering, drop the comment.

### 3.3 Same-contig liftover must merge — `snpdb/liftover.py:263`

`_liftover_using_existing_contig` (`:289`) only checks whether the *allele* already has a link in the
destination build, not whether the destination-build *variant* already has an allele of its own (e.g. an MT
variant imported directly into GRCh38 before the GRCh37 copy was lifted over). Mirror
`BulkAlleleLinkingVCFProcessor`:

```python
existing_allele_id_by_variant_id = dict(
    VariantAllele.objects.filter(genome_build=liftover.genome_build,
                                 variant__in=[v for _, v in av_tuples]).values_list("variant_id", "allele_id"))
for allele, variant in av_tuples:
    if existing_allele_id := existing_allele_id_by_variant_id.get(variant.pk):
        if existing_allele_id != allele.pk:
            keep, other = sorted([existing_allele_id, allele.pk])   # lowest pk wins, as merge_alleles does
            Allele.objects.get(pk=keep).merge(AlleleConversionTool.SAME_CONTIG, Allele.objects.get(pk=other))
        status = ProcessingStatus.SKIPPED   # nothing inserted - matches the bulk processor's convention
    else:
        variant_alleles.append(VariantAllele(...))
        status = ProcessingStatus.SUCCESS
    allele_liftovers.append(AlleleLiftover(allele=allele, liftover=liftover, status=status))
```

### 3.4 Fake tag generator — `analysis/fake_variant_tags.py:243`

`_create_alleles` should look up existing `VariantAllele`s for the picked variants first, assign
`fake_variant.allele_id` from those, and only create `Allele`+`VariantAllele` for variants that have none.
This is the writer behind 57 of the dev DB's 66 groups.

### 3.5 Extend `Allele.merge()` to move the FKs it currently strands — `snpdb/models/models_variant.py:146`

Production merges (`link_allele_to_existing_variants`, liftover) leave `VariantTag.allele` pointing at an
allele with no `VariantAllele`, which is the same "tag not visible in the other build" symptom as this
issue. Add, inside the `if can_merge:` block, using reverse accessors (no cross-app imports needed):

```python
other_allele.varianttag_set.update(allele=self)
other_allele.classification_set.update(allele=self)
other_allele.importedalleleinfo_set.update(allele=self)
other_allele.alleleliftover_set.update(allele=self)
other_allele.clinvarrecordcollection_set.update(allele=self)
```

Leave `ClinVarAllele` and `AlleleGrouping` alone: both are rebuilt from classifications
(`AlleleGrouping.objects.get_or_create(allele=allele)`, `classification_grouping.py:325`) and both carry
uniqueness that a blind update could violate. Leave `AlleleMergeLog` as history.

The migration in §3.6 does the same moves against historical models; keeping the two lists identical is
the point of writing both in this plan.

**Savepoint.** The existing `VariantAllele` move (`models_variant.py:182-191`) does `try: va.save()
except IntegrityError: va.delete()` with no `transaction.atomic()` around the save. Under the new
constraint that collision is the *normal* case whenever both alleles already link the same variant in the
same build (exactly the §3.3 scenario), and an uncaught-at-DB-level `IntegrityError` inside an outer
transaction (tests, `ATOMIC_REQUESTS`, a celery task wrapped in `atomic`) leaves the connection in
"current transaction is aborted" state for everything after it. Wrap the save:

```python
try:
    with transaction.atomic():
        va.save()
except IntegrityError:
    ...
```

**Classifications moved by merge need re-homing.** After `classification_set.update(allele=self)` the
moved classifications still point at a `ClinicalContext` on `other_allele` and sit in `other_allele`'s
`AlleleGrouping`. snpdb can't import classification code, so add a signal next to
`liftover_run_complete_signal` in `snpdb/models`:

```python
allele_merged_signal = django.dispatch.Signal()   # kwargs: old_allele, new_allele
```

sent at the end of a successful `merge()`. A receiver in `classification/signals/` (new
`classification_hooks_allele_merge.py`) does, for
`Classification.objects.filter(allele=new_allele).exclude(clinical_context__allele=new_allele)` plus any
whose grouping entry is under another allele:

```python
update_clinical_contexts(classifications, force_recalc_text="allele merged")
for classification in classifications:
    ClassificationGrouping.assign_grouping_for_classification(classification)
ClassificationGrouping.update_all_dirty()
```

`update_clinical_contexts` (`clinical_context_utils.py:126`) already knows how to move a classification
to the matching context on its new allele and recalc discordance on both sides. The same code is what
§3.9's command runs, so put the "re-home these classifications" body in one function in
`clinical_context_utils` (or the command module) and call it from both.

### 3.6 Migration — `snpdb/migrations/0224_variantallele_unique_variant_build.py`

Change `VariantAllele.Meta.unique_together` first, run `makemigrations --check --dry-run` to confirm the
only pending change is the `AlterUniqueTogether`, then write the migration by hand.

Dependencies (current heads — re-check at implementation time):

```python
dependencies = [
    ("snpdb", "0223_alter_nodecountsettings_built_in_filter"),
    ("analysis", "0124_clinvarnode_and_more"),
    ("classification", "0175_ekey_gene_fusion_options"),
    ("annotation", "0177_one_off_backfill_clinvar_somatic"),
]
```

Operations, in order:

1. `RunPython(_dedupe_variant_alleles, RunPython.noop)`
2. `AlterUniqueTogether(name="variantallele", unique_together={("variant", "genome_build")})`
3. `ManualOperation(task_id=ManualOperation.task_id_manage(["classification_fix_allele_links"]),
   note="Re-home clinical contexts and groupings of classifications moved off duplicate alleles (#1361)",
   test=_has_classifications_on_wrong_allele)` — see §3.9. Placed *after* the `RunPython` so its
   `test` sees the post-dedupe state; on most deployments it registers nothing.

The migration is atomic (Postgres), so a `RuntimeError` from step 1 rolls back everything. Step 2 drops
the 3-column unique index and builds a 2-column one under an `ACCESS EXCLUSIVE` lock on
`snpdb_variantallele` — seconds on a 100k-row dev DB, but on SA Path-scale tables it's a visible pause
during `migrate`. Normal for this project's deploys; just don't be surprised.

`_dedupe_variant_alleles(apps, schema_editor)`:

```python
VariantAllele = apps.get_model("snpdb", "VariantAllele")
Allele = apps.get_model("snpdb", "Allele")
VariantTag = apps.get_model("analysis", "VariantTag")
Classification = apps.get_model("classification", "Classification")
ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
ClinicalContext = apps.get_model("classification", "ClinicalContext")
ClinVarRecordCollection = apps.get_model("annotation", "ClinVarRecordCollection")
AlleleLiftover = apps.get_model("snpdb", "AlleleLiftover")

dupe_keys = (VariantAllele.objects.values("variant_id", "genome_build_id")
             .annotate(n=Count("pk")).filter(n__gt=1))
unresolvable = []
for key in dupe_keys:
    vas = list(VariantAllele.objects.filter(variant_id=key["variant_id"], genome_build_id=key["genome_build_id"])
               .select_related("allele").order_by("pk"))
    keeper_va = _pick_keeper(vas)               # see below
    keeper = keeper_va.allele
    for va in vas:
        if va.pk == keeper_va.pk:
            continue
        other = va.allele
        if other.clingen_allele_id and keeper.clingen_allele_id and other.clingen_allele_id != keeper.clingen_allele_id:
            unresolvable.append((key, keeper.pk, other.pk))
            continue
        if other.flag_collection_id and not keeper.flag_collection_id:
            keeper.flag_collection_id = other.flag_collection_id
            keeper.save(update_fields=["flag_collection_id"])
        VariantTag.objects.filter(allele=other).update(allele=keeper)
        Classification.objects.filter(allele=other).update(allele=keeper)
        ImportedAlleleInfo.objects.filter(allele=other).update(allele=keeper)
        AlleleLiftover.objects.filter(allele=other).update(allele=keeper)
        ClinVarRecordCollection.objects.filter(allele=other).update(allele=keeper)
        # ClinicalContext is unique on (allele, allele_origin_bucket, name). Move the ones that don't
        # collide; leave colliding ones on the surplus allele, exactly as Allele.merge() does - §3.9
        # re-homes their classifications afterwards. Never delete a context: discordance history hangs off it.
        keeper_cc_keys = set(ClinicalContext.objects.filter(allele=keeper).values_list("allele_origin_bucket", "name"))
        for cc in ClinicalContext.objects.filter(allele=other):
            if (cc.allele_origin_bucket, cc.name) not in keeper_cc_keys:
                cc.allele = keeper
                cc.save(update_fields=["allele"])
        # Other-build links on the surplus allele move across too (mirrors Allele.merge); a link the keeper
        # already has in that build is dropped rather than moved
        keeper_builds = set(VariantAllele.objects.filter(allele=keeper).values_list("genome_build_id", flat=True))
        for other_va in VariantAllele.objects.filter(allele=other).exclude(pk=va.pk):
            if other_va.genome_build_id in keeper_builds:
                AlleleLiftover.objects.filter(variant_allele=other_va).update(variant_allele=None)
                other_va.delete()
            else:
                other_va.allele = keeper
                other_va.save(update_fields=["allele"])
        AlleleLiftover.objects.filter(variant_allele=va).update(variant_allele=None)
        va.delete()
if unresolvable:
    raise RuntimeError(f"VariantAllele duplicates with two different ClinGen alleles need a human: {unresolvable}")
```

`_pick_keeper(vas)`: first VA whose allele has `clingen_allele_id`; else the one whose allele is referenced
by a `Classification` or `VariantTag`; else lowest pk. Lowest pk matches what `Variant.allele` and
`get_variant_allele_for_variant` have been returning, so the tag created in the issue keeps its allele.
Because a ClinGen-bearing allele is always chosen as keeper when one exists, the loop never needs to move
a `clingen_allele` across — the only ClinGen case left is the "two different CA ids" one, which raises.

Note the moved-VA handling above is what §3.5's savepoint fix does in live code; the migration is written
without `try/except IntegrityError` because a failed statement inside an atomic migration is fatal.

The "both have ClinGen" case never occurred in the dev DB (max one per group) and shouldn't be possible —
one variant can't register two CA ids — so raising is right: it stops the migration with the ids listed
rather than guessing.

Moving `VariantAllele` rows of the surplus allele to the keeper can itself collide with the new key if the
keeper already has a row in that other build; use `ignore_conflicts`-style handling by iterating and
deleting on collision, as `Allele.merge()` does (`models_variant.py:182`).

The surplus `Allele` rows stay behind as orphans (§2.3). If you'd rather delete them: after the loop,
`Allele.objects.filter(pk__in=surplus_ids, variantallele__isnull=True, varianttag__isnull=True,
classification__isnull=True, ...).delete()` — but `AlleleMergeLog` rows referencing them cascade away, and
nothing needs the space, so the plan leaves them. Any `ClinicalContext`s left on an orphan allele because
they collided with the keeper's are re-homed by §3.9; `ClinVarAllele`/`AlleleGrouping` rows on orphans
are inert.

Write a log line per group (`variant_id`, build, keeper allele, surplus alleles, counts moved) — this is the
only record of what the migration did besides `AlleleMergeLog`, which the migration does *not* write (it's
a merge in effect, so writing one row per surplus allele with `allele_linking_tool=""` and a
"#1361 dedupe" message is a reasonable optional extra).

Update the `VariantAllele` docstring to describe the new invariant ("one allele per variant per build;
several variants in a build may share an allele").

### 3.7 Backfill `VariantTag.allele` — `analysis/migrations/0125_backfill_variant_tag_allele.py`

Depends on `snpdb.0224`. `RunPython`:

```python
VariantTag = apps.get_model("analysis", "VariantTag")
VariantAllele = apps.get_model("snpdb", "VariantAllele")
for tag in VariantTag.objects.filter(allele__isnull=True).iterator(chunk_size=2000):
    va = VariantAllele.objects.filter(variant_id=tag.variant_id, genome_build_id=tag.genome_build_id).first()
    if va:
        tag.allele_id = va.allele_id
        tag.save(update_fields=["allele_id"])
```

Dev DB has 0 such tags; the issue's deployment has at least one. Tags whose variant has no `VariantAllele`
at all (task died before `populate` ran) are left for the next `variant_tag_created_task` /
`fix_liftover_existing_variant_tags` — this migration only repairs what it can prove.

Backfilled tags become visible in their own build immediately (`q_native_build`). The crashed task also
never called `create_liftover_pipelines` for them, so they won't show in the *other* build until their
allele is lifted over. For the variants this bug hits that's mostly moot — >10 kb symbolic variants can't
go through ClinGen and bcftools symbolic liftover is off in the deployments seen
(`LIFTOVER_BCFTOOLS_SYMBOLIC=False` in the dev `AlleleLiftover` errors). If wanted, from a shell after
deploy: `liftover_alleles(Allele.objects.filter(varianttag__in=<backfilled>))`. Not worth a migration.

### 3.8 `_liftover_variant_tag` — `analysis/tasks/variant_tag_tasks.py:46`

Keep the `.get()`. After §3.6 it's correct and fails loudly if the invariant is ever broken again.

### 3.9 `classification_fix_allele_links` management command (new, `classification/management/commands/`)

Re-homes classifications whose derived records disagree with `Classification.allele`:

```python
mismatched = (Classification.objects.filter(allele__isnull=False)
              .filter(Q(clinical_context__isnull=False) & ~Q(clinical_context__allele=F("allele"))
                      | ~Q(classificationgroupingentry__grouping__allele_origin_grouping__allele_grouping__allele=F("allele"))))
```

(check the exact reverse path from `Classification` to `ClassificationGroupingEntry` →
`ClassificationGrouping` → `AlleleOriginGrouping` → `AlleleGrouping.allele` when implementing). For each:
`update_clinical_contexts([...], force_recalc_text="allele merged (#1361)")`,
`ClassificationGrouping.assign_grouping_for_classification(c)`, then one
`ClassificationGrouping.update_all_dirty()`. Print counts; `--dry-run` flag.

Registered by the `ManualOperation` in §3.6 with `test=_has_classifications_on_wrong_allele` (same
clinical-context predicate on historical models — the grouping path is optional in the test). Dev DB:
exactly one classification sits on a duplicate-group allele, so the task would register there. The command
is also the body of the §3.5 signal receiver, so it earns its keep beyond the one-off.

---

## 4. Tests

Add to `snpdb/tests/test_clingen_allele.py` (uses `get_fake_annotation_version` + `MockClinGenAlleleRegistryAPI`):

- **`test_populate_never_clingen_variant_is_idempotent`** — build a GRCh37 `<DEL>` with
  `svlen=-(CLINGEN_ALLELE_MAX_ALLELE_SIZE + 1)` (`VariantCoordinate(chrom=, position=, ref=, alt='<DEL>',
  svlen=)` as in `snpdb/tests/test_variant.py:68`, then `Variant.objects.get_or_create` the way
  `slowly_create_test_variant` does). Call `populate_clingen_alleles_for_variants` three times; assert one
  `VariantAllele` and one `Allele` for the variant. This is the regression test for the issue.
- **`test_populate_creates_one_allele_for_new_never_clingen_variant`** — first call on a fresh variant
  creates exactly one VA with `clingen_error=None` (covers the `else:` branch of the rewritten loop).
- **`test_create_variant_allele_race_reuses_existing`** — pre-create a VA for the variant, call
  `get_variant_allele_for_variant` with `CLINGEN_ALLELE_REGISTRY_LOGIN` unset; assert it returns the
  existing VA and `Allele.objects.count()` is unchanged. (Exercises the helper's happy path; the
  `IntegrityError` branch is only reachable concurrently, so mock `VariantAllele.objects.create` to raise
  once if you want it covered.)

Add to `snpdb/tests/test_allele.py`:

- **`test_merge_moves_variant_tags_and_classification_fks`** — allele B has a `VariantTag` (and an
  `AlleleLiftover`); merge into A; assert they now point at A. Covers §3.5.
- **`test_merge_when_both_alleles_link_same_variant_and_build`** — A and B each have a VA for the same
  variant/build; merge B into A inside `TestCase`'s transaction; assert one VA remains on A and the test
  can still query afterwards. This is the savepoint fix in §3.5 — without it the test errors with
  "current transaction is aborted".

Add to `classification/tests/` (needs a classification fixture — `ClassificationTestDataMixin` as used in
`analysis/tests/test_classifications_node.py:143`):

- **`test_allele_merge_rehomes_classification`** — classification on allele B with a clinical context on B;
  merge B into A; assert `classification.allele == A`, `classification.clinical_context.allele == A`, and
  its grouping entry is under A's `AlleleGrouping`. Covers the §3.5 receiver / §3.9 body.

**Fixture audit.** Every `VariantAllele` created in tests (`create_mock_allele` and the nine direct
`VariantAllele.objects.create` calls) was checked: the "same variant, two builds" cases
(`test_variant_zygosity_preview_extra.py:65`, `test_variant_sample_information.py:33`,
`test_variant.py:251`) use distinct `genome_build`s, and no test links one variant/build twice. No fixture
changes are needed for the new constraint.

Add to `snpdb/tests/test_liftover.py` next to the existing SAME_CONTIG test (`:50`):

- **`test_same_contig_liftover_merges_when_dest_variant_already_linked`** — variant already linked to
  allele B in the destination build; run `_run_liftover_using_same_contig` for allele A; assert the
  variant has one VA, A and B are merged, and the `AlleleLiftover` is `SKIPPED`.

Skip tests that only restate the constraint (a second `VariantAllele.objects.create` raising
`IntegrityError` is Django's behaviour, not ours).

Run: `python3 manage.py test --keepdb snpdb.tests.test_clingen_allele snpdb.tests.test_allele
snpdb.tests.test_liftover analysis.tests variantopedia.tests.test_tagged_variant_grid`. Note `--keepdb`
keeps the test DB schema — the first run after adding 0224 needs the migration applied to it, which
`--keepdb` does automatically; if the test DB somehow holds duplicates from earlier runs the `RunPython`
will clean them.

---

## 5. Order of work

1. §3.1 + §3.2 + §3.3 + §3.4 (code fixes; each is independently safe under the *old* constraint too).
2. §3.5 (`Allele.merge` FK moves, savepoint, `allele_merged_signal`) + §3.9 command/receiver + tests.
3. §3.6 migration, verify against the dev DB (`vg-test2` has 66 groups of both flavours — a good fixture;
   run the dupe-finder from the issue's second comment, or Appendix A, before and after).
4. §3.7 backfill migration.
5. Tests in §4; full `snpdb`, `analysis`, `classification.tests.utils.test_urls`,
   `variantopedia.tests` runs.
6. Post-deploy check on the affected deployment: the issue's dupe-finder script returns 0 groups;
   the tag on `NC_000007.14:g.140753336_140763336del` shows in both builds.

## 6. Out of scope (noted for later)

- `populate_clingen_alleles_for_variants` line 63 excludes VAs with any `clingen_error`, so a
  `ServerError` is never retried by the bulk path even though `VariantAllele.needs_clingen_call()` says it
  should be. Separate issue.
- `Variant.allele` uses an unordered `.first()`. After 0224 a variant has at most one allele per build,
  and its cross-build rows (shared contigs) share an allele, so this is fine; no change.
- `VariantAlleleSource` / `VariantAlleleCollectionSource` / `VariantAlleleCollectionRecord` have no
  writers anywhere, and `AlleleLiftover.variant_allele` is never assigned. Candidates for removal in a
  separate tidy-up.
- `Allele.merge()` moves `ClinicalContext`s keyed on `name` only, but the model is unique on
  `(allele, allele_origin_bucket, name)` — a germline/somatic pair with the same name would raise. Not
  triggered by this issue; fix when touching the merge (§3.5) if convenient.

---

## Appendix A — repro / verification script (`manage.py shell < script.py`)

Reproduces the leak inside a rolled-back transaction and counts existing duplicate groups. Safe to run on
any DB.

```python
from django.db import transaction
from django.db.models import Count
from snpdb.models import GenomeBuild, Variant, VariantAllele
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI
from snpdb.clingen_allele import populate_clingen_alleles_for_variants

dupes = VariantAllele.objects.values("variant", "genome_build").annotate(n=Count("id")).filter(n__gt=1)
print("existing dupe groups:", dupes.count())

class Rollback(Exception):
    pass

grch38 = GenomeBuild.grch38()
v = Variant.objects.filter(Variant.get_contigs_q(grch38), svlen__lte=-10001, variantallele__isnull=True).first()
print("variant:", v, "skip reason:", v.clingen_allele_skip_reason())
try:
    with transaction.atomic():
        for i in range(3):
            populate_clingen_alleles_for_variants(grch38, [v], clingen_api=MockClinGenAlleleRegistryAPI())
            print(f"after call {i + 1}:", VariantAllele.objects.filter(variant=v, genome_build=grch38).count())
        raise Rollback()
except Rollback:
    print("rolled back")
```

Before the fix this prints 1, 2, 3; after §3.1 it prints 1, 1, 1.
