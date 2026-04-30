# Sample/Cohort Stats Refactor Plan

Closes:
- SACGF/variantgrid#853 — lazy-regenerate stats for current annotation version
- SACGF/variantgrid_private#1038 — cache node counts for CohortNode/TrioNode/PedigreeNode

## Goals

1. Replace the per-Sample stats objects with per-Cohort stats objects. Every Sample already has an auto-created 1-sample Cohort (via VCF), so this is a strict generalisation — no special case for single samples.
2. Fix multi-sample analysis nodes (Cohort/Trio/Pedigree) being slow by giving them the same cached-count fast path SampleNode has today.
3. Lazily regenerate stats when the active annotation version moves on, instead of leaving them stale forever.
4. Drop the `*PassingFilter` model duplication in favour of a `passing_filter` boolean.
5. Invalidate caches when `Cohort.version` is bumped — handled automatically via `on_delete=CASCADE` from the existing `CohortGenotypeCollection` lifecycle.

Out of scope (deferred): caching arbitrary node outputs across analyses, materialised filter-result caches, splitting display vs node-count storage.

## Current shape (for reference)

9 stats models, all keyed by Sample:

| Model | Keyed by | Source |
| --- | --- | --- |
| `SampleStats` / `SampleStatsPassingFilter` | sample | `snpdb/models/models_vcf.py:698-703` |
| `SampleVariantAnnotationStats` (+PF) | sample, variant_annotation_version | `annotation/models/models_sample_stats.py:153-158` |
| `SampleGeneAnnotationStats` (+PF) | sample, gene_annotation_version | `annotation/models/models_sample_stats.py:161-166` |
| `SampleClinVarAnnotationStats` (+PF) | sample, clinvar_version | `annotation/models/models_sample_stats.py:169-174` |

Populated by `annotation/tasks/calculate_sample_stats.py:_actually_calculate_vcf_stats` post-import. Read by `analysis/models/nodes/sources/sample_node.py:_get_cached_label_count`. CohortNode/TrioNode/PedigreeNode have no equivalent — they always run the live count query.

## Target shape

Same number of model classes, but reparented to Cohort and with `passing_filter` collapsed onto the row:

| Model | Keyed by | Notes |
| --- | --- | --- |
| `CohortGenotypeStats` | cohort_genotype_collection, sample (nullable), filter_key (nullable), passing_filter | Replaces `SampleStats(+PF)`. Stats derived from packed CohortGenotype data. |
| `CohortGenotypeVariantAnnotationStats` | cohort_genotype_collection, sample (nullable), variant_annotation_version, filter_key (nullable), passing_filter | Replaces `SampleVariantAnnotationStats(+PF)` |
| `CohortGenotypeGeneAnnotationStats` | cohort_genotype_collection, sample (nullable), gene_annotation_version, filter_key (nullable), passing_filter | Replaces `SampleGeneAnnotationStats(+PF)` |
| `CohortGenotypeClinVarAnnotationStats` | cohort_genotype_collection, sample (nullable), clinvar_version, filter_key (nullable), passing_filter | Replaces `SampleClinVarAnnotationStats(+PF)` |

Row dimensions:

- `sample IS NULL` → aggregate row for the whole cohort (used by CohortNode/TrioNode/PedigreeNode).
- `sample IS NOT NULL` → per-sample row, scoped under `cohort` (typically the sample's VCF cohort; used by SampleNode and the Sample detail page).
- `filter_key IS NULL` → no extra filter applied beyond the cohort/sample identity (the raw count). Always populated.
- `filter_key IS NOT NULL` → a deterministic canonical string identifying a node-side filter config the calc function pre-populated. Trios pre-populate keys for all 5 inheritance modes (denovo / autosomal_recessive / autosomal_dominant / x_linked_recessive / compound_het); future node types can register their own keys without a schema change.

Per-sample rows always have `filter_key IS NULL` — extra filtering is a property of the consuming node, not the sample.

### `filter_key` encoding

Single canonicalisation helper (lives in a new module shared by writers and readers):

```python
import json

def canonical_filter_key(filter_dict: dict | None) -> str | None:
    """ Deterministic string encoding of a filter config. None ↔ None.
        Sorted keys, no whitespace; round-trip stable. """
    if not filter_dict:
        return None
    return json.dumps(filter_dict, sort_keys=True, separators=(",", ":"))
```

Examples:

- TrioNode with denovo: `{"inheritance":"denovo"}`
- TrioNode with autosomal recessive + min_ad=10 (if we ever pre-compute that): `{"inheritance":"autosomal_recessive","min_ad":10}`
- Hypothetical CohortNode with per-sample zygosity: `{"per_sample_zygosity":[["s1","HET"],["s2","HOM"]]}`

The exact filter dicts each node knows how to build are owned by the node class; what the cache cares about is exact string equality on `filter_key`. We never query *into* the key.

### Cacheable-key registry

A small registry in `analysis/models/nodes/sources/_stats_cache.py`:

```python
class FilterKeyHandler:
    """ A node-class-specific helper that knows:
          1. Given an instance of the node, what filter_key (or 'uncacheable') describes its current config.
          2. Given a cohort, what filter_keys the calc function should pre-populate. """
    def filter_key_for_node(self, node) -> str | None | _Uncacheable: ...
    def filter_keys_to_precompute(self, cohort) -> Iterable[str | None]: ...

# Registered per node class:
FILTER_KEY_HANDLERS: dict[type, FilterKeyHandler] = {
    SampleNode: NoFilterHandler(),
    CohortNode: NoFilterHandler(),  # quality filters defeat cache via the existing path
    TrioNode:   TrioInheritanceHandler(),
    PedigreeNode: NoFilterHandler(),
}
```

Adding a new cacheable shape later = add a handler, no schema migration.

Stats FK directly to `CohortGenotypeCollection` (CGC), the per-`(cohort, cohort_version)` pointer that already exists in the codebase. This:

- **Matches the lifecycle precisely.** A VCF cohort is immutable, so its CGC is created once at VCF import and lives forever — per-sample stats hanging off it live as long as the genotype data does. A custom (multi-VCF) cohort bumps `cohort.version`, creates a new CGC, and the old CGC is `marked_for_deletion=True` and async-swept by `delete_old_cohort_genotypes_task`. CASCADE on the FK wipes stats keyed to the old CGC at the same time. No separate `delete_old_counts` clause for stats needed.
- **Drops the standalone `cohort_version` column.** The CGC carries its own `cohort_version`, reachable via the FK. Reads filter on the cohort's *current* CGC, the same way today's CohortNode does.
- **Justifies the `CohortGenotype*Stats` family name.** The stats are explicitly the count summary of one CGC.

We keep the count fields concrete (no JSON) — schema is stable and column queries are easier to reason about.

### Sample vs aggregate rows

Each model has a nullable `sample` FK alongside the (non-null) `cohort` FK:

- `sample IS NOT NULL` → per-sample stats. `cohort` is the parent the sample was imported under (typically `sample.vcf.cohort`). Read by SampleNode and the Sample detail page.
- `sample IS NULL` → aggregate stats for the whole cohort (counts across all members). Read by CohortNode / TrioNode / PedigreeNode.

No fake cohorts, no extra rows in Cohort grids/dropdowns. Sample-row count matches today's `Sample*Stats` row count one-for-one; aggregate rows are new and only created where there's an actual multi-sample consumer.

An N-sample VCF import writes **N+1 stats rows per stat class per filter state** in a single variant-iteration pass:
- 1 aggregate row keyed by `(cohort_genotype_collection=vcf.cohort.cohort_genotype_collection, sample=NULL, filter_key=NULL)`
- N per-sample rows keyed by `(cohort_genotype_collection=vcf.cohort.cohort_genotype_collection, sample=s, filter_key=NULL)` for each s in `vcf.sample_set`

For a custom Cohort built from samples that already have their per-sample rows from prior imports, `calculate_cohort_stats` only writes the new aggregate row(s).

For a Trio cohort, the registered `TrioInheritanceHandler.filter_keys_to_precompute` returns 5 canonical strings — `{"inheritance":"denovo"}`, `{"inheritance":"autosomal_recessive"}`, `{"inheritance":"autosomal_dominant"}`, `{"inheritance":"x_linked_recessive"}`, `{"inheritance":"compound_het"}` — and `calculate_cohort_stats` writes one extra aggregate row per stat class per filter state per key (against the trio cohort's CGC).

The first four are pure per-variant predicates over mom/dad/proband zygosity in the packed CohortGenotype array — each variant is evaluated and added to whichever inheritance buckets it qualifies for (potentially multiple, e.g. an AR variant on chrX is also XLR).

Compound het needs gene-level context: a variant counts as chet if its gene has at least one het variant inherited from each parent in the proband. The calc function handles this in the same iteration with a per-gene accumulator: as it visits each variant it tags `(parent_source, label)` per gene; after the iteration, post-processing sums the chet contribution per label/passing_filter from genes that have hets from both parents. Per-gene memory is bounded (a few small int counters), so even for 30k+ genes it's trivial.

### Unique constraints

Postgres treats NULL as distinct in unique indexes, so each model needs three partial unique constraints to cover the {sample IS NULL/NOT NULL} × {filter_key IS NULL/NOT NULL} matrix (the per-sample × filter-keyed combination is intentionally invalid). Per-sample rows always have `filter_key IS NULL`:

```python
class CohortGenotypeStats(AbstractCohortGenotypeStats):
    cohort_genotype_collection = models.ForeignKey(
        CohortGenotypeCollection, on_delete=CASCADE, related_name="genotype_stats")
    sample = models.ForeignKey(Sample, null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)  # canonical JSON; NULL for raw aggregate / per-sample
    passing_filter = models.BooleanField(default=False)

    class Meta:
        constraints = [
            # per-sample, no filter (the only valid sample-keyed shape)
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_stats_per_sample_uniq",
            ),
            # aggregate, no filter
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_stats_aggregate_uniq",
            ),
            # aggregate, filter-keyed (e.g. trio inheritance)
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_stats_aggregate_filter_uniq",
            ),
        ]
        indexes = [models.Index(fields=["sample", "passing_filter"])]
```

(Same shape on the three annotation-keyed models, with the corresponding `*_annotation_version` field added to all three constraints.)

## Model definitions

New file: `snpdb/models/models_cohort_stats.py` (Cohort lives in snpdb, so the base + genotype-level model goes here to keep the dependency direction).

```python
from django.db import models
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel

from snpdb.models import Cohort, SampleStatsCodeVersion
from snpdb.models.models_enums import ImportStatus
from patients.models_enums import Sex


class AbstractCohortGenotypeStats(TimeStampedModel):
    """ Base for genotype-level (no annotation version) cohort stats.
        Counts derived from packed CohortGenotype data — zygosity, variant class, x_het/x_hom.

        sample IS NOT NULL → per-sample stats; cohort is the parent the sample was imported under.
        sample IS NULL     → aggregate stats across cohort.get_samples(). """
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="genotype_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices,
                                     default=ImportStatus.CREATED)

    variant_count = models.IntegerField(default=0)
    snp_count = models.IntegerField(default=0)
    insertions_count = models.IntegerField(default=0)
    deletions_count = models.IntegerField(default=0)
    ref_count = models.IntegerField(default=0)
    het_count = models.IntegerField(default=0)
    hom_count = models.IntegerField(default=0)
    unk_count = models.IntegerField(default=0)
    x_hom_count = models.IntegerField(default=0)
    x_het_count = models.IntegerField(default=0)
    x_unk_count = models.IntegerField(default=0)

    class Meta:
        abstract = True

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        count = 0
        if zygosity_ref:
            count += self.ref_count
        if zygosity_het:
            count += self.het_count
        if zygosity_hom:
            count += self.hom_count
        if zygosity_unk:
            count += self.unk_count
        return count


class CohortGenotypeStats(AbstractCohortGenotypeStats):
    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_stats_aggregate_filter_uniq",
            ),
        ]

    @property
    def chrx_sex_guess(self):
        sex = Sex.UNKNOWN
        if self.x_het_count and self.x_hom_count:
            ratio = self.x_hom_count / self.x_het_count
            if ratio < 0.2:
                sex = Sex.FEMALE
            elif ratio > 0.8:
                sex = Sex.MALE
        return Sex(sex).label
```

New file: `annotation/models/models_cohort_stats.py` (annotation-version-keyed stats stay in annotation, same dependency direction as today):

```python
from django.db import models
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel

from annotation.models.models import (
    VariantAnnotationVersion, GeneAnnotationVersion, ClinVarVersion,
)
from snpdb.models import Cohort, SampleStatsCodeVersion
from snpdb.models.models_enums import BuiltInFilters


class CohortGenotypeVariantAnnotationStats(TimeStampedModel):
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="variant_annotation_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)

    variant_dbsnp_count = models.IntegerField(default=0)
    insertions_dbsnp_count = models.IntegerField(default=0)
    snp_dbsnp_count = models.IntegerField(default=0)
    deletions_dbsnp_count = models.IntegerField(default=0)
    ref_high_or_moderate_count = models.IntegerField(default=0)
    het_high_or_moderate_count = models.IntegerField(default=0)
    hom_high_or_moderate_count = models.IntegerField(default=0)
    unk_high_or_moderate_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "variant_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_var_ann_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "variant_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_var_ann_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "variant_annotation_version", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_var_ann_stats_aggregate_filter_uniq",
            ),
        ]

    @classmethod
    def load_for(cls, cohort, sample, annotation_version, passing_filter, filter_key=None):
        return cls.objects.get(
            cohort_genotype_collection=cohort.cohort_genotype_collection,
            sample=sample,  # None for aggregate
            variant_annotation_version=annotation_version.variant_annotation_version,
            filter_key=filter_key,
            passing_filter=passing_filter,
        )

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label):
        # Only IMPACT_HIGH_OR_MODERATE is a real stat here today.
        counts = {
            'ref': self.ref_high_or_moderate_count,
            'het': self.het_high_or_moderate_count,
            'hom': self.hom_high_or_moderate_count,
            'unk': self.unk_high_or_moderate_count,
        }
        total = 0
        if zygosity_ref:
            total += counts['ref']
        if zygosity_het:
            total += counts['het']
        if zygosity_hom:
            total += counts['hom']
        if zygosity_unk:
            total += counts['unk']
        return total


class CohortGenotypeGeneAnnotationStats(TimeStampedModel):
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="gene_annotation_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    gene_annotation_version = models.ForeignKey(GeneAnnotationVersion, on_delete=CASCADE)

    gene_count = models.IntegerField(default=0)
    ref_omim_phenotype_count = models.IntegerField(default=0)
    het_omim_phenotype_count = models.IntegerField(default=0)
    hom_omim_phenotype_count = models.IntegerField(default=0)
    unk_omim_phenotype_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "gene_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_gene_ann_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "gene_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_gene_ann_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "gene_annotation_version", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_gene_ann_stats_aggregate_filter_uniq",
            ),
        ]

    @classmethod
    def load_for(cls, cohort, sample, annotation_version, passing_filter, filter_key=None):
        return cls.objects.get(
            cohort_genotype_collection=cohort.cohort_genotype_collection,
            sample=sample,
            gene_annotation_version=annotation_version.gene_annotation_version,
            filter_key=filter_key,
            passing_filter=passing_filter,
        )

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        total = 0
        if zygosity_ref:
            total += self.ref_omim_phenotype_count
        if zygosity_het:
            total += self.het_omim_phenotype_count
        if zygosity_hom:
            total += self.hom_omim_phenotype_count
        if zygosity_unk:
            total += self.unk_omim_phenotype_count
        return total


class CohortGenotypeClinVarAnnotationStats(TimeStampedModel):
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="clinvar_annotation_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    clinvar_version = models.ForeignKey(ClinVarVersion, on_delete=CASCADE)

    clinvar_count = models.IntegerField(default=0)
    ref_clinvar_pathogenic_count = models.IntegerField(default=0)
    het_clinvar_pathogenic_count = models.IntegerField(default=0)
    hom_clinvar_pathogenic_count = models.IntegerField(default=0)
    unk_clinvar_pathogenic_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "clinvar_version", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_clinvar_ann_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "clinvar_version", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_clinvar_ann_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "clinvar_version", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_clinvar_ann_stats_aggregate_filter_uniq",
            ),
        ]

    @classmethod
    def load_for(cls, cohort, sample, annotation_version, passing_filter, filter_key=None):
        return cls.objects.get(
            cohort_genotype_collection=cohort.cohort_genotype_collection,
            sample=sample,
            clinvar_version=annotation_version.clinvar_version,
            filter_key=filter_key,
            passing_filter=passing_filter,
        )

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        total = 0
        if zygosity_ref:
            total += self.ref_clinvar_pathogenic_count
        if zygosity_het:
            total += self.het_clinvar_pathogenic_count
        if zygosity_hom:
            total += self.hom_clinvar_pathogenic_count
        if zygosity_unk:
            total += self.unk_clinvar_pathogenic_count
        return total
```

## Calculation pipeline

Rename the entry point to reflect that it now operates on a cohort, but keep `calculate_vcf_stats` as a thin wrapper for the import path that already calls it.

`annotation/tasks/calculate_sample_stats.py`:

```python
@celery.shared_task
def calculate_vcf_stats(vcf_id, annotation_version_id):
    """ Compute stats for the VCF's auto-cohort, and any cohorts that pre-existed
        but reference samples in this VCF (rare). """
    vcf = VCF.objects.get(pk=vcf_id)
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
    assert annotation_version.genome_build == vcf.genome_build
    try:
        calculate_cohort_stats(vcf.cohort, annotation_version)
        vcf.sample_set.update(import_status=ImportStatus.SUCCESS)
    except Exception:
        create_event(None, 'calculate_vcf_stats', details=get_traceback(), severity=LogLevel.ERROR)
        vcf.sample_set.update(import_status=ImportStatus.ERROR)
        raise


@celery.shared_task
def calculate_cohort_stats_task(cohort_id, annotation_version_id):
    cohort = Cohort.objects.get(pk=cohort_id)
    annotation_version = AnnotationVersion.objects.get(pk=annotation_version_id)
    calculate_cohort_stats(cohort, annotation_version)


def calculate_cohort_stats(cohort: Cohort, annotation_version: AnnotationVersion):
    """ Iterate the cohort's variant queryset once; populate all 4 stats classes
        for every member sample (sample IS NOT NULL rows) AND for the cohort
        as a whole (sample IS NULL aggregate row), for both passing_filter=False
        and (if filters exist) True.

        All rows are written against `cohort.cohort_genotype_collection` (the
        current CGC for this cohort version). When the cohort is later
        version-bumped, a new CGC replaces the old one; ON DELETE CASCADE
        wipes these stats when the old CGC is deleted.

        Filter keys to pre-populate are gathered from registered handlers:
          cgc = cohort.cohort_genotype_collection
          extra_keys = []
          for handler in registered_handlers_for(cohort):
              extra_keys.extend(handler.filter_keys_to_precompute(cohort))

        Each variant is dispatched in a single pass to:
          * the per-sample accumulator (sample=s, filter_key=NULL) for every member
          * the raw aggregate accumulator (sample=NULL, filter_key=NULL)
          * each (sample=NULL, filter_key=key) bucket the variant qualifies for,
            as decided by the handler that registered that key

        For a Trio cohort the registered TrioInheritanceHandler returns five
        canonical-JSON keys:
            {"inheritance":"denovo"}, {"inheritance":"autosomal_recessive"},
            {"inheritance":"autosomal_dominant"}, {"inheritance":"x_linked_recessive"},
            {"inheritance":"compound_het"}.
        Compound het uses a per-gene accumulator built up during the same
        variant iteration, then post-processed into per-label counts after the
        iteration completes.

        For a custom multi-sample Cohort built from samples that already have
        per-sample rows from prior imports, only the aggregate row(s) are written. """
    # Implementation mirrors today's _actually_calculate_vcf_stats but:
    #   * iterates cohort.get_variant_qs() once
    #   * inner accumulator keyed by (sample_id_or_None, filter_key_or_None, passing_filter)
    #   * handler-driven dispatch decides which filter_key buckets each variant lands in
    #   * bulk_create one row per accumulator key per stat class, FK'd to cgc
    ...
```

Sample-level stats are read by `cohort_genotype_collection=sample.vcf.cohort.cohort_genotype_collection, sample=sample`. Aggregate cohort stats are read with `cohort_genotype_collection=cohort.cohort_genotype_collection, sample=None`.

## Reader changes

`analysis/models/nodes/sources/sample_node.py:_get_cached_label_count`:

```python
def _get_cached_label_count(self, label):
    if self._has_filters_that_affect_label_counts():
        return None
    return _get_cached_label_count_for_cohort(
        cohort=self.sample.vcf.cohort,
        sample=self.sample,
        filter_key=None,  # SampleNode never applies a node-side filter
        annotation_version=self.analysis.annotation_version,
        passing_filter=bool(self.get_filter_code()),  # 0 or 1 → False/True
        zygosities=self._effective_zygosities(),
        label=label,
    )
```

New helper, lives in `analysis/models/nodes/sources/_stats_cache.py` (shared by SampleNode/CohortNode/TrioNode/PedigreeNode):

```python
from django.core.exceptions import ObjectDoesNotExist
from snpdb.models import CohortGenotypeStats
from annotation.models import (
    CohortGenotypeVariantAnnotationStats, CohortGenotypeGeneAnnotationStats, CohortGenotypeClinVarAnnotationStats,
)
from snpdb.models.models_enums import BuiltInFilters

LABEL_TO_MODEL = {
    BuiltInFilters.TOTAL: CohortGenotypeStats,
    BuiltInFilters.IMPACT_HIGH_OR_MODERATE: CohortGenotypeVariantAnnotationStats,
    BuiltInFilters.OMIM: CohortGenotypeGeneAnnotationStats,
    BuiltInFilters.CLINVAR: CohortGenotypeClinVarAnnotationStats,
}


_UNCACHEABLE = object()  # sentinel returned by node handlers when the node is in a config we don't precompute


def _get_cached_label_count_for_cohort(cohort, sample, filter_key, annotation_version,
                                       passing_filter, zygosities, label):
    """ sample=None for aggregate cohort stats; sample=Sample(...) for per-sample stats.
        filter_key=None for raw counts; canonical JSON string for filter-keyed counts.
        Returns None on cache miss (caller falls back to live count). """
    if filter_key is _UNCACHEABLE:
        return None  # handler decided this config isn't pre-populated
    model = LABEL_TO_MODEL.get(label)
    if model is None:
        return None
    cgc = cohort.cohort_genotype_collection
    try:
        if model is CohortGenotypeStats:
            obj = CohortGenotypeStats.objects.get(
                cohort_genotype_collection=cgc, sample=sample, filter_key=filter_key,
                passing_filter=passing_filter,
            )
        else:
            obj = model.load_for(cohort, sample, annotation_version, passing_filter,
                                 filter_key=filter_key)
    except ObjectDoesNotExist:
        _enqueue_lazy_recompute(cohort, annotation_version)
        return None  # caller falls back to live count this time
    return obj.count_for_zygosity(*zygosities, label=label)


def _enqueue_lazy_recompute(cohort, annotation_version):
    """ Fire-and-forget. Idempotent — task no-ops if the row now exists. """
    from annotation.tasks.calculate_sample_stats import calculate_cohort_stats_task
    calculate_cohort_stats_task.apply_async(
        args=[cohort.pk, annotation_version.pk],
        queue='annotation_workers',
    )
```

`analysis/models/nodes/sources/cohort_node.py` (and trio/pedigree) gain an override that delegates to the helper:

```python
class CohortNode(AbstractCohortBasedNode, AbstractZygosityCountNode):
    ...
    def _get_cached_label_count(self, label):
        if self._has_filters_that_affect_label_counts():
            return None
        if self.cohort is None:
            return None
        handler = FILTER_KEY_HANDLERS[type(self)]
        filter_key = handler.filter_key_for_node(self)
        return _get_cached_label_count_for_cohort(
            cohort=self.cohort,
            sample=None,  # aggregate row
            filter_key=filter_key,
            annotation_version=self.analysis.annotation_version,
            passing_filter=bool(self.get_filter_code()),
            zygosities=self._effective_zygosities(),
            label=label,
        )
```

The TrioNode handler builds the canonical key from the node's `inheritance` field. All 5 inheritance modes are cached:

```python
class TrioInheritanceHandler(FilterKeyHandler):
    _CACHED_MODES = {
        Inheritance.DENOVO: "denovo",
        Inheritance.RECESSIVE: "autosomal_recessive",
        Inheritance.DOMINANT: "autosomal_dominant",
        Inheritance.X_LINKED_RECESSIVE: "x_linked_recessive",
        Inheritance.COMPOUND_HET: "compound_het",
    }

    def filter_key_for_node(self, node):
        if node.inheritance == Inheritance.NONE:
            return None  # raw aggregate row
        mode = self._CACHED_MODES.get(node.inheritance)
        if mode is None:
            return _UNCACHEABLE  # unknown / new mode not yet registered
        return canonical_filter_key({"inheritance": mode})

    def filter_keys_to_precompute(self, cohort):
        if not hasattr(cohort, "trio"):
            return []
        return [canonical_filter_key({"inheritance": m}) for m in self._CACHED_MODES.values()]
```

TrioNode's `_has_filters_that_affect_label_counts` should *not* return True for the inheritance setting — it's now a cache dimension, not a defeating filter. Quality filters (AD/DP/GQ/PL) on `AbstractCohortBasedNode` still return True and defeat the cache as before.
```

(PedigreeNode reads via `self.pedigree.cohort` with a `NoFilterHandler`, returning `filter_key=None` for raw counts and `_UNCACHEABLE` for any pedigree-specific filtering — those configurations are too varied to enumerate.)

The base class `AnalysisNode.get_node_count` call site already understands a `None` return as "fall back to live count" — no caller changes needed.

## Invalidation on cohort version bump

No code changes needed. `Cohort.increment_version` already marks old `CohortGenotypeCollection` rows for deletion and queues `delete_old_cohort_genotypes_task` to actually delete them (see `snpdb/models/models_cohort.py:121-128`). Stats FK to `CohortGenotypeCollection` with `on_delete=CASCADE`, so when the sweep deletes an old CGC, the stats keyed to it go with it.

Two consequences worth noting:

- **Reads must filter on the *current* CGC.** A consumer with a Cohort grabs `cohort.cohort_genotype_collection` (which always returns the current-version CGC) before looking up stats. Stats for older versions still exist between `marked_for_deletion=True` and the async sweep, but they're FK'd to a different (older) CGC and never seen.
- **VCF cohorts never bump version**, so their CGC and stats live as long as the VCF does. CASCADE is silent in that case — no churn.

If we ever want to wipe stats *eagerly* on version bump (before the async sweep), we can add an explicit clause to `Cohort.increment_version`. Not needed for correctness — reads scope to the current CGC — but a future tightening if old rows pile up.

## Sample detail page (`snpdb/views/views.py:_sample_stats`)

Reads switch to the new models keyed by `(cohort=sample.vcf.cohort, sample=sample)`. Functionally identical: pulls the four stats objects + their passing-filter twins, builds the same DataFrames, renders the same template.

```python
def _sample_stats(sample):
    annotation_version = AnnotationVersion.latest(sample.genome_build)
    cgc = sample.vcf.cohort.cohort_genotype_collection
    base_kwargs = {"cohort_genotype_collection": cgc, "sample": sample}
    cls_kwargs_pf = [
        (CohortGenotypeStats, base_kwargs),
        (CohortGenotypeVariantAnnotationStats, {**base_kwargs,
            "variant_annotation_version": annotation_version.variant_annotation_version}),
        (CohortGenotypeGeneAnnotationStats, {**base_kwargs,
            "gene_annotation_version": annotation_version.gene_annotation_version}),
        (CohortGenotypeClinVarAnnotationStats, {**base_kwargs,
            "clinvar_version": annotation_version.clinvar_version}),
    ]
    ...  # existing dataframe assembly logic, parameterised on passing_filter
```

## Data migration

One Django migration per app, run in order. The biggest deployment is ~30k samples × 4 stat classes × 2 filter states ≈ 240k rows max — backfill in a single transaction is fine. No fake cohorts created.

`snpdb/migrations/NNNN_cohort_genotype_stats.py` (sketch):

```python
def forwards(apps, schema_editor):
    SampleStats = apps.get_model("snpdb", "SampleStats")
    SampleStatsPassingFilter = apps.get_model("snpdb", "SampleStatsPassingFilter")
    CohortGenotypeStats = apps.get_model("snpdb", "CohortGenotypeStats")

    rows = []
    qs = SampleStats.objects.select_related("sample__vcf__cohort__cohort_genotype_collection")
    for ss in qs.iterator():
        cgc = ss.sample.vcf.cohort.cohort_genotype_collection
        rows.append(CohortGenotypeStats(
            cohort_genotype_collection=cgc,
            sample=ss.sample,
            passing_filter=False,
            code_version=ss.code_version,
            import_status=ss.import_status,
            variant_count=ss.variant_count,
            snp_count=ss.snp_count,
            ...
        ))
    pf_qs = SampleStatsPassingFilter.objects.select_related("sample__vcf__cohort__cohort_genotype_collection")
    for ss in pf_qs.iterator():
        rows.append(CohortGenotypeStats(..., passing_filter=True, ...))
    CohortGenotypeStats.objects.bulk_create(rows, batch_size=2000)
```

Note: relies on `Cohort.cohort_genotype_collection` returning the current-version CGC. For VCF cohorts (the only ones with existing per-sample stats to migrate) this is always set, since CGCs are created at VCF import.

Same shape for the three annotation-keyed models in an `annotation/migrations/NNNN_cohort_annotation_stats.py`. No aggregate (sample=NULL) rows are written by the migration — those will be filled lazily by `calculate_cohort_stats` when a CohortNode/TrioNode/PedigreeNode first asks for them.

After the migration is applied and verified in production, a follow-up migration drops the old `Sample*Stats*` tables.

## Touch sites

The full punch list of places that import or use the 9 old models, beyond the already-mentioned `_get_cached_label_count` / `_sample_stats` / `calculate_sample_stats` task / `calculate_sample_stats` mgmt command:

| File:line | What | Refactor action |
| --- | --- | --- |
| `snpdb/models/models_vcf.py:372,377,382-388` | `Sample.delete_internal_data()` deletes all 9 reverse accessors (`samplestats`, `samplestatspassingfilter`, `samplevariantannotationstats_set`, etc.) | Replace with deletion of new models filtered by `sample=self` (per-sample rows only — aggregate rows die when CGC dies) |
| `snpdb/views/views.py:48,62-63,73,252-255` | `view_vcf()` / `_get_vcf_sample_stats()` reads `samplestats` per-sample for the VCF samples table; comment on `prefetch_related_objects("sample_set__samplestats")` | Switch to prefetching `sample_set` joined with new model on `cohort_genotype_collection=vcf.cohort.cohort_genotype_collection, sample__in=...` |
| `snpdb/templates/snpdb/data/view_sample.html:109,163,197` | Reverse OneToOne `sample.samplestats` used for `import_status`, `chrx_sex_guess`, `x_hom_count`, `x_het_count` | Pass the resolved row into the template context from the view; remove the reverse accessor uses |
| `snpdb/management/commands/find_duplicate_samples.py:14-15` | `Sample.objects.filter(samplestats__variant_count__isnull=False).values_list("samplestats__variant_count")` (reverse-relation query) | Rewrite as a join through the new model with `sample=OuterRef('pk'), filter_key__isnull=True, sample__isnull=False` (or similar) |
| `snpdb/models.py:23`, `annotation/models.py:15-17`, `analysis/models/nodes/sources/sample_node.py:12-14,17` | Rollup imports of the old class names | Replace with new class names |

No admin pages, DRF serializers, signal handlers, exports, or test fixtures reference these models directly. Clean surface area beyond what's listed above.

Recent migrations to be aware of as we sequence the schema change:
- `snpdb/migrations/0070`, `0071` — added/altered `code_version` FK on `SampleStats(+PF)`
- `annotation/migrations/0049`, `0050` (one-off backfill), `0051` — added/altered `code_version` FK on the 4 annotation stats models

The schema-add migration for the new models slots in after these. The `0050` backfill pattern (one-off historical fixup) is a good template for the data migration.

## Step-by-step

1. **Add new stats models** (`CohortGenotypeStats`, `CohortGenotypeVariantAnnotationStats`, `CohortGenotypeGeneAnnotationStats`, `CohortGenotypeClinVarAnnotationStats`) with FK to `CohortGenotypeCollection`, nullable `sample`, nullable `filter_key`, partial unique constraints. Schema migration only.
2. **Add `calculate_cohort_stats(cohort, annotation_version)`** in `annotation/tasks/calculate_sample_stats.py`, plus `calculate_cohort_stats_task`. Fork from `_actually_calculate_vcf_stats`; iterate variants once and write per-sample rows + an aggregate row + any handler-supplied filter-keyed rows in one pass, all FK'd to the cohort's CGC.
3. **Data migration** copying existing Sample-level stats rows → per-sample (`sample IS NOT NULL`) rows on the new models, FK'd to `sample.vcf.cohort.cohort_genotype_collection`.
4. **Switch readers**:
   - `analysis/models/nodes/sources/sample_node.py` → call `_get_cached_label_count_for_cohort(cohort=sample.vcf.cohort, sample=sample, ...)` (helper resolves the CGC internally).
   - `snpdb/views/views.py:_sample_stats` → read new models keyed by CGC + sample.
   - `snpdb/views/views.py:view_vcf` / `_get_vcf_sample_stats` → switch the VCF samples table to prefetch from the new model.
   - `snpdb/templates/snpdb/data/view_sample.html` → remove `sample.samplestats` reverse accessor uses; pass resolved stats row into context from `_sample_stats`.
   - `snpdb/management/commands/find_duplicate_samples.py` → rewrite the reverse-relation query.
   - `snpdb/models/models_vcf.py:Sample.delete_internal_data()` → delete from new model filtered by `sample=self` (per-sample rows; aggregate rows die with CGC).
   - Add `_get_cached_label_count` to `cohort_node.py`, `trio_node.py`, `pedigree_node.py`, with their respective `FilterKeyHandler`s.
5. **Switch writer.** Rewire `upload/tasks/vcf/genotype_vcf_tasks.py:CalculateVCFStatsTask` to call `calculate_cohort_stats(vcf.cohort, annotation_version)`. Remove the old per-sample writer.
6. **Wire lazy recompute.** `_get_cached_label_count_for_cohort` enqueues `calculate_cohort_stats_task` on miss. Task is idempotent.
7. **Verify CASCADE invalidation.** No `Cohort.increment_version` change needed — CGC sweep already handles it via `on_delete=CASCADE`. Add a test that confirms stats are wiped when `delete_old_cohort_genotypes_task` runs.
8. **Drop old models.** Once new path is verified in production, remove `SampleStats`, `SampleStatsPassingFilter`, `SampleVariantAnnotationStats(+PF)`, `SampleGeneAnnotationStats(+PF)`, `SampleClinVarAnnotationStats(+PF)` and the `AbstractPassingFilter` mixin. `annotation/management/commands/calculate_sample_stats.py` gets renamed/updated to operate on cohorts.

## Testing

- New `analysis/tests/test_cohort_stats_cache.py`: SampleNode hits per-sample (`sample IS NOT NULL`, `filter_key IS NULL`) rows; CohortNode/PedigreeNode hit aggregate (`sample IS NULL`, `filter_key IS NULL`) rows; TrioNode hits filter-keyed rows for all 5 inheritance modes including compound_het; cache miss falls back to live count and enqueues recompute.
- `canonical_filter_key` round-trip test: encoding the same dict twice produces the same string; encoding with different key orderings produces the same string. Reader and writer must use the same helper.
- CGC CASCADE test: build a custom cohort, populate stats, bump version + run `delete_old_cohort_genotypes_task`, confirm old-CGC stats are wiped and new-CGC stats are unaffected. Confirm VCF cohort stats are untouched (no version bump there).
- Annotation version upgrade test: read with new annotation version misses, fires task, second read hits.
- Trio calc test: a trio cohort import populates 6 aggregate rows per stat class per filter state (NULL filter_key + 5 inheritance keys) plus per-sample rows, all in one variant-iteration pass, all FK'd to the trio cohort's CGC.
- Compound het calc test: the per-gene accumulator correctly counts variants only when the gene has hets inherited from both parents in the proband; matches the live-recompute count.
- Sample detail page renders correctly post-migration (existing tests should pass once readers are switched).
- `URLTestCase` smoke for the sample detail and analysis editor pages.

## Invariants

The plan depends on these holding. If a future change breaks one, the design needs a revisit.

- **Per-sample rows always live on the VCF cohort's CGC.** A sample is in exactly one VCF cohort (its `sample.vcf.cohort`); per-sample stats rows are FK'd to that VCF cohort's CGC and nowhere else. A sample appearing in N custom cohorts doesn't get N per-sample rows — those custom cohorts only get aggregate (`sample=NULL`) rows. Single rule: **sample-keyed rows = VCF cohort CGC. Aggregate rows = whichever cohort's CGC**.
- **All cohorts use the same calc shape.** VCF cohorts and custom multi-VCF cohorts both store packed genotype arrays via `CohortGenotype` / `CohortGenotypeCollection`. The calc function works the same for both — there's no separate "non-VCF cohort" code path. Cost scales with cohort size, not cohort type.
- **Trios are immutable.** Mom/dad/proband assignment doesn't change after creation, so trio inheritance rows never go stale via reassignment. (If we ever did add trio editing, it would bump `cohort.version` and CASCADE would handle it.)
- **Filter-key encoding is centralised.** The cache only hits if the writer's `filter_keys_to_precompute` and the reader's `filter_key_for_node` produce byte-identical strings. Use the single `canonical_filter_key()` helper in both directions and the round-trip test in the testing section. A silent mismatch produces cache misses forever — it falls back gracefully but defeats the optimisation.

## Not in this plan (deferred)

- Cross-analysis caching of node outputs (filter-result counts, variant-id sets). Worth revisiting after this lands — the `(cgc, annotation_version, filter_key, ...)` key shape already supports arbitrary node-config strings, so the same table could be reused for richer caching by registering more handlers.
- Splitting "rich display breakdowns" from "node count cache" into separate tables. Current concrete-column model serves both.
- Materialised views, bitmap-style variant-id caches.
