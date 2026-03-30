# Quad Analysis Plan

## Issue
[SACGF/variantgrid_private#1263](https://github.com/SACGF/variantgrid_private/issues/1263) — Add Quad (2 parents + proband + sibling) analysis with the same clean, structured interface as auto-trio.

Clinical motivation: genomic autopsy cases often include both parents and two children. The trio interface is much clearer than PED-based pedigree setup.

---

## What Is a Quad?

A **Quad** is a 4-member family:
- **Mother** (parent, may be affected)
- **Father** (parent, may be affected)
- **Proband** (index case child, always affected)
- **Sibling** (second child, `sibling_affected` flag — usually unaffected)

The sibling adds filtering power: an unaffected sibling sharing a variant is evidence against pathogenicity.

---

## Design: New `Quad` Model + Shared Code with Trio

Trio and Quad are nearly identical in logic. The implementation plan therefore:

1. **Write Trio unit tests first** — safety net before any refactoring
2. **Refactor `trio_node.py`** — extract shared code into reusable functions/helpers
3. **Add `Quad` model and `QuadNode`** — re-using the extracted code
4. **Add views/forms/templates** — close parallel to Trio

The key shared pieces extracted from Trio:

| Shared component | Current location | What it does |
|-----------------|-----------------|-------------|
| `_build_family_zyg_q(cgc, triples)` | Extracted from `AbstractTrioInheritance._get_zyg_q` | Builds a zygosity Q from any list of `(sample, zyg_set, require_zygosity)` triples |
| `family_dominant_errors(mother_affected, father_affected)` | Extracted from `get_trio_inheritance_errors` | Checks that at least one parent is affected |
| `family_xlinked_recessive_errors(proband_sample, mother_affected)` | Extracted from `get_trio_inheritance_errors` | Checks proband is male and mother is unaffected |

---

## Phase 0: Trio Unit Tests (write FIRST)

**New file: `analysis/tests/test_trio_node.py`**

These tests establish a regression baseline. All must pass before Phase 1 refactoring starts.

```python
# analysis/tests/test_trio_node.py

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, TrioNode
from analysis.models.enums import TrioInheritance
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTrioNodeInheritance(TestCase):
    """
    Regression tests for TrioNode inheritance filtering.
    These must all pass before and after any refactoring.
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        # create_fake_trio sets mother_affected=True, father_affected=False
        cls.trio = create_fake_trio(user, cls.grch37)

        # A second trio with no affected parents (for dominant/xlinked negative tests)
        cls.trio_no_affected = create_fake_trio(user, cls.grch37)
        cls.trio_no_affected.mother_affected = False
        cls.trio_no_affected.father_affected = False
        cls.trio_no_affected.save()

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

    def _make_node(self, trio, inheritance, **kwargs):
        return TrioNode.objects.create(
            analysis=self.analysis, trio=trio,
            inheritance=inheritance, **kwargs
        )

    def _q_str(self, node):
        return str(node._get_node_arg_q_dict())

    # ── Source node modes produce non-empty Q ────────────────────────────────

    def test_recessive_produces_q(self):
        node = self._make_node(self.trio, TrioInheritance.RECESSIVE)
        self.assertNotEqual(self._q_str(node), '{}')

    def test_denovo_produces_q(self):
        node = self._make_node(self.trio, TrioInheritance.DENOVO)
        self.assertNotEqual(self._q_str(node), '{}')

    def test_dominant_produces_q(self):
        # trio has mother_affected=True so dominant is valid
        node = self._make_node(self.trio, TrioInheritance.DOMINANT)
        self.assertNotEqual(self._q_str(node), '{}')

    def test_xlinked_recessive_produces_q(self):
        # trio has mother_affected=True which triggers an error, but
        # the node still produces a Q (errors are separate from execution)
        node = self._make_node(self.trio_no_affected, TrioInheritance.XLINKED_RECESSIVE)
        self.assertNotEqual(self._q_str(node), '{}')

    # ── All simple modes produce distinct Q dicts ─────────────────────────────

    def test_all_modes_produce_distinct_q(self):
        simple_modes = [
            TrioInheritance.RECESSIVE,
            TrioInheritance.DENOVO,
            TrioInheritance.DOMINANT,
            TrioInheritance.XLINKED_RECESSIVE,
        ]
        q_strings = []
        for mode in simple_modes:
            node = self._make_node(self.trio, mode)
            q_strings.append(self._q_str(node))
        self.assertEqual(
            len(q_strings), len(set(q_strings)),
            "Each inheritance mode should produce a distinct Q dict"
        )

    # ── X-linked adds contig restriction ──────────────────────────────────────

    def test_xlinked_adds_contig_q(self):
        node = self._make_node(self.trio_no_affected, TrioInheritance.XLINKED_RECESSIVE)
        arg_q_dict = node._get_node_arg_q_dict()
        self.assertIn(None, arg_q_dict,
                      "X-linked recessive should add a contig=X filter under the None alias")

    def test_xlinked_get_contigs_returns_x(self):
        node = self._make_node(self.trio_no_affected, TrioInheritance.XLINKED_RECESSIVE)
        contigs = node._get_node_contigs()
        self.assertIsNotNone(contigs)
        contig_names = {c.name for c in contigs}
        self.assertIn('X', contig_names)

    # ── require_zygosity flag changes parent Q ────────────────────────────────

    def test_require_zygosity_false_changes_q(self):
        strict = self._make_node(self.trio, TrioInheritance.RECESSIVE, require_zygosity=True)
        loose = self._make_node(self.trio, TrioInheritance.RECESSIVE, require_zygosity=False)
        self.assertNotEqual(
            self._q_str(strict), self._q_str(loose),
            "require_zygosity=False should allow UNKNOWN_ZYGOSITY for parents"
        )

    # ── CompHet node type (filter, not source) ────────────────────────────────

    def test_compound_het_requires_one_parent_input(self):
        node = self._make_node(self.trio, TrioInheritance.COMPOUND_HET)
        self.assertEqual(node.min_inputs, 1)
        self.assertEqual(node.max_inputs, 1)

    def test_simple_modes_are_source_nodes(self):
        for mode in [TrioInheritance.RECESSIVE, TrioInheritance.DENOVO,
                     TrioInheritance.DOMINANT, TrioInheritance.XLINKED_RECESSIVE]:
            node = self._make_node(self.trio, mode)
            self.assertEqual(node.max_inputs, 0, f"{mode} should be a source node (max_inputs=0)")

    # ── Validation: dominant requires affected parent ─────────────────────────

    def test_dominant_with_affected_mother_no_errors(self):
        # trio has mother_affected=True
        errors = TrioNode.get_trio_inheritance_errors(self.trio, TrioInheritance.DOMINANT)
        self.assertEqual(errors, [])

    def test_dominant_no_affected_parent_raises_error(self):
        errors = TrioNode.get_trio_inheritance_errors(
            self.trio_no_affected, TrioInheritance.DOMINANT
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("affected parent", errors[0].lower())

    # ── Validation: X-linked requires unaffected mother ───────────────────────

    def test_xlinked_unaffected_mother_no_errors(self):
        errors = TrioNode.get_trio_inheritance_errors(
            self.trio_no_affected, TrioInheritance.XLINKED_RECESSIVE
        )
        self.assertEqual(errors, [])

    def test_xlinked_affected_mother_raises_error(self):
        # trio has mother_affected=True — should trigger X-linked error
        errors = TrioNode.get_trio_inheritance_errors(
            self.trio, TrioInheritance.XLINKED_RECESSIVE
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("mother", errors[0].lower())

    # ── Clone produces same Q ─────────────────────────────────────────────────

    def test_clone_produces_same_q(self):
        node = self._make_node(self.trio, TrioInheritance.RECESSIVE)
        clone = node.save_clone()
        self.assertEqual(
            self._q_str(node), self._q_str(clone),
            "Cloned TrioNode must produce identical Q dict"
        )
```

Run before refactoring:
```bash
python3 manage.py test --keepdb analysis.tests.test_trio_node
```

---

## Phase 1: Refactor `trio_node.py` — Extract Shared Helpers

**File: `analysis/models/nodes/sources/trio_node.py`** (modify in-place)

The goal is to extract two things so `quad_node.py` can re-use them without duplicating logic.

### 1a. Extract `_build_family_zyg_q`

Currently `AbstractTrioInheritance._get_zyg_q` hardcodes 3 samples. Extract the generic part:

```python
# At module level in trio_node.py (above AbstractTrioInheritance)

def _build_family_zyg_q(
    cohort_genotype_collection,
    sample_zyg_require: list[tuple[Sample, set, bool]]
) -> Q:
    """Build a zygosity Q for any number of family members.

    Args:
        cohort_genotype_collection: The CGC to build the Q for.
        sample_zyg_require: List of (sample, zygosity_set, require_zygosity).
            Entries with an empty zygosity_set are silently skipped
            (used for father in X-linked recessive).
    """
    sample_zygosities_dict = {}
    sample_require_zygosity_dict = {}
    for sample, zyg_set, require in sample_zyg_require:
        if zyg_set:
            sample_zygosities_dict[sample] = zyg_set
            sample_require_zygosity_dict[sample] = require
    return cohort_genotype_collection.get_zygosity_q(
        sample_zygosities_dict, sample_require_zygosity_dict
    )
```

Then update `AbstractTrioInheritance._get_zyg_q` to delegate:

```python
def _get_zyg_q(self, cohort_genotype_collection, trio_zyg_data) -> Q:
    """ trio_zyg_data = tuple of mum_zyg_set, dad_zyg_set, proband_zyg_set """
    trio = self.node.trio
    return _build_family_zyg_q(cohort_genotype_collection, [
        (trio.mother.sample, trio_zyg_data[0], self.node.require_zygosity),
        (trio.father.sample, trio_zyg_data[1], self.node.require_zygosity),
        (trio.proband.sample, trio_zyg_data[2], True),  # proband always required
    ])
```

### 1b. Extract inheritance validation helpers

Currently `get_trio_inheritance_errors` contains the validation logic inline. Extract as module-level helpers:

```python
# At module level in trio_node.py

def _dominant_requires_affected_parent_error(mother_affected: bool, father_affected: bool) -> str | None:
    if not (mother_affected or father_affected):
        return "Dominant inheritance requires an affected parent"
    return None


def _xlinked_recessive_errors(proband_sample, mother_affected: bool) -> list[str]:
    errors = []
    proband_is_female = proband_sample.patient and proband_sample.patient.sex == Sex.FEMALE
    if proband_is_female:
        errors.append("X-linked recessive inheritance doesn't currently work with female proband")
    elif mother_affected:
        errors.append("X-linked recessive inheritance requires an unaffected mother")
    return errors
```

Then `get_trio_inheritance_errors` becomes:

```python
@staticmethod
def get_trio_inheritance_errors(trio: Trio, inheritance) -> list[str]:
    errors = []
    if trio:
        if inheritance == TrioInheritance.DOMINANT:
            if err := _dominant_requires_affected_parent_error(trio.mother_affected, trio.father_affected):
                errors.append(err)
        elif inheritance == TrioInheritance.XLINKED_RECESSIVE:
            errors.extend(_xlinked_recessive_errors(trio.proband.sample, trio.mother_affected))
    return errors
```

Re-run Trio tests to confirm zero regressions:
```bash
python3 manage.py test --keepdb analysis.tests.test_trio_node
```

---

## Phase 2: `Quad` Model + Fixture

### 2a. `Quad` model

**File: `snpdb/models/models_cohort.py`** — add after the `Trio` class (~line 683)

```python
class Quad(GuardianPermissionsAutoInitialSaveMixin, SortByPKMixin, TimeStampedModel):
    """Mother + Father + Proband + Sibling.

    Extends the Trio concept to 4 family members. The sibling (typically
    unaffected) narrows down candidate variants because they share the same
    parental genome without sharing the proband's phenotype.
    """
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    mother = models.ForeignKey(CohortSample, related_name='quad_mother', on_delete=CASCADE)
    mother_affected = models.BooleanField(default=False)
    father = models.ForeignKey(CohortSample, related_name='quad_father', on_delete=CASCADE)
    father_affected = models.BooleanField(default=False)
    proband = models.ForeignKey(CohortSample, related_name='quad_proband', on_delete=CASCADE)
    sibling = models.ForeignKey(CohortSample, related_name='quad_sibling', on_delete=CASCADE)
    sibling_affected = models.BooleanField(default=False)

    @classmethod
    def get_permission_class(cls):
        return Cohort

    def get_permission_object(self):
        return self.cohort

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(cohort__in=queryset)

    @property
    def genome_build(self):
        return self.cohort.genome_build

    def get_cohort_samples(self):
        return [self.mother, self.father, self.proband, self.sibling]

    def get_samples(self):
        return Sample.objects.filter(cohortsample__in=self.get_cohort_samples()).order_by("pk")

    def get_absolute_url(self):
        return reverse('view_quad', kwargs={"pk": self.pk})

    def get_listing_url(self):
        return reverse('quads')

    @property
    def mother_details(self):
        affected = "affected" if self.mother_affected else "unaffected"
        return f"{self.mother} ({affected})"

    @property
    def father_details(self):
        affected = "affected" if self.father_affected else "unaffected"
        return f"{self.father} ({affected})"

    @property
    def sibling_details(self):
        affected = "affected" if self.sibling_affected else "unaffected"
        return f"{self.sibling} ({affected})"

    def __str__(self):
        return self.name or f"Quad {self.pk}"
```

**Migration:** `snpdb/migrations/XXXX_add_quad.py`

### 2b. `create_fake_quad()` fixture

**File: `snpdb/tests/utils/fake_cohort_data.py`** — add after `create_fake_trio`

```python
def create_fake_quad(user: User, genome_build: GenomeBuild,
                     sibling_affected: bool = False) -> 'Quad':
    """4-sample Cohort (proband, mother, father, sibling) + a Quad."""
    from snpdb.models import Quad

    vcf = VCF.objects.create(
        name="test_quad_vcf", genotype_samples=1, genome_build=genome_build,
        import_status=ImportStatus.SUCCESS, user=user, date=timezone.now()
    )
    proband_sample = Sample.objects.create(name="proband", vcf=vcf,
                                           import_status=ImportStatus.SUCCESS)
    mother_sample = Sample.objects.create(name="mother", vcf=vcf)
    father_sample = Sample.objects.create(name="father", vcf=vcf)
    sibling_sample = Sample.objects.create(name="sibling", vcf=vcf)

    assign_permission_to_user_and_groups(user, vcf)
    assign_permission_to_user_and_groups(user, proband_sample)

    cohort = Cohort.objects.create(
        name="test_quad_cohort", user=user, vcf=vcf,
        genome_build=genome_build, import_status=ImportStatus.SUCCESS
    )
    for i, sample in enumerate([proband_sample, mother_sample, father_sample, sibling_sample]):
        CohortSample.objects.create(
            cohort=cohort, sample=sample,
            cohort_genotype_packed_field_index=i, sort_order=i
        )
    assign_permission_to_user_and_groups(user, cohort)

    CohortGenotypeCollection.objects.create(
        cohort=cohort, cohort_version=cohort.version,
        num_samples=cohort.cohortsample_set.count()
    )

    proband_cs = cohort.cohortsample_set.get(sample__name='proband')
    mother_cs = cohort.cohortsample_set.get(sample__name='mother')
    father_cs = cohort.cohortsample_set.get(sample__name='father')
    sibling_cs = cohort.cohortsample_set.get(sample__name='sibling')

    return Quad.objects.create(
        name="test_quad",
        user=user,
        cohort=cohort,
        mother=mother_cs, mother_affected=False,
        father=father_cs, father_affected=False,
        proband=proband_cs,
        sibling=sibling_cs, sibling_affected=sibling_affected,
    )
```

---

## Phase 3: `QuadNode` Analysis Node

### 3a. `QuadInheritance` enum

**File: `analysis/models/enums.py`** — add after `TrioInheritance`

```python
class QuadInheritance(models.TextChoices):
    RECESSIVE         = 'R', 'Recessive'
    COMPOUND_HET      = 'C', 'C. Het'
    DOMINANT          = 'D', 'Dominant'
    DENOVO            = 'N', 'Denovo'
    XLINKED_RECESSIVE = 'X', 'X-Linked Recessive'
```

### 3b. `quad_node.py`

**New file: `analysis/models/nodes/sources/quad_node.py`**

Imports the shared helpers from `trio_node` to avoid duplicating logic:

```python
from abc import ABC, abstractmethod
from typing import Optional

from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from django.db import models
from django.db.models import Count
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import QuadInheritance, NodeErrorSource, AnalysisTemplateType
from analysis.models.nodes.sources import AbstractCohortBasedNode
# Re-use the shared helpers extracted from trio_node in Phase 1
from analysis.models.nodes.sources.trio_node import (
    _build_family_zyg_q,
    _dominant_requires_affected_parent_error,
    _xlinked_recessive_errors,
)
from annotation.models.models import VariantTranscriptAnnotation
from library.constants import DAY_SECS
from patients.models_enums import Zygosity, Sex
from snpdb.models import Quad, Sample, Contig


class AbstractQuadInheritance(ABC):
    NO_VARIANT = {Zygosity.MISSING, Zygosity.HOM_REF}
    HAS_VARIANT = {Zygosity.HET, Zygosity.HOM_ALT}
    UNAFFECTED_AND_AFFECTED_ZYGOSITIES = [NO_VARIANT, HAS_VARIANT]

    def __init__(self, node: 'QuadNode'):
        self.node = node

    def _get_zyg_q(self, cgc, quad_zyg_data) -> Q:
        """quad_zyg_data = (mum_zyg, dad_zyg, proband_zyg, sibling_zyg)"""
        quad = self.node.quad
        return _build_family_zyg_q(cgc, [
            (quad.mother.sample,  quad_zyg_data[0], self.node.require_zygosity),
            (quad.father.sample,  quad_zyg_data[1], self.node.require_zygosity),
            (quad.proband.sample, quad_zyg_data[2], True),   # proband always required
            (quad.sibling.sample, quad_zyg_data[3], self.node.require_zygosity),
        ])

    @abstractmethod
    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        pass

    @abstractmethod
    def get_method(self) -> str:
        pass

    def get_contigs(self) -> Optional[set[Contig]]:
        return None


class SimpleQuadInheritance(AbstractQuadInheritance):
    @abstractmethod
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        pass

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.quad.cohort.cohort_genotype_collection
        alias = cgc.cohortgenotype_alias
        q = self._get_zyg_q(cgc, self._get_mum_dad_proband_sibling_zygosities())
        return {alias: {str(q): q}}

    def get_method(self) -> str:
        mum_z, dad_z, prob_z, sib_z = self._get_mum_dad_proband_sibling_zygosities()
        return f"Mother:{mum_z} Father:{dad_z} Proband:{prob_z} Sibling:{sib_z}"


class QuadRecessive(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        # Unaffected sibling must NOT be HOM_ALT — they can be a carrier (HET) or ref
        if self.node.quad.sibling_affected:
            sibling_zyg = self.HAS_VARIANT
        else:
            sibling_zyg = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}
        return {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT}, sibling_zyg


class QuadDominant(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        quad = self.node.quad
        mother_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.mother_affected)]
        father_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.father_affected)]
        sibling_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.sibling_affected)]
        return mother_zyg, father_zyg, self.HAS_VARIANT, sibling_zyg


class QuadDenovo(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        # Both parents AND sibling are NO_VARIANT — strongest de novo confirmation
        return self.NO_VARIANT, self.NO_VARIANT, self.HAS_VARIANT, self.NO_VARIANT


class QuadXLinkedRecessive(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        # Father not checked (set()); sibling filtered by affected status
        if self.node.quad.sibling_affected:
            sibling_zyg = {Zygosity.HOM_ALT}
        else:
            sibling_zyg = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}
        return {Zygosity.HET}, set(), {Zygosity.HOM_ALT}, sibling_zyg

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = super().get_arg_q_dict()
        q = Q(locus__contig__name='X')
        arg_q_dict[None] = {str(q): q}
        return arg_q_dict

    def get_contigs(self) -> Optional[set[Contig]]:
        return set(self.node.quad.genome_build.contigs.filter(name='X'))


class QuadCompHet(AbstractQuadInheritance):
    """Compound Het for Quad. Same two-pass gene logic as TrioCompHet.

    TODO: An unaffected sibling having BOTH comp-het hits (one from mum AND one from dad)
    in the same gene is strong evidence against pathogenicity and should be excluded.
    This would require a third pass in the two-hit gene query. Left for a future iteration.
    (Tracked: SACGF/variantgrid_private#1263)
    """

    def _mum_but_not_dad(self):
        return {Zygosity.HET}, self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}

    def _dad_but_not_mum(self):
        return self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HET}

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.node.pk, s.node.version))
    def _get_comp_het_q_and_two_hit_genes(self):
        cgc = self.node.quad.cohort.cohort_genotype_collection
        parent = self.node.get_single_parent()
        mum_but_not_dad = self._get_zyg_q(cgc, self._mum_but_not_dad())
        dad_but_not_mum = self._get_zyg_q(cgc, self._dad_but_not_mum())
        comp_het_q = mum_but_not_dad | dad_but_not_mum

        annotation_kwargs = self.node.get_annotation_kwargs()

        def get_parent_genes(q):
            qs = parent.get_queryset(q, extra_annotation_kwargs=annotation_kwargs)
            return qs.values_list("varianttranscriptannotation__gene", flat=True).distinct()

        common_genes = set(get_parent_genes(mum_but_not_dad)) & set(get_parent_genes(dad_but_not_mum))
        vav = self.node.analysis.annotation_version.variant_annotation_version
        q_in_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, common_genes)
        parent_genes_qs = parent.get_queryset(q_in_genes, extra_annotation_kwargs=annotation_kwargs)
        parent_genes_qs = parent_genes_qs.values_list("varianttranscriptannotation__gene")
        two_hits = parent_genes_qs.annotate(gene_count=Count("pk")).filter(gene_count__gte=2)
        two_hit_genes = set(two_hits.values_list("varianttranscriptannotation__gene", flat=True).distinct())
        return comp_het_q, two_hit_genes

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        comp_het_q, two_hit_genes = self._get_comp_het_q_and_two_hit_genes()
        vav = self.node.analysis.annotation_version.variant_annotation_version
        comp_het_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, two_hit_genes)
        cgc = self.node.quad.cohort.cohort_genotype_collection
        q_hash = str(comp_het_q)
        return {
            cgc.cohortgenotype_alias: {q_hash: comp_het_q},
            None: {q_hash: comp_het_genes},
        }

    def get_method(self) -> str:
        return "Proband: HET, >=2 hits in gene from (mum OR dad), sibling HET for each hit"

    def get_contigs(self) -> Optional[set[Contig]]:
        _, two_hit_genes = self._get_comp_het_q_and_two_hit_genes()
        contig_qs = Contig.objects.filter(
            transcriptversion__genome_build=self.node.quad.genome_build,
            transcriptversion__gene_version__gene__in=two_hit_genes
        )
        return set(contig_qs.distinct())


class QuadNode(AbstractCohortBasedNode):
    quad = models.ForeignKey(Quad, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=QuadInheritance.choices,
                                   default=QuadInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        return 1 if self.inheritance == QuadInheritance.COMPOUND_HET else 0

    @staticmethod
    def get_quad_inheritance_errors(quad: Quad, inheritance) -> list[str]:
        """Reuses the same validation helpers as TrioNode."""
        errors = []
        if quad:
            if inheritance == QuadInheritance.DOMINANT:
                if err := _dominant_requires_affected_parent_error(
                    quad.mother_affected, quad.father_affected
                ):
                    errors.append(err)
            elif inheritance == QuadInheritance.XLINKED_RECESSIVE:
                errors.extend(
                    _xlinked_recessive_errors(quad.proband.sample, quad.mother_affected)
                )
        return errors

    def get_errors(self, include_parent_errors=True, flat=False):
        errors = super().get_errors(include_parent_errors=include_parent_errors)
        if self.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            if quad_errors := self.get_quad_inheritance_errors(self.quad, self.inheritance):
                errors.extend(((NodeErrorSource.CONFIGURATION, e) for e in quad_errors))
        if flat:
            errors = self.flatten_errors(errors)
        return errors

    def _get_cohort(self):
        return self.quad.cohort if self.quad else None

    def modifies_parents(self):
        return self.quad is not None

    def _inheritance_factory(self):
        classes = {
            QuadInheritance.COMPOUND_HET:    QuadCompHet,
            QuadInheritance.RECESSIVE:       QuadRecessive,
            QuadInheritance.DOMINANT:        QuadDominant,
            QuadInheritance.DENOVO:          QuadDenovo,
            QuadInheritance.XLINKED_RECESSIVE: QuadXLinkedRecessive,
        }
        return classes[QuadInheritance(self.inheritance)](self)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cohort, arg_q_dict = self.get_cohort_and_arg_q_dict()
        if cohort:
            inheritance = self._inheritance_factory()
            self.merge_arg_q_dicts(arg_q_dict, inheritance.get_arg_q_dict())
            self.merge_arg_q_dicts(arg_q_dict, self.get_vcf_locus_filters_arg_q_dict())
        return arg_q_dict

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        if self.quad:
            return self._inheritance_factory().get_contigs()
        return None

    def _get_method_summary(self):
        if self._get_cohort():
            return self._inheritance_factory().get_method()
        return "No cohort selected"

    def get_node_name(self):
        name_parts = [QuadInheritance(self.inheritance).label]
        if not self.require_zygosity:
            name_parts.append('(non strict)')
        if desc := self.get_filter_description():
            name_parts.append(f"({desc})")
        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return "Mother/Father/Proband/Sibling - filter for recessive/dominant/denovo inheritance"

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not self.quad:
            errors.append("No quad selected")
        else:
            errors.extend(self._get_genome_build_errors("quad", self.quad.genome_build))
        return errors

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts, visibility = [], {}
        if self.quad:
            cohort = self.quad.cohort
            cohorts = [cohort]
            visibility = {s: cohort.has_genotype for s in self.quad.get_samples()}
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        return self.quad.proband.sample if self.quad else None

    @staticmethod
    def get_node_class_label():
        return 'Quad'

    def __str__(self):
        return f"QuadNode: {self.pk}"


auditlog.register(QuadNode)
```

**Register in `analysis/models/nodes/sources/__init__.py`:**
```python
from .quad_node import *
```

**Migration:** `analysis/migrations/XXXX_add_quad_node.py`

---

## Phase 4: QuadNode Unit Tests

**New file: `analysis/tests/test_quad_node.py`**

```python
# analysis/tests/test_quad_node.py

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, QuadNode
from analysis.models.enums import QuadInheritance
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_quad


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestQuadNodeInheritance(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.quad = create_fake_quad(user, cls.grch37, sibling_affected=False)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

    def _make_node(self, inheritance, quad=None, **kwargs):
        return QuadNode.objects.create(
            analysis=self.analysis,
            quad=quad or self.quad,
            inheritance=inheritance,
            **kwargs
        )

    def _q_str(self, node):
        return str(node._get_node_arg_q_dict())

    # ── Basic: each mode produces a Q ─────────────────────────────────────────

    def test_recessive_produces_q(self):
        node = self._make_node(QuadInheritance.RECESSIVE)
        self.assertNotEqual(self._q_str(node), '{}')

    def test_denovo_produces_q(self):
        node = self._make_node(QuadInheritance.DENOVO)
        self.assertNotEqual(self._q_str(node), '{}')

    def test_xlinked_recessive_produces_q(self):
        node = self._make_node(QuadInheritance.XLINKED_RECESSIVE)
        self.assertNotEqual(self._q_str(node), '{}')

    def test_all_modes_produce_distinct_q(self):
        simple_modes = [
            QuadInheritance.RECESSIVE,
            QuadInheritance.DENOVO,
            QuadInheritance.XLINKED_RECESSIVE,
        ]
        q_strings = [self._q_str(self._make_node(m)) for m in simple_modes]
        self.assertEqual(len(q_strings), len(set(q_strings)),
                         "Each mode should produce a distinct Q dict")

    # ── Sibling constraint changes the Q ─────────────────────────────────────

    def test_affected_sibling_changes_recessive_q(self):
        """Affected sibling should have HAS_VARIANT, unaffected should not."""
        user = self.analysis.user
        affected_sib_quad = create_fake_quad(user, self.grch37, sibling_affected=True)
        unaffected = self._make_node(QuadInheritance.RECESSIVE)
        affected = self._make_node(QuadInheritance.RECESSIVE, quad=affected_sib_quad)
        self.assertNotEqual(
            self._q_str(unaffected), self._q_str(affected),
            "Affected vs unaffected sibling should change recessive Q"
        )

    def test_affected_sibling_changes_denovo_q(self):
        user = self.analysis.user
        affected_sib_quad = create_fake_quad(user, self.grch37, sibling_affected=True)
        unaffected = self._make_node(QuadInheritance.DENOVO)
        affected = self._make_node(QuadInheritance.DENOVO, quad=affected_sib_quad)
        # De Novo: sibling NO_VARIANT when unaffected, so an "affected sibling"
        # quad is a different (unusual) case — just check the Q differs
        self.assertNotEqual(self._q_str(unaffected), self._q_str(affected))

    # ── X-linked adds contig restriction ──────────────────────────────────────

    def test_xlinked_adds_contig_q(self):
        node = self._make_node(QuadInheritance.XLINKED_RECESSIVE)
        arg_q_dict = node._get_node_arg_q_dict()
        self.assertIn(None, arg_q_dict)

    def test_xlinked_get_contigs_returns_x(self):
        node = self._make_node(QuadInheritance.XLINKED_RECESSIVE)
        contigs = node._get_node_contigs()
        self.assertIsNotNone(contigs)
        self.assertIn('X', {c.name for c in contigs})

    # ── require_zygosity ──────────────────────────────────────────────────────

    def test_require_zygosity_false_changes_q(self):
        strict = self._make_node(QuadInheritance.RECESSIVE, require_zygosity=True)
        loose = self._make_node(QuadInheritance.RECESSIVE, require_zygosity=False)
        self.assertNotEqual(self._q_str(strict), self._q_str(loose))

    # ── Validation ────────────────────────────────────────────────────────────

    def test_dominant_no_affected_parent_raises_error(self):
        errors = QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.DOMINANT)
        self.assertGreater(len(errors), 0)

    def test_dominant_with_affected_mother_no_errors(self):
        self.quad.mother_affected = True
        self.quad.save()
        errors = QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.DOMINANT)
        self.assertEqual(errors, [])
        self.quad.mother_affected = False
        self.quad.save()

    def test_xlinked_affected_mother_raises_error(self):
        self.quad.mother_affected = True
        self.quad.save()
        errors = QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.XLINKED_RECESSIVE)
        self.assertGreater(len(errors), 0)
        self.quad.mother_affected = False
        self.quad.save()

    # ── CompHet node type ─────────────────────────────────────────────────────

    def test_compound_het_requires_one_parent_input(self):
        node = self._make_node(QuadInheritance.COMPOUND_HET)
        self.assertEqual(node.min_inputs, 1)
        self.assertEqual(node.max_inputs, 1)

    def test_simple_modes_are_source_nodes(self):
        for mode in [QuadInheritance.RECESSIVE, QuadInheritance.DENOVO,
                     QuadInheritance.DOMINANT, QuadInheritance.XLINKED_RECESSIVE]:
            node = self._make_node(mode)
            self.assertEqual(node.max_inputs, 0, f"{mode} should be a source node")

    # ── Clone ─────────────────────────────────────────────────────────────────

    def test_clone_produces_same_q(self):
        node = self._make_node(QuadInheritance.RECESSIVE)
        clone = node.save_clone()
        self.assertEqual(self._q_str(node), self._q_str(clone))

    # ── Trio vs Quad: same inheritance mode produces different Q ──────────────

    def test_quad_recessive_q_differs_from_trio_recessive_q(self):
        """The sibling in a Quad should make the Recessive Q tighter."""
        from analysis.models import TrioNode
        from analysis.models.enums import TrioInheritance
        from snpdb.tests.utils.fake_cohort_data import create_fake_trio

        trio = create_fake_trio(self.analysis.user, self.grch37)
        trio_node = TrioNode.objects.create(
            analysis=self.analysis, trio=trio,
            inheritance=TrioInheritance.RECESSIVE
        )
        quad_node = self._make_node(QuadInheritance.RECESSIVE)
        self.assertNotEqual(
            str(trio_node._get_node_arg_q_dict()),
            str(quad_node._get_node_arg_q_dict()),
            "Quad recessive Q must differ from Trio recessive Q (extra sibling constraint)"
        )
```

Run both test suites together:
```bash
python3 manage.py test --keepdb analysis.tests.test_trio_node analysis.tests.test_quad_node
```

---

## Phase 5: Forms, Views, Templates

### 5a. `QuadNodeForm`

**File: `analysis/forms/forms_nodes.py`** — add after `TrioNodeForm`

```python
class QuadNodeForm(GenomeBuildAutocompleteForwardMixin, VCFSourceNodeForm):
    genome_build_fields = ["quad"]

    class Meta:
        model = QuadNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "quad": ModelSelect2(url='quad_autocomplete',
                                 attrs={'data-placeholder': 'Quad...'}),
            "min_ad": WIDGET_INTEGER_MIN_0,
            "min_dp": WIDGET_INTEGER_MIN_0,
            "min_gq": WIDGET_INTEGER_MIN_0,
            "max_pl": WIDGET_INTEGER_MIN_0,
        }

    def clean(self):
        cleaned_data = super().clean()
        quad = cleaned_data.get("quad")
        inheritance = cleaned_data.get("inheritance")
        if self.instance.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            if quad and inheritance:
                for error in QuadNode.get_quad_inheritance_errors(quad, inheritance):
                    self.add_error("inheritance", error)
```

### 5b. `QuadNodeView`

**File: `analysis/views/nodes/node_views.py`** — add after `TrioNodeView`

```python
class QuadNodeView(NodeView):
    model = QuadNode
    form_class = QuadNodeForm

    def get_form_kwargs(self):
        form_kwargs = super().get_form_kwargs()
        form_kwargs["genome_build"] = self.object.analysis.genome_build
        return form_kwargs
```

### 5c. Node editor template

**New file: `analysis/templates/analysis/node_editors/quadnode_editor.html`**

Clone of `trionode_editor.html` with:
- `trio` → `quad` throughout
- `api_view_trio` → `api_view_quad`
- Extra sibling row in the info table
- `sibling_affected` added to the `AFFECTED_FIELDS` JS array

### 5d. Snpdb views, templates, URLs

**`snpdb/views/views.py`** — add after `view_trio`:

```python
def quads(request):
    show_group_data = UserGridConfig.get(request.user, 'Quads').show_group_data
    context = {"show_group_data": show_group_data}
    return render(request, 'snpdb/patients/quads.html', context)


def view_quad(request, pk):
    quad = Quad.get_for_user(request.user, pk)
    context = {"quad": quad,
               "has_write_permission": quad.cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_quad.html', context)
```

**`snpdb/urls.py`** — add alongside trio URLs:
```python
path('quads', views.quads, name='quads'),
path('view_quad/<int:pk>', views.view_quad, name='view_quad'),
path('api/quad/<pk>', views_rest.QuadView.as_view(), name='api_view_quad'),
path('autocomplete/Quad/', views_autocomplete.QuadAutocompleteView.as_view(), name='quad_autocomplete'),
path('quad/datatable/', DatabaseTableView.as_view(column_class=QuadsListColumns), name='quad_datatable'),
```

**Templates:**

- `snpdb/templates/snpdb/patients/quads.html` — mirrors `trios.html`
- `snpdb/templates/snpdb/patients/view_quad.html` — mirrors `view_trio.html` (uses `quad_table` tag)
- `snpdb/templates/snpdb/tags/quad_table.html` — mirrors `trio_table.html` with 4-row table (adds Sibling row)

**`snpdb/templatetags/model_tags.py`** — add after `trio_table`:

```python
@register.inclusion_tag("snpdb/tags/quad_table.html", takes_context=False)
def quad_table(quad: Quad):
    return {"quad": quad}
```

### 5e. Serializer + REST view

**`snpdb/serializers.py`**:

```python
class QuadSerializer(serializers.ModelSerializer):
    mother  = serializers.SerializerMethodField()
    father  = serializers.SerializerMethodField()
    proband = serializers.SerializerMethodField()
    sibling = serializers.SerializerMethodField()

    class Meta:
        model = Quad
        fields = '__all__'

    def get_mother(self, obj):  return obj.mother.name
    def get_father(self, obj):  return obj.father.name
    def get_proband(self, obj): return obj.proband.name
    def get_sibling(self, obj): return obj.sibling.name
```

**`snpdb/views/views_rest.py`**:

```python
@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class QuadView(RetrieveAPIView):
    serializer_class = QuadSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    lookup_fields = ('pk',)

    def get_queryset(self):
        return Quad.filter_for_user(self.request.user)
```

### 5f. Grids/Autocomplete

**`snpdb/grids.py`** — `QuadsListColumns` (mirrors `TriosListColumns`):

```python
class QuadsListColumns(DatatableConfig[Quad]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='user__username', label='User', orderable=True),
            RichColumn(key='modified', client_renderer='TableFormat.timestamp',
                       orderable=True, default_sort=SortOrder.DESC),
        ]

    def get_initial_queryset(self) -> QuerySet[Quad]:
        return Quad.filter_for_user(self.user)

    def filter_queryset(self, qs: QuerySet[Quad]) -> QuerySet[Quad]:
        if not UserGridConfig.get(self.user, 'Quads').show_group_data:
            qs = qs.filter(user=self.user)
        return qs
```

**`snpdb/views/views_autocomplete.py`** — `Quad` has no `import_status` of its own; filter via the cohort's status (same approach as `TrioAutocompleteView`):

```python
@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class QuadAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = Quad.filter_for_user(user).filter(cohort__import_status=ImportStatus.SUCCESS)
        return self.filter_to_genome_build(qs, "cohort__genome_build")
```

### 5g. `related_analyses` for Quad

**`analysis/related_analyses.py`** — add after `get_related_analysis_details_for_trio`:

```python
def get_related_analysis_details_for_quad(user, quads) -> list[tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = Analysis.filter_for_user(user).values_list("pk", flat=True)
    for quad_node in QuadNode.objects.filter(analysis__in=analyses_ids,
                                             quad__in=quads).select_related("analysis", "quad"):
        all_analysis_details[quad_node.analysis].add(str(quad_node.quad))
    return sort_analyses_by_date_and_merge_details(all_analysis_details)
```

**`analysis/templatetags/related_analyses_tags.py`** — add `related_analyses_for_quad` tag mirroring `related_analyses_for_trio`:

```python
@register.inclusion_tag("analysis/tags/related_analyses_for_quad.html", takes_context=True)
def related_analyses_for_quad(context, quad):
    pedigrees = quad.cohort.pedigree_set.all()
    update_context_with_related_analysis(context, quad.get_samples(), [quad.cohort], [quad], pedigrees)
    context["quad"] = quad
    return context
```

Also add `"quad"` to `single_model_args` in `analysis_templates_tag`, and add a proband-as-sample hack parallel to the existing trio hack:

```python
single_model_args = {"sample", "cohort", "trio", "quad", "pedigree"}
...
# existing trio hack:
if trio := kwargs.get("trio"):
    hidden_inputs["sample"] = trio.proband.sample_id
# new quad hack:
if quad := kwargs.get("quad"):
    hidden_inputs["sample"] = quad.proband.sample_id
```

Create `analysis/templates/analysis/tags/related_analyses_for_quad.html` mirroring `related_analyses_for_trio.html`.

---

## File Change Summary

| File | Change |
|------|--------|
| `analysis/tests/test_trio_node.py` | **NEW** — regression tests for Trio |
| `analysis/models/nodes/sources/trio_node.py` | Refactor: extract `_build_family_zyg_q`, `_dominant_requires_affected_parent_error`, `_xlinked_recessive_errors` |
| `analysis/models/enums.py` | Add `QuadInheritance` enum |
| `analysis/models/nodes/sources/quad_node.py` | **NEW** — QuadNode + all inheritance classes |
| `analysis/models/nodes/sources/__init__.py` | Add `from .quad_node import *` |
| `analysis/migrations/XXXX_add_quad_node.py` | New migration |
| `analysis/forms/forms_nodes.py` | Add `QuadNodeForm` |
| `analysis/views/nodes/node_views.py` | Add `QuadNodeView` |
| `analysis/urls.py` | Register `QuadNodeView` URL |
| `analysis/templates/analysis/node_editors/quadnode_editor.html` | **NEW** — node editor template |
| `analysis/related_analyses.py` | Add `get_related_analysis_details_for_quad` |
| `analysis/templatetags/related_analyses_tags.py` | Add `related_analyses_for_quad` tag; add `"quad"` to `single_model_args`; add proband-as-sample hack |
| `analysis/templates/analysis/tags/related_analyses_for_quad.html` | **NEW** — mirrors `related_analyses_for_trio.html` |
| `analysis/tests/test_quad_node.py` | **NEW** — QuadNode tests |
| `snpdb/models/models_cohort.py` | Add `Quad` model |
| `snpdb/migrations/XXXX_add_quad.py` | New migration |
| `snpdb/serializers.py` | Add `QuadSerializer` |
| `snpdb/views/views.py` | Add `quads()`, `view_quad()` |
| `snpdb/views/views_rest.py` | Add `QuadView` |
| `snpdb/views/views_autocomplete.py` | Add `QuadAutocompleteView` |
| `snpdb/urls.py` | Add quad URL patterns |
| `snpdb/grids.py` | Add `QuadsListColumns` |
| `snpdb/templatetags/model_tags.py` | Add `quad_table` tag |
| `snpdb/templates/snpdb/patients/quads.html` | **NEW** — list template |
| `snpdb/templates/snpdb/patients/view_quad.html` | **NEW** — detail template |
| `snpdb/templates/snpdb/tags/quad_table.html` | **NEW** — detail table tag |
| `snpdb/tests/utils/fake_cohort_data.py` | Add `create_fake_quad()` |

---

## Suggested Implementation Order

```
1. analysis/tests/test_trio_node.py         ← write, run, all green
2. trio_node.py refactor                    ← extract helpers, run test_trio_node again
3. Quad model + snpdb migration             ← snpdb model tests
4. create_fake_quad()                       ← fixture
5. enums.py QuadInheritance
6. quad_node.py + analysis migration
7. analysis/tests/test_quad_node.py         ← run both test suites
8. Forms, views, templates, URLs (snpdb + analysis)
9. Serializer, REST view, autocomplete, grids
10. related_analyses + template tags
```

**Out of scope for this plan (separate tasks):**
- Auto-analysis JSON template for Quad
- Quad wizard view (4-sample equivalent of `trio_wizard`)

---

## Resolved Decisions

1. **Autocomplete cohort status filter** — `Quad` has no `import_status` field. Use `cohort__import_status=ImportStatus.SUCCESS` in `QuadAutocompleteView`, same as how `TrioAutocompleteView` handles it. `get_initial_queryset` in the grid calls `Quad.filter_for_user(self.user)` with no `success_status_only` arg.

2. **Sibling CompHet exclusion** — Known limitation: an unaffected sibling having both comp-het hits isn't excluded. Add a `TODO` comment in `QuadCompHet` pointing to issue #1263. Will be tackled separately.

3. **Quad wizard** — Separate task. Not in scope here.

4. **`analysis_templates_tag` integration** — Add `"quad"` to `single_model_args` and add a proband-as-sample hack parallel to the existing trio one. See Phase 5g above for full details.
