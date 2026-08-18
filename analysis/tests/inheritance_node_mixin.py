"""
Shared inheritance-mode coverage for the pedigree source nodes (TrioNode / QuadNode).

The two nodes run the same inheritance filters over a different sample layout, so the
assertions below are identical - only the zygosity strings, the enum and the
"require parent zygosity" field name change. Subclasses supply those via the class
attributes and hooks at the top of InheritanceNodeTestsMixin.

samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
"""
from snpdb.models import Variant
from snpdb.models.models_cohort import CohortGenotype
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

# attribute, chrom, position - the zygosity string per pattern comes from the subclass
SHARED_VARIANTS = [
    ("recessive_v", "3", 1000),
    ("denovo_v", "3", 2000),
    ("control_v", "3", 4000),
    ("unknown_father_v", "3", 5000),
    ("xlinked_v", "X", 1000),
    ("xlinked_unknown_mother_v", "X", 2000),
]


def make_cohort_genotype(cgc, variant, samples_zygosity: str):
    n = len(samples_zygosity)
    CohortGenotype.objects.create(
        collection=cgc,
        variant=variant,
        ref_count=samples_zygosity.count('R'),
        het_count=samples_zygosity.count('E'),
        hom_count=samples_zygosity.count('O'),
        samples_zygosity=samples_zygosity,
        samples_allele_depth=[20] * n,
        samples_allele_frequency=[100] * n,
        samples_read_depth=[30] * n,
        samples_genotype_quality=[30] * n,
        samples_phred_likelihood=[0] * n,
    )


class InheritanceNodeTestsMixin:
    INHERITANCE = None  # TrioInheritance / QuadInheritance
    NODE_CLASS = None  # TrioNode / QuadNode
    REQUIRE_PARENT_ZYGOSITY = None  # field name - "require_zygosity" / "require_parent_zygosity"
    # samples_zygosity per SHARED_VARIANTS attribute, in the subclass's own sample order
    SHARED_ZYGOSITIES = {}

    @classmethod
    def create_shared_variants(cls, cgc):
        """ Call from setUpTestData once cls.cgc exists - sets recessive_v, denovo_v, ... on cls """
        for attribute, chrom, position in SHARED_VARIANTS:
            variant = slowly_create_test_variant(chrom, position, "A", "T", cls.grch37)
            make_cohort_genotype(cgc, variant, cls.SHARED_ZYGOSITIES[attribute])
            setattr(cls, attribute, variant)

    def _make_node(self, inheritance, **kwargs):
        raise NotImplementedError()

    def _get_inheritance_errors(self, inheritance):
        raise NotImplementedError()

    def _filter_variants(self, node):
        """Apply node's inheritance filter against real Variant data; return matching PKs."""
        arg_q_dict = node._get_node_arg_q_dict()
        cgc = node._get_cohort().cohort_genotype_collection
        qs = Variant.objects.annotate(**cgc.get_annotation_kwargs())
        for alias in (cgc.cohortgenotype_alias, None):
            for q in arg_q_dict.get(alias, {}).values():
                qs = qs.filter(q)
        return set(qs.values_list('pk', flat=True))

    # ── Recessive ─────────────────────────────────────────────────────────────

    def test_recessive_matches_recessive_variant(self):
        node = self._make_node(self.INHERITANCE.RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_recessive_excludes_other_variants(self):
        node = self._make_node(self.INHERITANCE.RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── Denovo ────────────────────────────────────────────────────────────────

    def test_denovo_matches_denovo_variant(self):
        node = self._make_node(self.INHERITANCE.DENOVO)
        self.assertIn(self.denovo_v.pk, self._filter_variants(node))

    def test_denovo_excludes_other_variants(self):
        node = self._make_node(self.INHERITANCE.DENOVO)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── X-linked recessive ────────────────────────────────────────────────────

    def test_xlinked_matches_x_chromosome_variant(self):
        node = self._make_node(self.INHERITANCE.XLINKED_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_xlinked_excludes_autosomal_variants(self):
        node = self._make_node(self.INHERITANCE.XLINKED_RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.denovo_v.pk, ids)

    # ── require parent zygosity ───────────────────────────────────────────────

    def test_require_parent_zygosity_true_excludes_unknown_parental_zygosity(self):
        node = self._make_node(self.INHERITANCE.RECESSIVE, **{self.REQUIRE_PARENT_ZYGOSITY: True})
        self.assertNotIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_require_parent_zygosity_false_allows_unknown_parental_zygosity(self):
        node = self._make_node(self.INHERITANCE.RECESSIVE, **{self.REQUIRE_PARENT_ZYGOSITY: False})
        self.assertIn(self.unknown_father_v.pk, self._filter_variants(node))

    # ── Node type ─────────────────────────────────────────────────────────────

    def test_compound_het_requires_one_parent_input(self):
        node = self._make_node(self.INHERITANCE.COMPOUND_HET)
        self.assertEqual(node.min_inputs, 1)
        self.assertEqual(node.max_inputs, 1)

    def test_simple_modes_are_source_nodes(self):
        for mode in [self.INHERITANCE.RECESSIVE, self.INHERITANCE.DENOVO,
                     self.INHERITANCE.DOMINANT, self.INHERITANCE.XLINKED_RECESSIVE]:
            node = self._make_node(mode)
            self.assertEqual(node.max_inputs, 0, f"{mode} should be a source node")

    def test_clone_matches_same_variants(self):
        node = self._make_node(self.INHERITANCE.RECESSIVE)
        clone = node.save_clone()
        self.assertEqual(self._filter_variants(node), self._filter_variants(clone))

    # ── ALL_RECESSIVE (AR ∪ XLR) ──────────────────────────────────────────────

    def test_all_recessive_matches_autosomal_recessive_variant(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_all_recessive_matches_xlinked_variant(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_all_recessive_excludes_unrelated_variants(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    def test_all_recessive_require_parent_zygosity_excludes_unknown_father_on_autosome(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE, **{self.REQUIRE_PARENT_ZYGOSITY: True})
        self.assertNotIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_all_recessive_no_require_parent_zygosity_includes_unknown_father_on_autosome(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE, **{self.REQUIRE_PARENT_ZYGOSITY: False})
        self.assertIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_all_recessive_require_parent_zygosity_excludes_unknown_mother_on_xlr(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE, **{self.REQUIRE_PARENT_ZYGOSITY: True})
        self.assertNotIn(self.xlinked_unknown_mother_v.pk, self._filter_variants(node))

    def test_all_recessive_no_require_parent_zygosity_includes_unknown_mother_on_xlr(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE, **{self.REQUIRE_PARENT_ZYGOSITY: False})
        self.assertIn(self.xlinked_unknown_mother_v.pk, self._filter_variants(node))

    def test_all_recessive_is_source_node(self):
        node = self._make_node(self.INHERITANCE.ALL_RECESSIVE)
        self.assertEqual(node.max_inputs, 0)

    # ── ANY_AFFECTED ──────────────────────────────────────────────────────────

    def test_any_affected_is_source_node(self):
        node = self._make_node(self.INHERITANCE.ANY_AFFECTED)
        self.assertEqual(node.max_inputs, 0)

    def test_any_affected_always_valid_no_errors(self):
        self.assertEqual(self._get_inheritance_errors(self.INHERITANCE.ANY_AFFECTED), [])

    # ── Zygosity table ────────────────────────────────────────────────────────

    def test_zygosity_table_all_recessive_two_line_cells(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        self.assertIn('AR:', data[self.INHERITANCE.ALL_RECESSIVE]['mother'])
        self.assertIn('XLR:', data[self.INHERITANCE.ALL_RECESSIVE]['mother'])

    def test_zygosity_table_all_recessive_other_filters_mentions_chr_x(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        self.assertIn('Chr X', data[self.INHERITANCE.ALL_RECESSIVE]['other_filters_mother'])

    def test_zygosity_table_any_affected_proband_has_variant(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        self.assertTrue(data[self.INHERITANCE.ANY_AFFECTED]['proband'])

    def test_zygosity_table_xlinked_has_chr_x_other_filter(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        self.assertEqual(data[self.INHERITANCE.XLINKED_RECESSIVE]['other_filters_mother'], "Chr X only")

    def test_zygosity_table_recessive_has_no_other_filter_keys(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        keys = [k for k in data[self.INHERITANCE.RECESSIVE] if k.startswith('other_filters_')]
        self.assertEqual(keys, [])

    def test_zygosity_table_compound_het_has_gene_constraint_in_other_filters(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        self.assertIn("gene", data[self.INHERITANCE.COMPOUND_HET]['other_filters_mother'])

    def test_zygosity_table_compound_het_has_no_note_key(self):
        data = self.NODE_CLASS.get_zygosity_table_data()
        self.assertNotIn('note', data[self.INHERITANCE.COMPOUND_HET])
