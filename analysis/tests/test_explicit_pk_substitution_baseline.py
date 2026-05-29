"""
Baseline regression test for the issue #546 explicit-PK substitution optimisation.

This locks in the CURRENT AnalysisNode query behaviour (before the small-parent
pk-substitution change) so that after the optimisation lands we can re-run it and
confirm node counts, result PK sets, and grid genotype display values are byte-for-byte
identical.

It deliberately exercises the parts the substitution touches or risks:

  * CohortGenotype data split across UNCOMMON + COMMON partitions (the rare/common
    split, #1119) so the FilteredRelation join must union both partitions.
  * A SampleNode source whose variants span both partitions.
  * ZygosityNode and AlleleFrequencyNode children that filter on the cohortgenotype
    alias (the only filter nodes that consume the alias for counts) - the case the
    substitution must keep correct.
  * The grid genotype display columns (zygosity, allele depth, read depth, allele
    frequency) resolved through the alias join, for variants in both partitions.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.enums import NodeStatus, ZygosityNodeZygosity
from analysis.models.nodes.filters.allele_frequency_node import AlleleFrequencyNode
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.filters.zygosity_node import ZygosityNode
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import AnnotationRun
from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from patients.models_enums import Zygosity
from snpdb.models import GenomeBuild, Variant
from snpdb.models.models_cohort import (
    CohortGenotype, CohortGenotypeCollection, CohortGenotypeCommonFilterVersion,
)
from snpdb.models.models_enums import CohortGenotypeCollectionType
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestExplicitPkSubstitutionBaseline(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_546_baseline')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.trio = create_fake_trio(cls.user, cls.grch37)
        cls.cohort = cls.trio.cohort  # 3 samples: proband(0), mother(1), father(2)
        cls.proband = cls.cohort.cohortsample_set.get(sample__name='proband').sample

        cls.uncommon_cgc = CohortGenotypeCollection.objects.get(
            cohort=cls.cohort, cohort_version=cls.cohort.version,
            collection_type=CohortGenotypeCollectionType.UNCOMMON)

        # Build the COMMON partition + link, so the FilteredRelation must union both
        # partitions (rare/common split, #1119).
        cfv = CohortGenotypeCommonFilterVersion.objects.create(
            gnomad_version="2.1.1", gnomad_af_min=0.05, genome_build=cls.grch37)
        cls.common_cgc = CohortGenotypeCollection.objects.create(
            cohort=cls.cohort, cohort_version=cls.cohort.version,
            num_samples=cls.uncommon_cgc.num_samples,
            collection_type=CohortGenotypeCollectionType.COMMON,
            common_filter=cfv)
        cls.uncommon_cgc.common_collection = cls.common_cgc
        cls.uncommon_cgc.save()

        # Variants. proband is packed index 0; zygosity strings are 3 chars (proband, mother, father).
        # Encoding: E=HET, O=HOM_ALT, R=HOM_REF, .=MISSING.
        # UNCOMMON (rare) partition:
        cls.v_u_het = slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)
        cls._add_genotype(cls.uncommon_cgc, cls.v_u_het, "E..", ad=15, af=0.10, dp=40)
        cls.v_u_hom = slowly_create_test_variant("1", 2000, "A", "T", cls.grch37)
        cls._add_genotype(cls.uncommon_cgc, cls.v_u_hom, "O..", ad=25, af=0.50, dp=50)
        cls.v_u_ref = slowly_create_test_variant("1", 3000, "A", "T", cls.grch37)  # HOM_REF -> not in SampleNode default
        cls._add_genotype(cls.uncommon_cgc, cls.v_u_ref, "R..", ad=35, af=0.95, dp=60)
        # COMMON partition:
        cls.v_c_het = slowly_create_test_variant("1", 4000, "A", "T", cls.grch37)
        cls._add_genotype(cls.common_cgc, cls.v_c_het, "E..", ad=45, af=0.80, dp=70)
        cls.v_c_hom = slowly_create_test_variant("1", 5000, "A", "T", cls.grch37)
        cls._add_genotype(cls.common_cgc, cls.v_c_hom, "O..", ad=55, af=0.90, dp=80)

        # SampleNode default keeps HET + HOM_ALT (zygosity_het/zygosity_hom default True).
        cls.source_pks = {cls.v_u_het.pk, cls.v_u_hom.pk, cls.v_c_het.pk, cls.v_c_hom.pk}

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    @classmethod
    def _add_genotype(cls, cgc, variant, samples_zygosity, ad, af, dp):
        """ Insert a CohortGenotype row into the given collection's partition table
            (mirrors production - ORM inserts target the right partition). Per-sample
            arrays put the proband at index 0. """
        n = len(samples_zygosity)
        old_db_table = CohortGenotype._meta.db_table
        try:
            CohortGenotype._meta.db_table = cgc.get_partition_table()
            CohortGenotype.objects.create(
                collection=cgc, variant=variant,
                ref_count=samples_zygosity.count('R'),
                het_count=samples_zygosity.count('E'),
                hom_count=samples_zygosity.count('O'),
                unk_count=samples_zygosity.count('U'),
                samples_zygosity=samples_zygosity,
                samples_allele_depth=[ad] + [0] * (n - 1),
                samples_allele_frequency=[af] + [0.0] * (n - 1),
                samples_read_depth=[dp] + [0] * (n - 1),
                samples_genotype_quality=[30] * n,
                samples_phred_likelihood=[0] * n,
            )
        finally:
            CohortGenotype._meta.db_table = old_db_table

    @staticmethod
    def _ready(node):
        """ Make a node 'ready' with its count set, as the load pipeline would, so child
            nodes can compose it (and so the future pk-substitution path can read count). """
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()
        return node

    def _source_node(self):
        return self._ready(SampleNode.objects.create(analysis=self.analysis, sample=self.proband))

    def _child(self, node):
        node.add_parent(self._source_node())
        node._cached_parents = None  # clear stale cache from create()'s save()
        node.save()
        return node

    # ── Source spanning both partitions ───────────────────────────────────────

    def test_source_node_spans_both_partitions(self):
        sample_node = self._source_node()
        pks = set(sample_node.get_queryset().values_list("pk", flat=True))
        self.assertEqual(pks, self.source_pks)
        self.assertEqual(sample_node.count, 4)

    # ── ZygosityNode child (filters on the cohortgenotype alias) ───────────────

    def test_zygosity_node_het(self):
        zyg = ZygosityNode.objects.create(analysis=self.analysis, sample=self.proband,
                                          zygosity=ZygosityNodeZygosity.HET)
        zyg = self._child(zyg)
        pks = set(zyg.get_queryset().values_list("pk", flat=True))
        self.assertEqual(pks, {self.v_u_het.pk, self.v_c_het.pk})

    def test_zygosity_node_hom_alt(self):
        zyg = ZygosityNode.objects.create(analysis=self.analysis, sample=self.proband,
                                          zygosity=ZygosityNodeZygosity.HOM_ALT)
        zyg = self._child(zyg)
        pks = set(zyg.get_queryset().values_list("pk", flat=True))
        self.assertEqual(pks, {self.v_u_hom.pk, self.v_c_hom.pk})

    # ── AlleleFrequencyNode child (filters on the cohortgenotype alias) ────────

    def test_allele_frequency_node_range(self):
        af = AlleleFrequencyNode.objects.create(analysis=self.analysis, sample=self.proband)
        af = self._child(af)
        # AF stored as 0-1 fractions (vcf.allele_frequency_percent=False). Keep 0.4 <= AF <= 0.85.
        af.nodeallelefrequencyfilter.nodeallelefrequencyrange_set.update(min=0.4, max=0.85)
        pks = set(af.get_queryset().values_list("pk", flat=True))
        # Of the source set, AF: v_u_het=0.10(no), v_u_hom=0.50(yes), v_c_het=0.80(yes), v_c_hom=0.90(no)
        self.assertEqual(pks, {self.v_u_hom.pk, self.v_c_het.pk})

    # ── Grid genotype display columns (alias join), across both partitions ─────

    def test_grid_genotype_display_columns(self):
        sample_node = self._source_node()
        zyg_field = self.proband.get_cohort_genotype_alias_and_field("zygosity")[1]
        ad_field = self.proband.get_cohort_genotype_alias_and_field("allele_depth")[1]
        dp_field = self.proband.get_cohort_genotype_alias_and_field("read_depth")[1]
        af_field = self.proband.get_cohort_genotype_alias_and_field("allele_frequency")[1]

        rows = {r["pk"]: r for r in
                sample_node.get_queryset().values("pk", zyg_field, ad_field, dp_field, af_field)}

        expected = {
            self.v_u_het.pk: (Zygosity.HET, 15, 40, 0.10),
            self.v_u_hom.pk: (Zygosity.HOM_ALT, 25, 50, 0.50),
            self.v_c_het.pk: (Zygosity.HET, 45, 70, 0.80),   # from the COMMON partition
            self.v_c_hom.pk: (Zygosity.HOM_ALT, 55, 80, 0.90),
        }
        self.assertEqual(set(rows), set(expected))
        for pk, (zyg, ad, dp, af) in expected.items():
            self.assertEqual(rows[pk][zyg_field], zyg, f"zygosity for {pk}")
            self.assertEqual(rows[pk][ad_field], ad, f"allele_depth for {pk}")
            self.assertEqual(rows[pk][dp_field], dp, f"read_depth for {pk}")
            self.assertAlmostEqual(rows[pk][af_field], af, places=5, msg=f"allele_frequency for {pk}")

    # ── MergeNode over small parents (the other caller of the shared helper) ───

    def _ready_zygosity_child(self, zygosity):
        """ A ZygosityNode on the proband, with its own ready SampleNode source, made ready
            so MergeNode can read its count and compose it. """
        zyg = ZygosityNode.objects.create(analysis=self.analysis, sample=self.proband,
                                          zygosity=zygosity)
        zyg.add_parent(self._source_node())
        zyg._cached_parents = None  # clear stale cache from create()'s save()
        zyg.save()
        # Reload fresh so parents_changed is cleared - otherwise the _ready() save() would
        # bump_version() again and reset status back to DIRTY.
        zyg = ZygosityNode.objects.get(pk=zyg.pk)
        return self._ready(zyg)

    def test_merge_node_equivalent_with_gate_on_vs_off(self):
        """ #546 - the shared small-parent helper also feeds the MergeNode path. A MergeNode
            merging two small ZygosityNode parents (HET + HOM_ALT) must yield the identical PK
            set whether the explicit-PK substitution gate is on (each parent collapses to a
            literal Q(pk__in=[...])) or off (each parent runs as a full subquery and is merged).
            The result is the union of the two parents. """
        def merge_pks():
            het = self._ready_zygosity_child(ZygosityNodeZygosity.HET)
            hom = self._ready_zygosity_child(ZygosityNodeZygosity.HOM_ALT)
            merge = MergeNode.objects.create(analysis=self.analysis)
            merge.add_parent(het)
            merge.add_parent(hom)
            merge._cached_parents = None
            merge.save()
            return set(merge.get_queryset().values_list("pk", flat=True))

        with override_settings(ANALYSIS_NODE_STORE_ID_SIZE_MAX=0):
            pks_gate_off = merge_pks()
        with override_settings(ANALYSIS_NODE_STORE_ID_SIZE_MAX=1000):
            pks_gate_on = merge_pks()

        expected = {self.v_u_het.pk, self.v_c_het.pk, self.v_u_hom.pk, self.v_c_hom.pk}
        self.assertEqual(pks_gate_off, expected)
        self.assertEqual(pks_gate_on, pks_gate_off)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestRarePopulationNodePruningEquivalence(TestCase):
    """ Issue #546 - the explicit-PK substitution must be equivalent to the non-substituted path
        even when a downstream rare PopulationNode prunes the parent cohortgenotype join to the
        UNCOMMON partition (common_variants=False, #1119).

        Without substitution (gate off): the rare PopulationNode pushes common_variants=False up the
        chain, so the parent's join only spans the UNCOMMON partition - the COMMON variant never
        enters the query.

        With substitution (gate on): get_parent_pks materialises the parent standalone with
        common_variants=True (BOTH partitions), so the COMMON variant IS in the literal pk__in list -
        but the downstream rare gnomAD filter then removes it via VariantAnnotation (common-partition
        variants are common in gnomAD, so their AF is above the rare cutoff). The final PK set matches.

        This pins that both gates yield the identical result PK set. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_546_rare_pop')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.annotation_run = AnnotationRun.objects.create()

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        # Variant queryset joins variantannotation on the analysis's resolved version - so the
        # rare gnomAD filter only matches rows created against that exact version.
        cls.vav = cls.analysis.annotation_version.variant_annotation_version

        cls.trio = create_fake_trio(cls.user, cls.grch37)
        cls.cohort = cls.trio.cohort
        cls.proband = cls.cohort.cohortsample_set.get(sample__name='proband').sample

        cls.uncommon_cgc = CohortGenotypeCollection.objects.get(
            cohort=cls.cohort, cohort_version=cls.cohort.version,
            collection_type=CohortGenotypeCollectionType.UNCOMMON)

        # COMMON partition + link, so the parent join must union both partitions unless pruned.
        cfv = CohortGenotypeCommonFilterVersion.objects.create(
            gnomad_version="2.1.1", gnomad_af_min=0.05, genome_build=cls.grch37)
        cls.common_cgc = CohortGenotypeCollection.objects.create(
            cohort=cls.cohort, cohort_version=cls.cohort.version,
            num_samples=cls.uncommon_cgc.num_samples,
            collection_type=CohortGenotypeCollectionType.COMMON,
            common_filter=cfv)
        cls.uncommon_cgc.common_collection = cls.common_cgc
        cls.uncommon_cgc.save()

        # UNCOMMON (rare) variant: HET, gnomAD AF below the 1% cutoff -> kept by rare PopulationNode.
        cls.v_rare = slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)
        cls._add_genotype(cls.uncommon_cgc, cls.v_rare, "E..")
        cls._add_annotation(cls.v_rare, gnomad_af=0.001, af_1kg=0.001, af_uk10k=0.001, max_af=0.001)

        # COMMON variant: HET, gnomAD AF well above the cutoff -> removed by rare PopulationNode.
        cls.v_common = slowly_create_test_variant("1", 4000, "A", "T", cls.grch37)
        cls._add_genotype(cls.common_cgc, cls.v_common, "E..")
        cls._add_annotation(cls.v_common, gnomad_af=0.20, af_1kg=0.20, af_uk10k=0.20, max_af=0.20)

    @classmethod
    def _add_genotype(cls, cgc, variant, samples_zygosity):
        n = len(samples_zygosity)
        old_db_table = CohortGenotype._meta.db_table
        try:
            CohortGenotype._meta.db_table = cgc.get_partition_table()
            CohortGenotype.objects.create(
                collection=cgc, variant=variant,
                ref_count=samples_zygosity.count('R'),
                het_count=samples_zygosity.count('E'),
                hom_count=samples_zygosity.count('O'),
                unk_count=samples_zygosity.count('U'),
                samples_zygosity=samples_zygosity,
                samples_allele_depth=[20] + [0] * (n - 1),
                samples_allele_frequency=[0.5] + [0.0] * (n - 1),
                samples_read_depth=[40] + [0] * (n - 1),
                samples_genotype_quality=[30] * n,
                samples_phred_likelihood=[0] * n,
            )
        finally:
            CohortGenotype._meta.db_table = old_db_table

    @classmethod
    def _add_annotation(cls, variant, gnomad_af, af_1kg, af_uk10k, max_af):
        """ Insert into the version's VariantAnnotation partition table directly (INHERITS
            partitioning doesn't route base-table inserts), so the annotation_version SQL
            transformer that rewrites the query to that partition can see the row. """
        partition_table = cls.vav.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        old_db_table = VariantAnnotation._meta.db_table
        try:
            VariantAnnotation._meta.db_table = partition_table
            VariantAnnotation.objects.create(
                version=cls.vav, variant=variant, annotation_run=cls.annotation_run,
                gnomad_af=gnomad_af, af_1kg=af_1kg, af_uk10k=af_uk10k, max_af=max_af,
                predictions_num_pathogenic=0, predictions_num_benign=0)
        finally:
            VariantAnnotation._meta.db_table = old_db_table

    def _ready_source(self):
        node = SampleNode.objects.create(analysis=self.analysis, sample=self.proband)
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()
        return node

    def _rare_population_child(self):
        # percent <= 5 with gnomad_af -> _has_common_variants() is False -> pushes common_variants=False up.
        pop = PopulationNode.objects.create(
            analysis=self.analysis, percent=1.0,
            gnomad_af=True, gnomad_popmax_af=False, af_1kg=True, af_uk10k=True, topmed_af=False)
        pop.add_parent(self._ready_source())
        pop._cached_parents = None
        pop.save()
        return pop

    def _result_pks(self):
        return set(self._rare_population_child().get_queryset().values_list("pk", flat=True))

    def test_rare_population_node_result_equivalent_with_gate_on_vs_off(self):
        with override_settings(ANALYSIS_NODE_STORE_ID_SIZE_MAX=0):
            pks_gate_off = self._result_pks()
        with override_settings(ANALYSIS_NODE_STORE_ID_SIZE_MAX=1000):
            pks_gate_on = self._result_pks()

        # The rare variant survives both paths; the common variant is excluded by both.
        self.assertEqual(pks_gate_off, {self.v_rare.pk})
        self.assertEqual(pks_gate_on, pks_gate_off)
