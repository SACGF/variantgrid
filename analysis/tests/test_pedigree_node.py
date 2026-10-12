"""
PedigreeNode inheritance filters against real CohortGenotype rows: which zygosities an unaffected
member may carry under each model, and that a member of unknown affection is unconstrained.

samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, PedigreeNode
from annotation.fake_data import get_fake_annotation_version
from pedigree.models import CohortSamplePedFileRecord, PedFileRecord, PedigreeInheritance
from snpdb.fake_data import create_fake_pedigree, make_cohort_genotype
from snpdb.models import GenomeBuild, Variant
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPedigreeNodeInheritance(TestCase):
    """ Pedigree over create_fake_cohort: [proband (affected), mother (unaffected), father (unknown)] """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser_pedigree_node')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.pedigree = create_fake_pedigree(user, cls.grch37)
        family = cls.pedigree.ped_file_family
        cohort = cls.pedigree.cohort
        for sample_name, affection in [("proband", True), ("mother", False), ("father", None)]:
            record = PedFileRecord.objects.create(family=family, sample=sample_name, affection=affection)
            cohort_sample = cohort.cohortsample_set.get(sample__name=sample_name)
            CohortSamplePedFileRecord.objects.create(pedigree=cls.pedigree, cohort_sample=cohort_sample,
                                                     ped_file_record=record)
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cohort)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        # [proband, mother, father]
        cls.dominant_mother_ref_v = cls._make_variant(1000, "ERR")
        cls.dominant_mother_missing_v = cls._make_variant(2000, "E.R")
        cls.dominant_mother_het_v = cls._make_variant(3000, "EER")
        cls.dominant_father_hom_v = cls._make_variant(4000, "ERO")
        cls.recessive_mother_het_v = cls._make_variant(5000, "OER")
        cls.recessive_mother_ref_v = cls._make_variant(6000, "ORR")
        cls.recessive_mother_hom_v = cls._make_variant(7000, "OOR")
        cls.recessive_father_hom_v = cls._make_variant(8000, "OEO")

    @classmethod
    def _make_variant(cls, position, zygosity):
        variant = slowly_create_test_variant("3", position, "A", "T", cls.grch37)
        make_cohort_genotype(cls.cgc, variant, zygosity)
        return variant

    def _filter_variants(self, inheritance_model):
        node = PedigreeNode.objects.create(analysis=self.analysis, pedigree=self.pedigree,
                                           inheritance_model=inheritance_model)
        arg_q_dict = node._get_node_arg_q_dict()
        qs = Variant.objects.annotate(**self.cgc.get_annotation_kwargs())
        for alias in (self.cgc.cohortgenotype_alias, None):
            for q in arg_q_dict.get(alias, {}).values():
                qs = qs.filter(q)
        return set(qs.values_list('pk', flat=True))

    def test_dominant(self):
        pks = self._filter_variants(PedigreeInheritance.AUTOSOMAL_DOMINANT)
        self.assertIn(self.dominant_mother_ref_v.pk, pks)
        self.assertIn(self.dominant_mother_missing_v.pk, pks)
        self.assertNotIn(self.dominant_mother_het_v.pk, pks)

    def test_dominant_unknown_affection_is_unconstrained(self):
        pks = self._filter_variants(PedigreeInheritance.AUTOSOMAL_DOMINANT)
        self.assertIn(self.dominant_father_hom_v.pk, pks)

    def test_recessive(self):
        pks = self._filter_variants(PedigreeInheritance.AUTOSOMAL_RECESSIVE)
        self.assertIn(self.recessive_mother_het_v.pk, pks)
        self.assertIn(self.recessive_mother_ref_v.pk, pks)
        self.assertNotIn(self.recessive_mother_hom_v.pk, pks)
        self.assertNotIn(self.dominant_mother_ref_v.pk, pks)

    def test_recessive_unknown_affection_is_unconstrained(self):
        pks = self._filter_variants(PedigreeInheritance.AUTOSOMAL_RECESSIVE)
        self.assertIn(self.recessive_father_hom_v.pk, pks)
