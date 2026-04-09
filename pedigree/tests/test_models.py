"""
Tests for pedigree.models — validate(), PedFileRecord.affected, create_automatch_pedigree.
"""
import django.test
from django.contrib.auth.models import User

from patients.models_enums import Sex
from pedigree.models import (
    CohortSamplePedFileRecord,
    PedFile,
    PedFileFamily,
    PedFileRecord,
    Pedigree,
    create_automatch_pedigree,
    validate,
)
from snpdb.models import ImportStatus
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort


def _make_family(user, name="FAM001"):
    ped_file = PedFile.objects.create(name="test", user=user, import_status=ImportStatus.SUCCESS)
    return PedFileFamily.objects.create(ped_file=ped_file, name=name)


def _make_record(family, sample, sex=None, affection=None, father=None, mother=None):
    return PedFileRecord.objects.create(
        family=family, sample=sample, sex=sex, affection=affection,
        father=father, mother=mother,
    )


class TestValidate(django.test.TestCase):

    def setUp(self):
        self.user = User.objects.create_user('validate_test')
        self.family = _make_family(self.user)

    def test_empty_records_fails(self):
        errors = validate([])
        self.assertTrue(errors)

    def test_single_unaffected_individual_fails(self):
        r = _make_record(self.family, 'proband', sex=Sex.MALE, affection=False)
        errors = validate([r])
        self.assertTrue(any('affected' in e.lower() for e in errors))

    def test_unknown_affection_none_not_counted_as_affected(self):
        # affection=None (unknown) must not satisfy the "at least 1 affected" requirement
        r = _make_record(self.family, 'proband', sex=Sex.MALE, affection=None)
        errors = validate([r])
        self.assertTrue(any('affected' in e.lower() for e in errors))

    def test_father_with_female_sex_fails(self):
        father = _make_record(self.family, 'father', sex=Sex.FEMALE, affection=False)
        proband = _make_record(self.family, 'proband', sex=Sex.MALE, affection=True, father=father)
        errors = validate([proband, father])
        self.assertTrue(any('not male' in e for e in errors))

    def test_mother_with_male_sex_fails(self):
        mother = _make_record(self.family, 'mother', sex=Sex.MALE, affection=False)
        proband = _make_record(self.family, 'proband', sex=Sex.FEMALE, affection=True, mother=mother)
        errors = validate([proband, mother])
        self.assertTrue(any('not female' in e for e in errors))

    def test_father_with_null_sex_fails(self):
        # Unknown sex on a parent is not sufficient — the role requires a known matching sex
        father = _make_record(self.family, 'father', sex=None, affection=False)
        proband = _make_record(self.family, 'proband', sex=Sex.MALE, affection=True, father=father)
        errors = validate([proband, father])
        self.assertTrue(errors)


class TestPedFileRecordAffected(django.test.TestCase):

    def setUp(self):
        user = User.objects.create_user('affected_test')
        self.family = _make_family(user)

    def test_none_should_not_return_unaffected(self):
        # BUG-2: affection=None (unknown) returns "unaffected" due to `if self.affection:`
        r = _make_record(self.family, 'proband', affection=None)
        self.assertNotEqual(r.affected, "unaffected",
                            "affection=None (unknown) should not display as 'unaffected'")


class TestPedFileFamilyIsValid(django.test.TestCase):

    def setUp(self):
        self.user = User.objects.create_user('family_valid_test')

    def test_family_with_no_records_is_invalid(self):
        family = _make_family(self.user)
        self.assertFalse(family.is_valid)
        self.assertTrue(family.errors)


class TestCreateAutomatchPedigree(django.test.TestCase):

    def setUp(self):
        self.user = User.objects.create_user('automatch_test')
        self.grch37 = GenomeBuild.get_name_or_alias('GRCh37')

    def test_matching_samples_are_linked(self):
        cohort = create_fake_cohort(self.user, self.grch37)
        family = _make_family(self.user)
        _make_record(family, 'proband', affection=True)
        _make_record(family, 'mother', affection=False)
        _make_record(family, 'father', affection=False)

        pedigree = create_automatch_pedigree(self.user, family, cohort)

        links = CohortSamplePedFileRecord.objects.filter(pedigree=pedigree)
        self.assertEqual(links.count(), 3)

    def test_no_matching_samples_creates_empty_pedigree(self):
        cohort = create_fake_cohort(self.user, self.grch37)
        family = _make_family(self.user)
        _make_record(family, 'no_match_sample', affection=True)

        pedigree = create_automatch_pedigree(self.user, family, cohort)

        links = CohortSamplePedFileRecord.objects.filter(pedigree=pedigree)
        self.assertEqual(links.count(), 0)

    def test_partial_match_links_only_matching(self):
        cohort = create_fake_cohort(self.user, self.grch37)
        family = _make_family(self.user)
        _make_record(family, 'proband', affection=True)
        _make_record(family, 'no_match', affection=False)

        pedigree = create_automatch_pedigree(self.user, family, cohort)

        links = CohortSamplePedFileRecord.objects.filter(pedigree=pedigree)
        self.assertEqual(links.count(), 1)


class TestPedigreeGetSamples(django.test.TestCase):

    def setUp(self):
        self.user = User.objects.create_user('getsamples_test')
        self.grch37 = GenomeBuild.get_name_or_alias('GRCh37')

    def _setup_pedigree(self):
        cohort = create_fake_cohort(self.user, self.grch37)
        family = _make_family(self.user)
        affected_rec   = _make_record(family, 'proband', affection=True)
        unaffected_rec = _make_record(family, 'mother',  affection=False)
        unknown_rec    = _make_record(family, 'father',  affection=None)

        pedigree = Pedigree.objects.create(
            user=self.user, name="test ped", cohort=cohort, ped_file_family=family
        )
        cs_proband = cohort.cohortsample_set.get(sample__name='proband')
        cs_mother  = cohort.cohortsample_set.get(sample__name='mother')
        cs_father  = cohort.cohortsample_set.get(sample__name='father')
        CohortSamplePedFileRecord.objects.create(
            pedigree=pedigree, cohort_sample=cs_proband, ped_file_record=affected_rec)
        CohortSamplePedFileRecord.objects.create(
            pedigree=pedigree, cohort_sample=cs_mother, ped_file_record=unaffected_rec)
        CohortSamplePedFileRecord.objects.create(
            pedigree=pedigree, cohort_sample=cs_father, ped_file_record=unknown_rec)
        return pedigree

    def test_get_samples_affected_only(self):
        pedigree = self._setup_pedigree()
        names = list(pedigree.get_samples(affected=True).values_list('name', flat=True))
        self.assertEqual(names, ['proband'])

    def test_get_samples_unaffected_only(self):
        pedigree = self._setup_pedigree()
        names = list(pedigree.get_samples(affected=False).values_list('name', flat=True))
        self.assertEqual(names, ['mother'])

    def test_get_samples_affected_excludes_unknown(self):
        # affection=None must not appear when filtering for affected=True
        pedigree = self._setup_pedigree()
        names = list(pedigree.get_samples(affected=True).values_list('name', flat=True))
        self.assertNotIn('father', names)
