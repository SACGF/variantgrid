"""
    X-linked recessive proband sex check + revealing nodes a template run hid (sapath#439)
"""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from guardian.shortcuts import assign_perm

from analysis.forms.forms_nodes import TrioNodeForm
from analysis.models import Analysis, QuadNode, TrioNode
from analysis.models.enums import AnalysisTemplateType, QuadInheritance, TrioInheritance, TrioSample
from annotation.fake_annotation import get_fake_annotation_version
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient
from patients.models_enums import Sex
from snpdb.models import (
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
    SampleStatsCodeVersion,
    Trio,
)
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort, create_fake_quad, create_fake_trio


class TestXLinkedRecessiveProbandSex(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_xlinked_sex')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.trio = create_fake_trio(cls.user, cls.genome_build)
        cls.trio.mother_affected = False  # XLR requires an unaffected mother
        cls.trio.save()
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.trio.cohort)
        cls.code_version = SampleStatsCodeVersion.objects.create(name="SampleStats", version=1, code_git_hash="test")

    def setUp(self):
        self.proband_sample = self.trio.proband.sample
        patient = Patient.objects.create(family_code="fam", sex=Sex.FEMALE)
        self.proband_sample.patient = patient
        self.proband_sample.save()
        self.trio.proband_sex = None
        self.trio.save()

    def _set_detected_male(self):
        CohortGenotypeStats.objects.create(cohort_genotype_collection=self.cgc, sample=self.proband_sample,
                                           code_version=self.code_version,
                                           x_hom_count=200, x_het_count=10)
        self.trio.refresh_from_db()  # clear cached detected_sex

    def _get_errors(self):
        return TrioNode.get_trio_inheritance_errors(self.trio, TrioInheritance.XLINKED_RECESSIVE)

    def test_recorded_female_errors(self):
        """ Historical behaviour - the patient record decides unless the trio says otherwise """
        errors = self._get_errors()
        self.assertEqual(1, len(errors))
        self.assertNotIn("chrX genotypes", errors[0])

    def test_recorded_female_with_detected_male_suggests_setting_proband_sex(self):
        self._set_detected_male()
        errors = self._get_errors()
        self.assertEqual(1, len(errors))
        self.assertIn("chrX genotypes detected male", errors[0])

    def test_trio_proband_sex_overrides_patient_record(self):
        """ Scientist resolved the mismatch in the trio wizard """
        self._set_detected_male()
        self.trio.proband_sex = Sex.MALE
        self.trio.save()
        self.assertEqual([], self._get_errors())


class TestTrioWizardProbandSex(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_trio_wizard_sex')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.samples = [cs.sample for cs in cls.cohort.cohortsample_set.order_by("sort_order")]
        for sample in cls.samples[1:]:
            assign_permission_to_user_and_groups(cls.user, sample)
        cls.url = reverse('trio_wizard', kwargs={"cohort_id": cls.cohort.pk,
                                                 "sample1_id": cls.samples[0].pk,
                                                 "sample2_id": cls.samples[1].pk,
                                                 "sample3_id": cls.samples[2].pk})

    def _post(self, proband_sex):
        self.client.force_login(self.user)
        return self.client.post(self.url, {"sample_1": TrioSample.PROBAND,
                                           "sample_2": TrioSample.MOTHER,
                                           "sample_3": TrioSample.FATHER,
                                           "proband_sex": proband_sex})

    def test_chosen_proband_sex_saved(self):
        self._post(Sex.MALE)
        trio = Trio.objects.get(cohort=self.cohort)
        self.assertEqual(Sex.MALE, trio.effective_proband_sex)

    def test_no_choice_leaves_patient_record_in_charge(self):
        self.samples[0].patient = Patient.objects.create(family_code="fam", sex=Sex.FEMALE)
        self.samples[0].save()
        self._post("")
        trio = Trio.objects.get(cohort=self.cohort)
        self.assertIsNone(trio.proband_sex)
        self.assertEqual(Sex.FEMALE, trio.effective_proband_sex)

    def test_rerunning_wizard_updates_chosen_sex(self):
        self._post(Sex.MALE)
        self._post(Sex.FEMALE)
        trio = Trio.objects.get(cohort=self.cohort)
        self.assertEqual(Sex.FEMALE, trio.effective_proband_sex)


class TestQuadProbandSex(TestCase):
    """ QuadNode shares the X-linked recessive check, so the quad carries its own proband sex """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_quad_sex')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.quad = create_fake_quad(cls.user, cls.genome_build)
        cls.quad.mother_affected = False  # XLR requires an unaffected mother
        cls.quad.save()

    def setUp(self):
        proband_sample = self.quad.proband.sample
        proband_sample.patient = Patient.objects.create(family_code="fam", sex=Sex.FEMALE)
        proband_sample.save()
        self.quad.proband_sex = None
        self.quad.save()

    def _get_errors(self):
        return QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.XLINKED_RECESSIVE)

    def test_recorded_female_errors(self):
        self.assertEqual(1, len(self._get_errors()))

    def test_quad_proband_sex_overrides_patient_record(self):
        self.quad.proband_sex = Sex.MALE
        self.quad.save()
        self.assertEqual([], self._get_errors())


class TestRevealHiddenNodes(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_reveal_owner')[0]
        cls.other_user = User.objects.get_or_create(username='testuser_reveal_other')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.analysis = Analysis(genome_build=cls.genome_build)
        cls.analysis.set_defaults_and_save(cls.user)

    def setUp(self):
        self.node = TrioNode.objects.create(analysis=self.analysis, visible=False)
        self.child = TrioNode.objects.create(analysis=self.analysis, visible=False)
        self.node.add_child(self.child)
        self.sibling = TrioNode.objects.create(analysis=self.analysis, visible=False)
        self.url = reverse('node_reveal_hidden', kwargs={"analysis_id": self.analysis.pk, "node_id": self.node.pk})

    def test_reveals_node_and_descendants(self):
        self.client.force_login(self.user)
        response = self.client.post(self.url)
        self.assertEqual(200, response.status_code)
        self.node.refresh_from_db()
        self.child.refresh_from_db()
        self.sibling.refresh_from_db()
        self.assertTrue(self.node.visible)
        self.assertTrue(self.child.visible)
        self.assertFalse(self.sibling.visible)

        data = response.json()
        self.assertEqual({self.node.pk, self.child.pk}, {n["attributes"]["node_id"] for n in data["nodes"]})
        self.assertEqual([{"source_id": self.node.get_css_id(), "target_id": self.child.get_css_id()}], data["edges"])

    def test_read_only_user_denied(self):
        assign_perm(self.analysis.get_read_perm(), self.other_user, self.analysis)
        self.client.force_login(self.other_user)
        response = self.client.post(self.url)
        self.assertEqual(403, response.status_code)
        self.node.refresh_from_db()
        self.assertFalse(self.node.visible)


class TestIgnoreFieldErrors(TestCase):
    """ ignore_field_errors turns the inheritance checks into warnings so the node can save and run """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_ignore_field_errors')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.trio = create_fake_trio(cls.user, cls.genome_build)  # mother affected - fails X-linked recessive
        cls.analysis = Analysis(genome_build=cls.genome_build)
        cls.analysis.set_defaults_and_save(cls.user)

    def _make_node(self, **kwargs):
        return TrioNode.objects.create(analysis=self.analysis, trio=self.trio,
                                       inheritance=TrioInheritance.XLINKED_RECESSIVE, **kwargs)

    def test_inheritance_error_blocks_node(self):
        node = self._make_node()
        self.assertEqual({"inheritance"}, set(node.get_field_errors()))
        self.assertEqual(1, len(node.get_errors(include_parent_errors=False)))
        self.assertEqual([], node.get_warnings())
        self.assertTrue(node.can_ignore_errors())

    def test_ignored_error_becomes_warning(self):
        node = self._make_node(ignore_field_errors=True)
        self.assertEqual([], node.get_errors(include_parent_errors=False))
        self.assertEqual(1, len(node.get_warnings()))
        self.assertIn("unaffected mother", node.get_warnings()[0])
        self.assertFalse(node.can_ignore_errors())

    def test_other_errors_cannot_be_ignored(self):
        node = TrioNode.objects.create(analysis=self.analysis, trio=None,
                                       inheritance=TrioInheritance.XLINKED_RECESSIVE)
        self.assertFalse(node.can_ignore_errors())

    def test_template_skips_checks(self):
        template_analysis = Analysis(genome_build=self.genome_build, template_type=AnalysisTemplateType.TEMPLATE)
        template_analysis.set_defaults_and_save(self.user)
        node = TrioNode.objects.create(analysis=template_analysis, trio=self.trio,
                                       inheritance=TrioInheritance.XLINKED_RECESSIVE)
        self.assertEqual({}, node.get_field_errors())
        self.assertEqual([], node.get_warnings())

    def _form(self, ignore):
        node = self._make_node()
        data = {"trio": self.trio.pk, "inheritance": TrioInheritance.XLINKED_RECESSIVE}
        if ignore:
            data["ignore_field_errors"] = "on"
        form = TrioNodeForm(data=data, instance=node)
        form.is_valid()
        return form

    def test_form_reports_inheritance_error(self):
        form = self._form(ignore=False)
        self.assertIn("inheritance", form.errors)
        self.assertFalse(form.instance.ignore_field_errors)

    def test_form_ignore_lets_inheritance_through(self):
        form = self._form(ignore=True)
        self.assertNotIn("inheritance", form.errors)
        self.assertTrue(form.instance.ignore_field_errors)


class TestRevealHiddenNodesIgnoresErrors(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_reveal_ignore')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        cls.trio = create_fake_trio(cls.user, cls.genome_build)
        cls.analysis = Analysis(genome_build=cls.genome_build)
        cls.analysis.set_defaults_and_save(cls.user)

    def _reveal(self, node):
        self.client.force_login(self.user)
        url = reverse('node_reveal_hidden', kwargs={"analysis_id": self.analysis.pk, "node_id": node.pk})
        self.assertEqual(200, self.client.post(url).status_code)
        node.refresh_from_db()
        return node

    def test_ignorable_errors_waived(self):
        node = TrioNode.objects.create(analysis=self.analysis, visible=False, trio=self.trio,
                                       inheritance=TrioInheritance.XLINKED_RECESSIVE)
        node = self._reveal(node)
        self.assertTrue(node.visible)
        self.assertTrue(node.ignore_field_errors)
        self.assertTrue(node.valid)

    def test_other_errors_left_alone(self):
        node = TrioNode.objects.create(analysis=self.analysis, visible=False, trio=None)
        node = self._reveal(node)
        self.assertTrue(node.visible)
        self.assertFalse(node.ignore_field_errors)
