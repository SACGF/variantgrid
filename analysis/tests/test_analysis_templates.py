"""
#1496 - saving a template version creates a draft only its writers can run. It goes live when
someone explicitly makes it active, which is also the rollback path for an older version.
"""
import json

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls.base import reverse
from django.utils import timezone
from guardian.shortcuts import assign_perm

from analysis.analysis_templates import get_sample_analysis
from analysis.models import (
    Analysis,
    AnalysisTemplate,
    AnalysisTemplateRun,
    AnalysisTemplateType,
    AnalysisTemplateVersion,
    AnalysisVariable,
)
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.unittest_utils import prevent_request_warnings
from library.guardian_utils import DjangoPermission, assign_permission_to_user_and_groups
from snpdb.models import VCF, GenomeBuild, ImportStatus, Sample

ANALYSIS_NAME_TEMPLATE = "%(template)s for %(input)s"


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class AnalysisTemplateDraftTestCase(TestCase):
    """ A template with a single sample source node, so it can be launched against self.sample """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.get_or_create(username="test_analysis_templates_owner")[0]
        cls.non_owner = User.objects.get_or_create(username="test_analysis_templates_non_owner")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.vcf = VCF.objects.create(name="template_vcf", genotype_samples=1, genome_build=cls.grch37,
                                     import_status=ImportStatus.SUCCESS, user=cls.owner, date=timezone.now())
        cls.sample = Sample.objects.create(name="template_sample", vcf=cls.vcf,
                                           import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.owner, cls.vcf)
        assign_permission_to_user_and_groups(cls.owner, cls.sample)
        cls._grant_read(cls.non_owner, cls.sample)

        cls.template = cls._create_template("sample template")

    @staticmethod
    def _grant_read(user, obj):
        assign_perm(DjangoPermission.perm(obj, DjangoPermission.READ), user, obj)

    @classmethod
    def _create_template(cls, name) -> AnalysisTemplate:
        analysis = Analysis(genome_build=cls.grch37, template_type=AnalysisTemplateType.TEMPLATE)
        analysis.set_defaults_and_save(cls.owner)
        node = SampleNode.objects.create(analysis=analysis)
        AnalysisVariable.objects.create(node=node, field="sample", class_name="snpdb.Sample")
        return AnalysisTemplate.objects.create(name=name, user=cls.owner, analysis=analysis)

    def _new_version(self) -> AnalysisTemplateVersion:
        # Fresh instance - Analysis.clone() re-points the analysis object it is handed at the snapshot
        template = AnalysisTemplate.objects.get(pk=self.template.pk)
        return template.new_version(ANALYSIS_NAME_TEMPLATE)

    def _grant_template_read(self, user):
        self._grant_read(user, AnalysisTemplate.objects.get(pk=self.template.pk).analysis)


class TestVersionLifecycle(AnalysisTemplateDraftTestCase):

    def test_new_version_creates_a_draft(self):
        first = self._new_version()
        self.assertFalse(first.active)
        self.assertIsNone(self.template.active)
        self.assertEqual(self.template.draft, first)

        second = self._new_version()
        self.assertIsNone(self.template.active)
        self.assertEqual(self.template.draft, second)
        # Only the latest is a draft - the earlier one is now history
        first.refresh_from_db()
        self.assertFalse(first.is_draft)

    def test_activating_an_older_version_leaves_the_latest_as_the_draft(self):
        first = self._new_version()
        second = self._new_version()

        first.activate()
        self.assertEqual(self.template.active, first)
        second.refresh_from_db()
        self.assertTrue(second.is_draft)

        second.activate()
        first.refresh_from_db()
        self.assertFalse(first.active)
        self.assertFalse(first.is_draft)
        self.assertEqual(self.template.active, second)

    def test_draft_is_named_as_one(self):
        draft = self._new_version()
        self.assertIn("(draft)", str(draft))
        draft.activate()
        self.assertNotIn("(draft)", str(draft))


class TestFilterForUser(AnalysisTemplateDraftTestCase):

    def setUp(self):
        super().setUp()
        self._grant_template_read(self.non_owner)

    def test_a_template_with_only_a_draft_is_invisible_to_others(self):
        draft = self._new_version()
        self.assertEqual(list(AnalysisTemplateVersion.filter_for_user(self.owner)), [draft])
        self.assertEqual(list(AnalysisTemplateVersion.filter_for_user(self.non_owner)), [])

    def test_others_only_see_the_active_version(self):
        active = self._new_version()
        active.activate()
        draft = self._new_version()

        self.assertEqual(set(AnalysisTemplateVersion.filter_for_user(self.owner)), {active, draft})
        self.assertEqual(list(AnalysisTemplateVersion.filter_for_user(self.non_owner)), [active])

    def test_history_is_offered_to_nobody(self):
        old = self._new_version()
        old.activate()
        active = self._new_version()
        active.activate()
        draft = self._new_version()

        self.assertEqual(set(AnalysisTemplateVersion.filter_for_user(self.owner)), {active, draft})


class TestTemplateRunVersion(AnalysisTemplateDraftTestCase):

    def test_no_active_version_is_an_error(self):
        self._new_version()  # draft only
        with self.assertRaises(ValueError):
            AnalysisTemplateRun.create(self.template, self.grch37, user=self.owner)

    def test_run_uses_the_version_it_was_given(self):
        active = self._new_version()
        active.activate()
        draft = self._new_version()

        template_run = AnalysisTemplateRun.create(self.template, self.grch37, user=self.owner,
                                                  template_version=draft)
        self.assertEqual(template_run.template_version, draft)
        self.assertEqual(template_run.analysis.analysistemplaterun.template_version, draft)

        template_run.populate_arguments({"sample": self.sample})
        template_run.populate_analysis_name()
        self.assertIn("(draft)", template_run.analysis.name)

    def test_get_sample_analysis_reuses_its_run_when_a_draft_sits_on_top(self):
        self._new_version().activate()
        self._new_version()  # draft above the active version

        analysis = get_sample_analysis(self.sample, self.template)
        self.assertEqual(get_sample_analysis(self.sample, self.template), analysis)


class TestLaunchViews(AnalysisTemplateDraftTestCase):

    def _create_from_template(self, user, atv):
        self.client.force_login(user)
        url = reverse("create_analysis_from_template", kwargs={"genome_build_name": self.grch37.name})
        return self.client.post(url, {"tag_uuid": "tag",
                                      "tag-analysis_template_version": atv.pk,
                                      "sample": self.sample.pk})

    def test_anyone_can_run_the_active_version(self):
        active = self._new_version()
        active.activate()
        self._grant_template_read(self.non_owner)

        for user in (self.owner, self.non_owner):
            with self.subTest(user=user.username):
                response = self._create_from_template(user, active)
                self.assertEqual(response.status_code, 302)

    def test_the_owner_can_run_a_draft(self):
        response = self._create_from_template(self.owner, self._new_version())
        self.assertEqual(response.status_code, 302)

    @prevent_request_warnings
    def test_a_viewer_cannot_run_a_draft(self):
        draft = self._new_version()
        self._grant_template_read(self.non_owner)
        self.assertEqual(self._create_from_template(self.non_owner, draft).status_code, 403)


class TestActivateView(AnalysisTemplateDraftTestCase):

    def _activate(self, user, atv):
        self.client.force_login(user)
        return self.client.post(reverse("analysis_template_version_activate", kwargs={"pk": atv.pk}))

    def test_writer_activates(self):
        draft = self._new_version()
        response = self._activate(self.owner, draft)
        self.assertEqual(response.status_code, 302)
        draft.refresh_from_db()
        self.assertTrue(draft.active)

    @prevent_request_warnings
    def test_viewer_cannot_activate(self):
        draft = self._new_version()
        self._grant_template_read(self.non_owner)
        self.assertEqual(self._activate(self.non_owner, draft).status_code, 403)
        draft.refresh_from_db()
        self.assertFalse(draft.active)


class TestSettingsPage(AnalysisTemplateDraftTestCase):

    def test_active_and_draft_are_both_editable(self):
        active = self._new_version()
        active.activate()
        draft = self._new_version()

        self.client.force_login(self.owner)
        url = reverse("analysis_template_settings", kwargs={"pk": self.template.pk})
        content = self.client.get(url).content.decode()

        # A settings form each, so a draft can be configured before it goes live
        for atv in (active, draft):
            self.assertIn(f'name="atv-{atv.pk}"', content)
        # ... and a make-active form for the draft alone
        self.assertIn(reverse("analysis_template_version_activate", kwargs={"pk": draft.pk}), content)
        self.assertNotIn(reverse("analysis_template_version_activate", kwargs={"pk": active.pk}), content)

    def test_saving_a_draft_leaves_the_active_versions_settings_alone(self):
        active = self._new_version()
        active.activate()
        draft = self._new_version()

        self.client.force_login(self.owner)
        url = reverse("analysis_template_settings", kwargs={"pk": self.template.pk})
        response = self.client.post(url, {f"atv-{draft.pk}": "save",
                                          f"atv-{draft.pk}-analysis_name_template": "draft name",
                                          f"atv-{draft.pk}-appears_in_links": "on"})
        self.assertEqual(response.status_code, 200)
        draft.refresh_from_db()
        active.refresh_from_db()
        self.assertEqual(draft.analysis_name_template, "draft name")
        self.assertTrue(draft.appears_in_links)
        self.assertEqual(active.analysis_name_template, ANALYSIS_NAME_TEMPLATE)


class TestTemplateSaveView(AnalysisTemplateDraftTestCase):

    def _save(self, **data):
        self.client.force_login(self.owner)
        url = reverse("analysis_template_save", kwargs={"pk": self.template.pk})
        response = self.client.post(url, {"analysis_name_template": ANALYSIS_NAME_TEMPLATE, **data})
        self.assertEqual(response.status_code, 200)
        return json.loads(response.content)

    def test_save_creates_a_draft_and_leaves_the_active_version_alone(self):
        active = self._new_version()
        active.activate()

        data = self._save()
        self.assertFalse(data["active"])
        atv = self.template.latest_version_obj
        self.assertEqual(atv.version, data["version"])
        self.assertEqual(atv.analysis_name_template, ANALYSIS_NAME_TEMPLATE)
        active.refresh_from_db()
        self.assertTrue(active.active)

    def test_save_and_make_active_replaces_the_previous_version(self):
        previous = self._new_version()
        previous.activate()

        data = self._save(activate=1)
        self.assertTrue(data["active"])
        self.assertEqual(data["replaced_version"], previous.version)
        previous.refresh_from_db()
        self.assertFalse(previous.active)
        self.assertEqual(self.template.active.version, data["version"])

    def test_first_activated_version_replaces_nothing(self):
        data = self._save(activate=1)
        self.assertTrue(data["active"])
        self.assertIsNone(data["replaced_version"])
