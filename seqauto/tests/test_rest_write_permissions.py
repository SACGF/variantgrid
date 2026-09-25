from django.conf import settings
from django.contrib.auth.models import Group, User
from django.urls import reverse
from rest_framework.test import APITestCase

from seqauto.models import EnrichmentKit, Experiment


class SeqAutoRESTWritePermissionTest(APITestCase):
    """ Any logged-in user reads the seqauto API, only superusers and SEQAUTO_API_WRITE_GROUP write it """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="seqauto_api_reader")
        cls.pipeline_user = User.objects.create_user(username="seqauto_api_pipeline")
        write_group = Group.objects.get_or_create(name=settings.SEQAUTO_API_WRITE_GROUP)[0]
        cls.pipeline_user.groups.add(write_group)
        cls.enrichment_kit = EnrichmentKit.objects.create(name="seqauto_api_kit")

    def _create_experiment(self, user, name):
        self.client.force_authenticate(user=user)
        return self.client.post(reverse("api_experiment-list"), {"name": name}, format="json")

    def test_user_outside_write_group_can_read_but_not_write(self):
        response = self._create_experiment(self.user, "exp_reader")
        self.assertEqual(response.status_code, 403)
        self.assertFalse(Experiment.objects.filter(name="exp_reader").exists())

        url = reverse("api_enrichment_kit-detail", kwargs={"pk": self.enrichment_kit.pk})
        self.assertEqual(self.client.get(url).status_code, 200)
        self.assertEqual(self.client.delete(url).status_code, 403)
        self.assertTrue(EnrichmentKit.objects.filter(pk=self.enrichment_kit.pk).exists())

        response = self.client.post(reverse("api_sequencing_files_bulk_create"), {}, format="json")
        self.assertEqual(response.status_code, 403)

    def test_write_group_member_can_write(self):
        response = self._create_experiment(self.pipeline_user, "exp_pipeline")
        self.assertEqual(response.status_code, 201)
        self.assertTrue(Experiment.objects.filter(name="exp_pipeline").exists())
