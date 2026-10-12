from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from guardian.shortcuts import assign_perm

from snpdb.models.models_genomic_interval import (
    GenomicIntervalsCategory,
    GenomicIntervalsCollection,
)


class TestViewGenomicIntervals(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        owner = User.objects.create_user(username="gic_owner", password="x")
        cls.viewer = User.objects.create_user(username="gic_viewer", password="x")
        category = GenomicIntervalsCategory.objects.get_or_create(name="test category")[0]
        cls.gic = GenomicIntervalsCollection.objects.create(name="original", category=category, user=owner)
        assign_perm(GenomicIntervalsCollection.get_read_perm(), cls.viewer, cls.gic)

    def test_view_permission_cannot_save(self):
        self.client.force_login(self.viewer)
        url = reverse("view_genomic_intervals", kwargs={"genomic_intervals_collection_id": self.gic.pk})
        self.assertEqual(200, self.client.get(url).status_code)

        response = self.client.post(url, {"name": "renamed", "user": self.viewer.pk})
        self.assertEqual(403, response.status_code)
        self.gic.refresh_from_db()
        self.assertEqual("original", self.gic.name)
