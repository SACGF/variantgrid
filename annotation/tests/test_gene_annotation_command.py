from unittest.mock import patch

from django.test import TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import (
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.management.commands.gene_annotation import Command
from annotation.models import AnnotationVersion, GeneAnnotationVersion
from snpdb.models import GenomeBuild


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class CreateAndPopulateGeneAnnotationVersionTests(TestCase):
    """ Saving a GeneAnnotationVersion bumps a new AnnotationVersion pointing at it, so a populate step
        that dies used to leave a live version with no gene annotation - and gene level filters
        (eg the OMIM built in filter) then silently match nothing """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)

    def _make_gene_annotation_version(self, *args, **kwargs) -> GeneAnnotationVersion:
        """ Stands in for _create_gene_annotation_version, which needs GnomADGeneConstraint data """
        gav = self.annotation_version.gene_annotation_version
        return GeneAnnotationVersion.objects.create(
            gene_annotation_release=gav.gene_annotation_release,
            ontology_version=gav.ontology_version,
            gnomad_import_date=gav.gnomad_import_date)

    def _run(self, populate_side_effect=None):
        command = Command()
        with patch.object(Command, "_create_gene_annotation_version",
                          side_effect=self._make_gene_annotation_version), \
             patch.object(Command, "_populate_gene_annotation_version",
                          side_effect=populate_side_effect):
            return command._create_and_populate_gene_annotation_version(
                self.annotation_version.gene_annotation_version.gene_annotation_release,
                self.annotation_version.ontology_version, None, set())

    def test_failed_populate_leaves_annotation_versions_untouched(self):
        original = set(AnnotationVersion.objects.values_list("pk", "gene_annotation_version_id"))
        with self.assertRaises(ValueError):
            self._run(populate_side_effect=ValueError("populate blew up"))
        self.assertEqual(set(AnnotationVersion.objects.values_list("pk", "gene_annotation_version_id")),
                         original)

    def test_failed_populate_leaves_the_empty_version_unreferenced(self):
        with self.assertRaises(ValueError):
            self._run(populate_side_effect=ValueError("populate blew up"))
        referenced = set(AnnotationVersion.objects.values_list("gene_annotation_version_id", flat=True))
        orphan = GeneAnnotationVersion.objects.order_by("pk").last()
        self.assertNotIn(orphan.pk, referenced)

    def test_successful_populate_points_annotation_version_at_the_new_version(self):
        gav = self._run()
        latest = AnnotationVersion.latest(self.grch37, validate=False)
        self.assertEqual(latest.gene_annotation_version_id, gav.pk)
