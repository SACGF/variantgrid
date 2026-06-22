"""
Tests for VAV-archive guards on annotation queryset builders.

@see claude/issue_1536_data_archive_plan.md §3
"""

from django.contrib.auth.models import User
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from analysis.models import Analysis
from annotation.annotation_version_querysets import (
    get_variant_queryset_for_annotation_version,
    get_variants_qs_for_annotation,
)
from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import AnnotationVersion, VariantAnnotationVersion, VariantTranscriptAnnotation
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from snpdb.archive import DataArchivedError
from snpdb.models import GenomeBuild


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationArchiveGuardTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        cls.vav = VariantAnnotationVersion.objects.create(
            **kwargs, status=VariantAnnotationVersion.Status.ACTIVE
        )
        cls.av = AnnotationVersion.latest(cls.grch37, validate=False)

    def _archive_vav(self):
        self.vav.data_archived_date = timezone.now()
        self.vav.save()

    def test_get_variant_queryset_raises_when_archived(self):
        self._archive_vav()
        av = AnnotationVersion.objects.select_related("variant_annotation_version").get(pk=self.av.pk)
        with self.assertRaises(DataArchivedError):
            get_variant_queryset_for_annotation_version(av)

    def test_get_variants_qs_for_annotation_raises_when_archived(self):
        self._archive_vav()
        av = AnnotationVersion.objects.select_related("variant_annotation_version").get(pk=self.av.pk)
        with self.assertRaises(DataArchivedError):
            get_variants_qs_for_annotation(av)

    def test_overlapping_genes_q_raises_when_archived(self):
        self._archive_vav()
        gene_qs = Gene.objects.none()
        with self.assertRaises(DataArchivedError):
            VariantTranscriptAnnotation.get_overlapping_genes_q(self.vav, gene_qs)

    def test_analysis_warns_and_errors_on_archived_vav(self):
        """ Analysis.get_warnings + get_errors surface VAV-archive. """
        user = User.objects.create_user(username="anal_user")
        analysis = Analysis.objects.create(
            user=user, genome_build=self.grch37, name="t", annotation_version=self.av
        )
        self.assertNotIn("archived", " ".join(analysis.get_warnings()).lower())
        self._archive_vav()
        analysis = Analysis.objects.select_related("annotation_version__variant_annotation_version").get(pk=analysis.pk)
        warnings = analysis.get_warnings()
        errors = analysis.get_errors()
        self.assertTrue(any("archived" in w.lower() for w in warnings),
                        f"warnings did not mention archive: {warnings}")
        self.assertTrue(any("archived" in e.lower() for e in errors),
                        f"errors did not mention archive: {errors}")
