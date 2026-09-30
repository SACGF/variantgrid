from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import VariantTag
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.fake_data import (
    create_fake_variants,
    get_fake_annotation_version,
    retire_seeded_annotation_version,
)
from annotation.models import AnnotationRun, VariantTranscriptAnnotation
from classification.enums import AlleleOriginBucket, ShareLevel
from classification.models import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.fake_data import create_fake_transcript_version
from snpdb.models import GenomeBuild, Tag, Variant, VariantZygosityCountCollection
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele
from variantopedia.interesting_nearby import (
    get_transcript_and_domains,
    get_transcript_and_exons,
    get_transcripts_and_codons,
    interesting_summary,
)


class InterestingNearbyTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.grch37()
        cls.old_vav = get_fake_annotation_version(cls.genome_build).variant_annotation_version
        retire_seeded_annotation_version(cls.genome_build)
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.vav = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)
        cls.variant = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk").first()
        cls.allele = create_mock_allele(cls.variant, cls.genome_build)

        ClassificationTestUtils.setUp()
        cls.lab, cls.user = ClassificationTestUtils.lab_and_user()
        cls.other_user = User.objects.get_or_create(username="interesting_nearby_other_user")[0]

        for tag_id in ("Artefact", "Follow-up"):
            tag = Tag.objects.create(pk=tag_id)
            VariantTag.objects.create(variant=cls.variant, allele=cls.allele, tag=tag, analysis=None,
                                      genome_build=cls.genome_build, user=cls.user)
        for i in range(2):
            Classification.objects.create(lab=cls.lab, user=cls.user, lab_record_id=f"nearby_{i}",
                                          variant=cls.variant, allele=cls.allele,
                                          allele_origin_bucket=AlleleOriginBucket.GERMLINE,
                                          share_level=ShareLevel.ALL_USERS)

    def _variant_qs(self):
        qs = get_variant_queryset_for_annotation_version(self.annotation_version)
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        return qs.filter(pk=self.variant.pk)

    def test_tag_counts_are_per_tagging_not_per_classification(self):
        summary, tag_counts = interesting_summary(self._variant_qs(), self.user, self.genome_build)
        self.assertEqual({"Artefact": 1, "Follow-up": 1}, tag_counts)
        self.assertIn("1 variants", summary)

    def test_tag_counts_only_include_tags_the_viewer_can_see(self):
        _, tag_counts = interesting_summary(self._variant_qs(), self.other_user, self.genome_build)
        self.assertEqual({}, tag_counts)

    def test_transcript_lookups_use_the_given_annotation_version(self):
        transcript = create_fake_transcript_version(self.genome_build).transcript
        annotation_run = AnnotationRun.objects.create()
        VariantTranscriptAnnotation.objects.create(version=self.old_vav, variant=self.variant,
                                                   annotation_run=annotation_run, transcript=transcript,
                                                   exon="3/10", hgvs_c="NM_1.1:c.100A>G",
                                                   interpro_domain="OldDomain")
        VariantTranscriptAnnotation.objects.create(version=self.vav, variant=self.variant,
                                                   annotation_run=annotation_run, transcript=transcript,
                                                   exon="5/10", hgvs_c="NM_1.1:c.200A>G",
                                                   interpro_domain="Domain1&Domain2")
        self.assertEqual({transcript.pk: ":c.200"}, get_transcripts_and_codons(self.variant, self.vav))
        self.assertEqual({transcript.pk: "5/10"}, get_transcript_and_exons(self.variant, self.vav))
        self.assertEqual({transcript.pk: {"Domain1", "Domain2"}},
                         get_transcript_and_domains(self.variant, self.vav))
