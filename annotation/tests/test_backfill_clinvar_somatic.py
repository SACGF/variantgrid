from datetime import datetime

from django.core.management import call_command
from django.test import TestCase
from django.utils.timezone import make_aware

from annotation.fake_annotation import create_fake_variants
from annotation.models.models import ClinVar, ClinVarVersion
from annotation.models.models_enums import ClinVarOncogenicity
from classification.enums import SomaticClinicalSignificance
from snpdb.models import GenomeBuild, Variant


def _utc(year: int) -> datetime:
    return make_aware(datetime(year, 1, 1))


class BackfillClinVarSomaticTest(TestCase):
    """ The derived somatic_tier / highest_oncogenicity columns for ClinVar versions loaded before
        they were written by the import """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_variants(cls.genome_build)
        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        cls.somatic_variant, cls.oncogenic_variant, cls.germline_variant = variants

        cls.old_version = cls._create_version("clinvar_20240101.vcf.gz", _utc(2024))
        cls.latest_version = cls._create_version("clinvar_20260101.vcf.gz", _utc(2026))
        for clinvar_version in (cls.old_version, cls.latest_version):
            cls._create_clinvar_rows(clinvar_version)

    @classmethod
    def _create_version(cls, filename: str, annotation_date) -> ClinVarVersion:
        clinvar_version = ClinVarVersion.objects.create(filename=filename, sha256_hash=filename,
                                                        genome_build=cls.genome_build)
        # annotation_date is auto_now_add, so the version order has to be set after the insert
        ClinVarVersion.objects.filter(pk=clinvar_version.pk).update(annotation_date=annotation_date)
        return clinvar_version

    @classmethod
    def _create_clinvar_rows(cls, clinvar_version: ClinVarVersion):
        common = {"version": clinvar_version, "clinvar_variation_id": 1, "clinvar_allele_id": 1}
        ClinVar.objects.create(variant=cls.somatic_variant,
                               somatic_clinical_significance="Tier_I_-_Strong", **common)
        ClinVar.objects.create(variant=cls.oncogenic_variant,
                               oncogenic_classification="Likely_oncogenic", **common)
        ClinVar.objects.create(variant=cls.germline_variant,
                               clinical_significance="Pathogenic", highest_pathogenicity=5, **common)

    def _get(self, clinvar_version: ClinVarVersion, variant: Variant) -> ClinVar:
        return ClinVar.objects.get(version=clinvar_version, variant=variant)

    def test_latest_version_gains_derived_values(self):
        call_command("one_off_backfill_clinvar_somatic")
        somatic = self._get(self.latest_version, self.somatic_variant)
        self.assertEqual(somatic.somatic_tier, SomaticClinicalSignificance.TIER_1)
        oncogenic = self._get(self.latest_version, self.oncogenic_variant)
        self.assertEqual(oncogenic.highest_oncogenicity, ClinVarOncogenicity.LIKELY_ONCOGENIC)

    def test_germline_only_row_untouched(self):
        call_command("one_off_backfill_clinvar_somatic")
        germline = self._get(self.latest_version, self.germline_variant)
        self.assertIsNone(germline.somatic_tier)
        self.assertEqual(germline.highest_oncogenicity, 0)
        self.assertEqual(germline.highest_pathogenicity, 5)

    def test_older_version_left_alone(self):
        call_command("one_off_backfill_clinvar_somatic")
        self.assertIsNone(self._get(self.old_version, self.somatic_variant).somatic_tier)

    def test_older_version_when_named(self):
        call_command("one_off_backfill_clinvar_somatic", clinvar_version=self.old_version.pk)
        self.assertEqual(self._get(self.old_version, self.somatic_variant).somatic_tier,
                         SomaticClinicalSignificance.TIER_1)
