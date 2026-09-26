"""Putting the space back in gene-level classification targets, and the 'splice' suffix onto
splice annotation and evidence - @see gene_level_strings_respace."""
from django.core.management import call_command
from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from annotation.gene_level_annotation import annotate_gene_level_run
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    VariantAnnotation,
    VariantAnnotationPipelineType,
)
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import Classification, ImportedAlleleInfo
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.models import HGNC, GeneSymbol, HGNCImport
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import create_splice_event_variant
from library.utils import md5sum_str
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion

AR_HGNC_ID = 644
FGF4_HGNC_ID = 3682


class GeneLevelStringsRespaceTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch38()
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(AR_HGNC_ID, "AR"), (FGF4_HGNC_ID, "FGF4")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _classification(self, data: dict) -> Classification:
        classification = Classification.create(user=self.user, lab=self.lab, data=data, save=True,
                                               source=SubmissionSource.API)
        classification.publish_latest(user=self.user)
        # What the import task does once the record exists
        classification.ensure_allele_info()
        classification.save()
        return classification

    def test_a_stripped_copy_number_target_gets_its_space_back(self):
        """ The imported value, its md5, and the evidence history all move together """
        classification = self._classification({SpecialEKeys.C_HGVS: {"value": "FGF4amplification"},
                                               SpecialEKeys.GENOME_BUILD: {"value": "GRCh38"}})
        allele_info = classification.allele_info
        self.assertEqual("FGF4amplification", allele_info.imported_c_hgvs)

        call_command("gene_level_strings_respace")

        allele_info.refresh_from_db()
        self.assertEqual("FGF4 amplification", allele_info.imported_c_hgvs)
        self.assertEqual(md5sum_str("FGF4 amplification"), allele_info.imported_md5_hash)
        classification.refresh_from_db()
        self.assertEqual("FGF4 amplification", classification.get(SpecialEKeys.C_HGVS))
        published = classification.last_published_version
        self.assertEqual("FGF4 amplification", published.published_evidence[SpecialEKeys.C_HGVS]["value"])
        self.assertEqual("FGF4 amplification", published.delta[SpecialEKeys.C_HGVS]["value"])

    def test_a_splice_target_is_left_as_the_lab_wrote_it(self):
        """ 'ARV7' and 'EGFRvIII' are both real spellings, so nothing says which one had a space """
        allele_info = ImportedAlleleInfo.objects.create(
            imported_c_hgvs="ARV7", imported_md5_hash=md5sum_str("ARV7"),
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                self.genome_build))
        call_command("gene_level_strings_respace")
        allele_info.refresh_from_db()
        self.assertEqual("ARV7", allele_info.imported_c_hgvs)

    def test_a_record_already_holding_the_spaced_value_is_reported_rather_than_merged(self):
        patch_version = GenomeBuildPatchVersion.get_unspecified_patch_version_for(self.genome_build)
        ImportedAlleleInfo.objects.create(imported_c_hgvs="FGF4 amplification",
                                          imported_md5_hash=md5sum_str("FGF4 amplification"),
                                          imported_genome_build_patch_version=patch_version)
        stripped = ImportedAlleleInfo.objects.create(imported_c_hgvs="FGF4amplification",
                                                     imported_md5_hash=md5sum_str("FGF4amplification"),
                                                     imported_genome_build_patch_version=patch_version)
        call_command("gene_level_strings_respace")
        stripped.refresh_from_db()
        self.assertEqual("FGF4amplification", stripped.imported_c_hgvs, "left for the user to merge")

    def test_splice_annotation_and_evidence_gain_the_suffix(self):
        vav = get_fake_annotation_version(self.genome_build).variant_annotation_version
        splice_variant = create_splice_event_variant("AR", "V7").variant
        range_lock = AnnotationRangeLock.objects.create(version=vav, min_variant=splice_variant,
                                                        max_variant=splice_variant, count=1)
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL)
        annotate_gene_level_run(annotation_run)
        # As the rows read before the suffix
        VariantAnnotation.objects.filter(version=vav, variant=splice_variant).update(hgvs_c="AR-V7", hgvs_g="AR-V7")
        classification = self._classification({SpecialEKeys.C_HGVS: {"value": "AR V7"},
                                               SpecialEKeys.G_HGVS: {"value": "AR-V7"},
                                               SpecialEKeys.SPLICE_LABEL: {"value": "AR-V7"},
                                               SpecialEKeys.GENOME_BUILD: {"value": "GRCh38"}})
        classification.allele_info.set_variant_and_save(matched_variant=splice_variant)

        call_command("gene_level_strings_respace")

        variant_annotation = VariantAnnotation.objects.get(version=vav, variant=splice_variant)
        self.assertEqual("AR-V7 splice", variant_annotation.hgvs_c)
        self.assertEqual("AR-V7 splice", variant_annotation.hgvs_g)
        classification.refresh_from_db()
        self.assertEqual("AR-V7 splice", classification.get(SpecialEKeys.G_HGVS))
        self.assertEqual("AR-V7 splice", classification.get(SpecialEKeys.SPLICE_LABEL))
        self.assertEqual("AR V7", classification.get(SpecialEKeys.C_HGVS), "the lab's own value")
