from django.test import TestCase, override_settings

from annotation.fake_data import create_fake_variant_annotation, get_fake_annotation_version
from classification.enums import ShareLevel, SpecialEKeys, SubmissionSource
from classification.models import Classification, ClassificationModification, ImportedAlleleInfo
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.enums import AnnotationConsortium
from genes.fake_data import create_fake_gene_version
from genes.signals.gene_symbol_search import gene_symbol_preview_classifications_extra
from snpdb.enums import AlleleConversionTool, AlleleOrigin
from snpdb.models import Allele, GenomeBuild, GenomeBuildPatchVersion, VariantAllele
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True)
class ClassificationGeneSymbolFilterTest(TestCase):
    """ A classification is in a gene by its allele's annotation, or by the gene symbol in its evidence """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.grch37()
        cls.vav = get_fake_annotation_version(cls.grch37).variant_annotation_version
        cls.gene_symbol = create_fake_gene_version(cls.grch37, "ENSG00000159216", "RUNX1",
                                                   AnnotationConsortium.ENSEMBL).gene_symbol

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _allele(self, position: int, gene=None) -> Allele:
        variant = slowly_create_test_variant("1", position, "A", "T", self.grch37)
        allele = Allele.objects.create()
        VariantAllele.objects.create(variant=variant, genome_build=self.grch37, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE,
                                     allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
        if gene:
            variant_annotation = create_fake_variant_annotation(variant, self.vav)
            variant_annotation.gene = gene
            variant_annotation.save()
        return allele

    def _classify(self, allele: Allele, gene_symbol: str, allele_origin: str = "germline") -> Classification:
        data = {SpecialEKeys.GENE_SYMBOL: {"value": gene_symbol}, SpecialEKeys.ALLELE_ORIGIN: {"value": allele_origin}}
        classification = Classification.create(user=self.user, lab=self.lab, data=data, save=True,
                                               source=SubmissionSource.API)
        classification.publish_latest(user=self.user, share_level=ShareLevel.ALL_USERS)
        allele_info = ImportedAlleleInfo.objects.create(
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(self.grch37),
            imported_c_hgvs=f"NM_001754.5:c.{classification.pk}A>T", allele=allele)
        Classification.objects.filter(pk=classification.pk).update(allele_info=allele_info)
        return classification

    def test_matches_on_allele_annotation_or_evidence_gene_symbol(self):
        gene = self.gene_symbol.geneversion_set.get().gene
        annotated_in_gene = self._classify(self._allele(1000, gene=gene), gene_symbol="NOT_RUNX1")
        named_in_evidence = self._classify(self._allele(1100), gene_symbol="runx1")
        self._classify(self._allele(1200), gene_symbol="BRCA1")

        q = classification_gene_symbol_filter(self.gene_symbol)
        matched = set(ClassificationModification.latest_for_user(user=self.user).filter(q)
                      .values_list("classification_id", flat=True))
        self.assertEqual(matched, {annotated_in_gene.pk, named_in_evidence.pk})

    def test_preview_counts_each_allele_origin(self):
        self._classify(self._allele(1000), gene_symbol="RUNX1")
        self._classify(self._allele(1100), gene_symbol="RUNX1")
        self._classify(self._allele(1200), gene_symbol="RUNX1", allele_origin="somatic")

        extras = gene_symbol_preview_classifications_extra(None, user=self.user, obj=self.gene_symbol)
        self.assertEqual({extra.key: extra.value for extra in extras},
                         {"Germline Classification Count": 2, "Somatic Classification Count": 1})
