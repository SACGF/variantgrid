from django.contrib.auth.models import User
from django.test import TestCase
from django.utils.timezone import localdate

from analysis.models.models_variant_tag import VariantTagsImport
from annotation.models.models import ClinVarVersion, ManualVariantEntryCollection
from classification.models.classification import ClassificationImport
from snpdb.models import ImportedWikiCollection, ImportSource, VCF
from snpdb.models.models_enums import AlleleConversionTool
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import LiftoverRun
from upload.models import FileUpload, UploadedClassificationImport, UploadedClinVarVersion, UploadedFileTypes, \
    UploadedLiftover, UploadedManualVariantEntryCollection, UploadedVariantTags, UploadedVCF, \
    UploadedWikiCollection, UploadPipeline


class UploadPipelineGenomeBuildTest(TestCase):
    """ UploadPipeline.genome_build resolves via UploadData - one case per UploadedFileTypes value """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="upload_data_user", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _pipeline(self, file_type) -> UploadPipeline:
        file_upload = FileUpload.objects.create(user=self.user, name=f"{file_type}.file",
                                                path=f"/tmp/{file_type}.file", file_type=file_type,
                                                import_source=ImportSource.COMMAND_LINE)
        return UploadPipeline.objects.create(file_upload=file_upload)

    def _uploaded_vcf_without_vcf(self, pipeline: UploadPipeline):
        """ The VCF-ish pipelines create an UploadedVCF with a null vcf alongside their own satellite """
        UploadedVCF.objects.create(file_upload=pipeline.file_upload, upload_pipeline=pipeline)

    def test_vcf(self):
        pipeline = self._pipeline(UploadedFileTypes.VCF)
        vcf = VCF.objects.create(name="build.vcf", date=localdate(), user=self.user,
                                 genotype_samples=0, genome_build=self.grch37)
        UploadedVCF.objects.create(file_upload=pipeline.file_upload, upload_pipeline=pipeline, vcf=vcf)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_vcf_before_vcf_created(self):
        pipeline = self._pipeline(UploadedFileTypes.VCF)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertIsNone(pipeline.genome_build)

    def test_manual_variant_entry(self):
        pipeline = self._pipeline(UploadedFileTypes.MANUAL_VARIANT_ENTRY)
        mvec = ManualVariantEntryCollection.objects.create(user=self.user, genome_build=self.grch37)
        UploadedManualVariantEntryCollection.objects.create(file_upload=pipeline.file_upload, collection=mvec)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_vcf_insert_variants_only_manual_entry(self):
        pipeline = self._pipeline(UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY)
        mvec = ManualVariantEntryCollection.objects.create(user=self.user, genome_build=self.grch37)
        UploadedManualVariantEntryCollection.objects.create(file_upload=pipeline.file_upload, collection=mvec)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_vcf_insert_variants_only_classification_import(self):
        pipeline = self._pipeline(UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY)
        classification_import = ClassificationImport.objects.create(user=self.user, genome_build=self.grch37)
        UploadedClassificationImport.objects.create(file_upload=pipeline.file_upload,
                                                    classification_import=classification_import)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_liftover(self):
        pipeline = self._pipeline(UploadedFileTypes.LIFTOVER)
        liftover = LiftoverRun.objects.create(user=self.user, genome_build=self.grch37,
                                              conversion_tool=AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        UploadedLiftover.objects.create(file_upload=pipeline.file_upload, liftover=liftover)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_variant_tags(self):
        pipeline = self._pipeline(UploadedFileTypes.VARIANT_TAGS)
        variant_tags_import = VariantTagsImport.objects.create(user=self.user, genome_build=self.grch37)
        UploadedVariantTags.objects.create(file_upload=pipeline.file_upload,
                                           variant_tags_import=variant_tags_import)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_clinvar(self):
        pipeline = self._pipeline(UploadedFileTypes.CLINVAR)
        clinvar_version = ClinVarVersion.objects.create(filename="clinvar_20240101.vcf.gz", sha256_hash="abc",
                                                        genome_build=self.grch37)
        UploadedClinVarVersion.objects.create(file_upload=pipeline.file_upload, clinvar_version=clinvar_version)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_wiki_variant(self):
        pipeline = self._pipeline(UploadedFileTypes.WIKI_VARIANT)
        wiki_collection = ImportedWikiCollection.objects.create(match_column_name="Variant",
                                                                genome_build=self.grch37)
        UploadedWikiCollection.objects.create(file_upload=pipeline.file_upload, wiki_collection=wiki_collection)
        self._uploaded_vcf_without_vcf(pipeline)
        self.assertEqual(pipeline.genome_build, self.grch37)

    def test_wiki_gene_has_no_build(self):
        pipeline = self._pipeline(UploadedFileTypes.WIKI_GENE)
        wiki_collection = ImportedWikiCollection.objects.create(match_column_name="Gene Symbol")
        UploadedWikiCollection.objects.create(file_upload=pipeline.file_upload, wiki_collection=wiki_collection)
        self.assertIsNone(pipeline.genome_build)

    def test_build_less_file_types(self):
        for file_type in (UploadedFileTypes.BED, UploadedFileTypes.GENE_LIST, UploadedFileTypes.GENE_COVERAGE,
                          UploadedFileTypes.PED, UploadedFileTypes.PATIENT_RECORDS,
                          UploadedFileTypes.VARIANT_CLASSIFICATIONS):
            with self.subTest(file_type=file_type):
                self.assertIsNone(self._pipeline(file_type).genome_build)
