import os

from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.fake_annotation import get_fake_annotation_settings_dict
from annotation.models import VariantAnnotation, VariantAnnotationPipelineType
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import AnnotationRun, VariantAnnotationVersion
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import get_vep_version_from_vcf, vep_dict_to_variant_annotation_version_kwargs, \
    VEPConfig
from library.genomics.vcf_enums import VariantClass
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf


# CNV stuff only comes in with columns version >= 3
@override_settings(**get_fake_annotation_settings_dict(columns_version=3))
class TestAnnotationVCFCNV(TestCase):
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
    TEST_ANNOTATION_VCF_GRCH37 = os.path.join(TEST_DATA_DIR, "test_columns_version3_grch37_sv.vep_annotated.vcf")
    TEST_ANNOTATION_VCF_GRCH38 = os.path.join(TEST_DATA_DIR, "test_columns_version3_grch38_sv.vep_annotated.vcf")

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.variant_annotation_versions_by_build = {}

        VCFS = {
            "GRCh37": cls.TEST_ANNOTATION_VCF_GRCH37,
            "GRCh38": cls.TEST_ANNOTATION_VCF_GRCH38,
        }

        for genome_build_name, vcf in VCFS.items():
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
            vep_config = VEPConfig(genome_build)
            vep_dict = get_vep_version_from_vcf(vcf)
            kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_dict)
            kwargs["genome_build"] = genome_build
            vav, created = VariantAnnotationVersion.objects.get_or_create(**kwargs)
            if not created:
                print("Truncating!")
                vav.truncate_related_objects()
            cls.variant_annotation_versions_by_build[genome_build_name] = vav

            slowly_create_loci_and_variants_for_vcf(genome_build, vcf, get_variant_id_from_info=True)

    def test_import_variant_annotations_grch37(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh37')
        vav = self.variant_annotation_versions_by_build[genome_build.name]

        print("Variants: ", Variant.objects.count())
        print("VariantAnnotation: ", VariantAnnotation.objects.count())
        print(f"VariantAnnotationVersion: {vav}")

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH37)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

        # Verify a few fields - got these from vep_reader --pick
        # 3	127535894	.	G	<INV>	.	.	SVTYPE=INV;SVLEN=1184482;variant_id=202
        va = VariantAnnotation.objects.get(variant_id=202)
        self.assertEqual(va.variant_class, VariantClass.INVERSION)
        self.assertEqual(va.impact, PathogenicityImpact.MODIFIER)

        # 17	41236500	.	G	<DEL>	.	.	SVTYPE=DEL;SVLEN=14500;variant_id=203
        va = VariantAnnotation.objects.get(variant_id=203)
        self.assertEqual(va.variant_class, VariantClass.DELETION)
        self.assertEqual(va.impact, PathogenicityImpact.HIGH)

        # 21	36418200	.	G	<DUP>	.	.	SVTYPE=DUP;SVLEN=6000;variant_id=205
        va = VariantAnnotation.objects.get(variant_id=205)
        self.assertEqual(va.variant_class, VariantClass.DUPLICATION)
        self.assertEqual(va.impact, PathogenicityImpact.MODIFIER)

    def test_import_variant_annotations_grch38(self):
        genome_build = GenomeBuild.get_name_or_alias('GRCh38')
        vav = self.variant_annotation_versions_by_build[genome_build.name]

        print("Variants: ", Variant.objects.count())
        print("VariantAnnotation: ", VariantAnnotation.objects.count())
        print(f"VariantAnnotationVersion: {vav}")

        annotation_range_lock, _ = get_annotation_range_lock_and_unannotated_count(vav)
        annotation_range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                      vcf_annotated_filename=self.TEST_ANNOTATION_VCF_GRCH38)
        import_vcf_annotations(annotation_run, delete_temp_files=False, vep_version_check=False)

        # Verify a few fields
        # 3	128486000	.	C	<DUP>	.	.	SVLEN=9951;SVTYPE=DUP;variant_id=101
        va = VariantAnnotation.objects.get(variant_id=101)
        self.assertEqual(va.variant_class, VariantClass.DUPLICATION)
        self.assertEqual(va.impact, PathogenicityImpact.MODIFIER)

        # 21	35041808	.	G	<DEL>	.	.	SVLEN=10000;SVTYPE=DEL;variant_id=102
        va = VariantAnnotation.objects.get(variant_id=102)
        self.assertEqual(va.variant_class, VariantClass.DELETION)
        self.assertEqual(va.impact, PathogenicityImpact.HIGH)
