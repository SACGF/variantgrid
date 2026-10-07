import os
import tempfile

import cyvcf2
from django.test import TestCase

from annotation.annotation_pipeline_routing import (
    pipeline_type_for_variant,
    pipeline_type_variant_q,
)
from annotation.annotation_run_files import write_qs_to_vcf
from annotation.fake_data import get_fake_annotation_version
from annotation.models.models_enums import VariantAnnotationPipelineType
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import GenomeBuild, Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)

STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


class AnnotationPipelineRoutingTest(TestCase):
    """ #1358 - a symbolic del/dup/inv under ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE is annotated as a small variant """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        def _symbolic(alt, svlen):
            vc = VariantCoordinate(chrom="3", position=128200125, ref="G", alt=alt, svlen=svlen)
            return slowly_create_test_variant_from_coordinate(vc, cls.grch37)

        cls.short_del = _symbolic(VCFSymbolicAllele.DEL, -60)
        cls.expected = {
            slowly_create_test_variant("3", 128198980, "A", "T", cls.grch37): STANDARD,
            cls.short_del: STANDARD,
            _symbolic(VCFSymbolicAllele.DUP, 999): STANDARD,
            _symbolic(VCFSymbolicAllele.DEL, -1000): STRUCTURAL,
            _symbolic(VCFSymbolicAllele.CNV, 60): STRUCTURAL,  # No explicit form, whatever the size
        }
        cls.long_del = _symbolic(VCFSymbolicAllele.DEL, -2000)

    def test_variant_and_queryset_routing_agree(self):
        qs = Variant.objects.filter(pk__in=[v.pk for v in self.expected])
        for pipeline_type in (STANDARD, STRUCTURAL):
            expected_ids = {v.pk for v, pt in self.expected.items() if pt == pipeline_type}
            self.assertEqual(expected_ids, set(qs.filter(pipeline_type_variant_q(pipeline_type))
                                               .values_list("pk", flat=True)), pipeline_type)
        for variant, pipeline_type in self.expected.items():
            self.assertEqual(pipeline_type, pipeline_type_for_variant(variant), variant)

    def test_dump_writes_short_symbolic_explicit(self):
        qs = Variant.objects.filter(pk__in=[self.short_del.pk, self.long_del.pk])
        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf_filename = os.path.join(tmp_dir, "dump.vcf")
            write_qs_to_vcf(vcf_filename, self.grch37, qs, small_symbolic_as_explicit=True)
            records = {r.INFO["variant_id"]: r for r in cyvcf2.VCF(vcf_filename)}

        short = records[self.short_del.pk]
        self.assertEqual((61, 1), (len(short.REF), len(short.ALT[0])))
        self.assertEqual([VCFSymbolicAllele.DEL], records[self.long_del.pk].ALT)
