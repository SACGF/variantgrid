import os
import tempfile

import cyvcf2
from django.test import TestCase

from annotation.annotation_pipeline_routing import (
    pipeline_type_for_variant,
    pipeline_type_variant_q,
)
from annotation.annotation_run_files import write_qs_to_vcf
from annotation.enums import VariantAnnotationPipelineType
from annotation.fake_data import get_fake_annotation_version
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils.hash_utils import sha256sum_str
from snpdb.models import GenomeBuild, Locus, Sequence, Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)

STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


class AnnotationPipelineRoutingTest(TestCase):
    """ #1358 - a symbolic del/dup/inv under structural_variant_min_size is annotated as a small variant """
    SV_MIN_SIZE = 1000

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)

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

        contig = cls.grch37.contigs.get(name="3")

        def _explicit(position, ref, alt):
            """ Stored as given - slowly_create_test_variant would make a long insertion a <DUP> """
            ref_seq = Sequence.objects.get_or_create(seq=ref, seq_sha256_hash=sha256sum_str(ref))[0]
            alt_seq = Sequence.objects.get_or_create(seq=alt, seq_sha256_hash=sha256sum_str(alt))[0]
            locus = Locus.objects.get_or_create(contig=contig, position=position, ref=ref_seq)[0]
            end = VariantCoordinate.from_explicit_no_svlen("3", position, ref, alt).end
            return Variant.objects.create(locus=locus, alt=alt_seq, end=end)

        # #1708 - an explicit variant is an SV on the length of its longer allele
        cls.explicit_expected = {
            _explicit(128300000, "G" * 999, "CT"): STANDARD,
            _explicit(128310000, "G" * 1000, "CT"): STRUCTURAL,
            _explicit(128320000, "G", "GC" * 499 + "T"): STANDARD,
            _explicit(128330000, "G", "GC" * 500): STRUCTURAL,
        }

    def _assert_routing_agrees(self, expected: dict, explicit_sv_min_size):
        qs = Variant.objects.filter(pk__in=[v.pk for v in expected])
        for pipeline_type in (STANDARD, STRUCTURAL):
            expected_ids = {v.pk for v, pt in expected.items() if pt == pipeline_type}
            q = pipeline_type_variant_q(pipeline_type, self.SV_MIN_SIZE, explicit_sv_min_size)
            self.assertEqual(expected_ids, set(qs.filter(q).values_list("pk", flat=True)), pipeline_type)
        for variant, pipeline_type in expected.items():
            self.assertEqual(pipeline_type,
                             pipeline_type_for_variant(variant, self.SV_MIN_SIZE, explicit_sv_min_size), variant)

    def test_variant_and_queryset_routing_agree(self):
        self._assert_routing_agrees(self.expected | self.explicit_expected, self.SV_MIN_SIZE)

    def test_explicit_stays_standard_on_version_without_cut_off(self):
        """ Versions made before #1708 have no explicit cut-off, so nothing they annotated moves """
        self._assert_routing_agrees(dict.fromkeys(self.explicit_expected, STANDARD), None)

    def test_get_for_variant_routes_by_each_versions_cut_off(self):
        """ A version pinned at a lower cut-off sent the short del to the SV pipeline """
        vav = self.annotation_version.variant_annotation_version
        VariantAnnotationVersion.objects.filter(pk=vav.pk).update(structural_variant_min_size=50)
        lock = AnnotationRangeLock.objects.create(version=vav, min_variant=self.short_del,
                                                  max_variant=self.short_del, count=1)
        for pipeline_type in (STANDARD, STRUCTURAL):
            AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type)
        self.assertEqual(STRUCTURAL, AnnotationRun.get_for_variant(self.short_del, self.grch37).pipeline_type)

    def test_dump_writes_short_symbolic_explicit(self):
        qs = Variant.objects.filter(pk__in=[self.short_del.pk, self.long_del.pk])
        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf_filename = os.path.join(tmp_dir, "dump.vcf")
            write_qs_to_vcf(vcf_filename, self.grch37, qs, structural_variant_min_size=self.SV_MIN_SIZE)
            records = {r.INFO["variant_id"]: r for r in cyvcf2.VCF(vcf_filename)}

        short = records[self.short_del.pk]
        self.assertEqual((61, 1), (len(short.REF), len(short.ALT[0])))
        self.assertEqual([VCFSymbolicAllele.DEL], records[self.long_del.pk].ALT)
