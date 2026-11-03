import copy
import os
from io import StringIO

from django.conf import settings
from django.core.management import call_command
from django.test import TestCase
from django.test.utils import override_settings

from annotation.enums import VariantAnnotationPipelineType
from annotation.fake_data import get_fake_vep_version, retire_seeded_annotation_version
from annotation.models.models import VariantAnnotation
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from genes.enums import AnnotationConsortium
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import GenomeBuild, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant_from_coordinate

TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
SAME = "gnomAD-SV_v2.1_DEL_1_SAME"  # 1:300000-310000, AF 0.01 - a match in both directions
INSIDE = "gnomAD-SV_v2.1_DEL_1_INSIDE"  # 1:300000-302000, AF 1e-05 - 20% of our deletion
BIGGER = "gnomAD-SV_v2.1_DEL_1_BIGGER"  # 1:299000-311000, AF 0.5 - we cover 83% of it
AROUND = "gnomAD-SV_v2.1_DEL_1_AROUND"  # 1:250000-350000, AF 0.3 - we cover 10% of it
COORDS = {SAME: "1:300000-310000", INSIDE: "1:300000-302000", BIGGER: "1:299000-311000",
          AROUND: "1:250000-350000"}
AF = {SAME: "0.01", INSIDE: "1e-05", BIGGER: "0.5", AROUND: "0.3"}


def _gnomad_sv_test_settings() -> dict:
    annotation = copy.deepcopy(settings.ANNOTATION)
    annotation[settings.BUILD_GRCH37]["vep_config"]["gnomad_sv"] = "one_off_sv_overlap_gnomad_sv.vcf"
    return {"ANNOTATION": annotation, "ANNOTATION_VEP_BASE_DIR": TEST_DATA_DIR}


@override_settings(**_gnomad_sv_test_settings())
class OneOffSVOverlapBothDirectionsTest(TestCase):
    """ A 10 kb deletion annotated under the one-way rule, with every overlapping record still on the row """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        retire_seeded_annotation_version(cls.grch37)
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        kwargs["sv_overlap_min_fraction"] = 0.8
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)

        def create_deletion(position):
            vc = VariantCoordinate(chrom="1", position=position, ref="A", alt=VCFSymbolicAllele.DEL, svlen=-10000)
            return slowly_create_test_variant_from_coordinate(vc, cls.grch37)

        cls.variants = {name: create_deletion(300000 + i) for i, name in
                        enumerate(["refilled", "unchanged", "blanked", "pruned"])}
        lock = AnnotationRangeLock.objects.create(version=cls.vav, min_variant=cls.variants["refilled"],
                                                  max_variant=cls.variants["pruned"], count=4)
        sv_run = AnnotationRun.objects.create(annotation_range_lock=lock,
                                              pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT)

        def annotate(variant, names: list[str], percents: list[int]):
            lowest = min(names, key=lambda n: float(AF[n]))
            VariantAnnotation.objects.create(
                version=cls.vav, variant=variant, annotation_run=sv_run,
                gnomad_sv_overlap_name="&".join(names),
                gnomad_sv_overlap_af="&".join(AF[n] for n in names),
                gnomad_sv_overlap_coords="&".join(COORDS[n] for n in names),
                gnomad_sv_overlap_percent="&".join(str(p) for p in percents),
                gnomad_af=float(AF[lowest]), gnomad_ac=1 if lowest == INSIDE else 42)

        annotate(cls.variants["refilled"], [SAME, INSIDE], [100, 20])  # INSIDE's AF was taken
        annotate(cls.variants["unchanged"], [SAME, BIGGER], [100, 100])
        annotate(cls.variants["blanked"], [INSIDE], [20])
        annotate(cls.variants["pruned"], [AROUND, SAME], [100, 100])  # SAME's AF was taken, AROUND only displayed

    def _annotation(self, name) -> VariantAnnotation:
        return VariantAnnotation.objects.get(version=self.vav, variant=self.variants[name])

    def test_repairs(self):
        out = StringIO()
        call_command("one_off_sv_overlap_both_directions", stdout=out)
        self.assertIn("4 SV rows checked, 1 pruned, 1 refilled, 1 blanked", out.getvalue())

        refilled = self._annotation("refilled")
        self.assertEqual(refilled.gnomad_sv_overlap_name, SAME)
        self.assertEqual(refilled.gnomad_sv_overlap_percent, "100")
        self.assertAlmostEqual(refilled.gnomad_af, 0.01)
        self.assertEqual(refilled.gnomad_ac, 42)  # from the VCF, not the row
        self.assertEqual(refilled.gnomad_an, 4200)
        self.assertAlmostEqual(refilled.gnomad_popmax_af, 0.02)

        unchanged = self._annotation("unchanged")
        self.assertEqual(unchanged.gnomad_sv_overlap_name, f"{SAME}&{BIGGER}")
        self.assertAlmostEqual(unchanged.gnomad_af, 0.01)

        blanked = self._annotation("blanked")
        self.assertIsNone(blanked.gnomad_sv_overlap_name)
        self.assertIsNone(blanked.gnomad_af)
        self.assertIsNone(blanked.gnomad_ac)

        pruned = self._annotation("pruned")
        self.assertEqual(pruned.gnomad_sv_overlap_name, SAME)
        self.assertEqual(pruned.gnomad_sv_overlap_af, "0.01")
        self.assertAlmostEqual(pruned.gnomad_af, 0.01)
        self.assertEqual(pruned.gnomad_ac, 42)

    def test_dry_run_writes_nothing(self):
        out = StringIO()
        call_command("one_off_sv_overlap_both_directions", dry_run=True, stdout=out)
        self.assertIn("DRY RUN", out.getvalue())
        self.assertIn("1 refilled", out.getvalue())
        self.assertAlmostEqual(self._annotation("refilled").gnomad_af, 1e-05)
