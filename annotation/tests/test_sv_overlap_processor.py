from django.test import TestCase

from annotation import vep_columns
from annotation.enums import VariantAnnotationPipelineType
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import SVOverlapProcessor
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import VariantCoordinate

# A 10 kb deletion, chr4:35146800-35156800. VEP has already dropped the gnomAD SVs it covers less than 80% of
OUR_SV = VariantCoordinate(chrom="4", position=35146800, ref="N", alt=VCFSymbolicAllele.DEL, svlen=-10000)
SAME = "4:35146800-35156800"
SLIGHTLY_BIGGER = "4:35146000-35157500"  # 11.5 kb: we cover 87% of it
INSIDE = "4:35146800-35148800"  # 2 kb: covers 20% of ours
AROUND = "4:35100000-35200000"  # 100 kb: we cover 10% of it


class SVOverlapProcessorTest(TestCase):
    """ An SV can overlap several gnomAD-SV records - one of them supplies the gnomAD columns the
        analysis PopulationNode filters on, while the gnomad_sv_overlap_* fields keep them all (#2589) """

    @staticmethod
    def _process(raw_db_data: dict, same_type: bool = False) -> dict:
        """ The version's snapshot, not settings, decides the filter and pick """
        cvf_list = vep_columns.filter_for(genome_build_name="GRCh37",
                                          pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT)
        variant_annotation_version = VariantAnnotationVersion(sv_overlap_same_type=same_type,
                                                              sv_overlap_single_value_method="lowest_af",
                                                              sv_overlap_min_fraction=0.8)
        SVOverlapProcessor(cvf_list, variant_annotation_version).process(raw_db_data, OUR_SV)
        return raw_db_data

    @staticmethod
    def _raw(names: str, afs: str, coords: str, **extra) -> dict:
        """ As the inserter sees it: VEP's coords custom lands in both the coords and percent columns """
        return {
            "variant_class": "deletion",
            "gnomad_sv_overlap_name": names,
            "gnomad_sv_overlap_af": afs,
            "gnomad_sv_overlap_percent": coords,
            "gnomad_sv_overlap_coords": coords,
            "gnomad_af": afs,
            **extra,
        }

    def test_lowest_af_compares_numerically(self):
        """ gnomAD-SV mixes '4.6e-05' and '0.006085' notation, which sort the wrong way as strings """
        raw_db_data = self._process(self._raw("gnomAD-SV_v2.1_DEL_4_44027&gnomAD-SV_v2.1_DEL_4_44079",
                                              "4.6e-05&0.006085", f"{SAME}&{SLIGHTLY_BIGGER}",
                                              gnomad_ac="1&132", gnomad_popmax_af="0.000131&0.015789"))

        self.assertEqual("4.6e-05", raw_db_data["gnomad_af"])
        self.assertEqual("1", raw_db_data["gnomad_ac"])
        self.assertEqual("0.000131", raw_db_data["gnomad_popmax_af"])
        # Every overlapping record is kept for display, with the share of our SV it covers
        self.assertEqual("4.6e-05&0.006085", raw_db_data["gnomad_sv_overlap_af"])
        self.assertEqual("100&100", raw_db_data["gnomad_sv_overlap_percent"])

    def test_same_type_filter_applies_before_pick(self):
        """ A DUP overlapping a rarer DEL takes the DUP's frequencies """
        raw_db_data = self._process({
            **self._raw("gnomAD-SV_v2.1_DEL_17_160435&gnomAD-SV_v2.1_DUP_17_43398", "4.6e-05&0.001383",
                        f"{SAME}&{SLIGHTLY_BIGGER}"),
            "variant_class": "duplication",
        }, same_type=True)

        self.assertEqual("0.001383", raw_db_data["gnomad_af"])
        self.assertEqual("0.001383", raw_db_data["gnomad_sv_overlap_af"])

    def test_small_gnomad_sv_inside_ours_dropped(self):
        """ A large SV covers every small gnomAD SV inside it - VEP's cutoff (the gnomAD side) lets them through """
        raw_db_data = self._process(self._raw("gnomAD-SV_v2.1_DEL_4_1&gnomAD-SV_v2.1_DEL_4_2", "0.01&8e-06",
                                              f"{SAME}&{INSIDE}"))

        self.assertEqual("0.01", raw_db_data["gnomad_af"])
        self.assertEqual("gnomAD-SV_v2.1_DEL_4_1", raw_db_data["gnomad_sv_overlap_name"])
        self.assertEqual("100", raw_db_data["gnomad_sv_overlap_percent"])

    def test_large_gnomad_sv_around_ours_dropped(self):
        """ The other way: a common gnomAD SV we are a small part of (what VEP reciprocal=1 lets through) """
        raw_db_data = self._process(self._raw("gnomAD-SV_v2.1_DEL_4_1&gnomAD-SV_v2.1_DEL_4_3", "0.01&0.5",
                                              f"{SAME}&{AROUND}"))

        self.assertEqual("0.01", raw_db_data["gnomad_af"])
        self.assertEqual("gnomAD-SV_v2.1_DEL_4_1", raw_db_data["gnomad_sv_overlap_name"])

    def test_nothing_passes_blanks_columns(self):
        raw_db_data = self._process(self._raw("gnomAD-SV_v2.1_DEL_4_2&gnomAD-SV_v2.1_DEL_4_3", "8e-06&0.5",
                                              f"{INSIDE}&{AROUND}"))

        self.assertIsNone(raw_db_data["gnomad_af"])
        self.assertIsNone(raw_db_data["gnomad_sv_overlap_name"])
        self.assertIsNone(raw_db_data["gnomad_sv_overlap_percent"])

    def test_no_matching_records_blanks_columns(self):
        raw_db_data = self._process(self._raw("", "", ""))

        self.assertIsNone(raw_db_data["gnomad_af"])
