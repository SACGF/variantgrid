from django.test import TestCase

from annotation import vep_columns
from annotation.enums import VariantAnnotationPipelineType
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import SVOverlapProcessor


class SVOverlapProcessorTest(TestCase):
    """ An SV can overlap several gnomAD-SV records - one of them supplies the gnomAD columns the
        analysis PopulationNode filters on, while the gnomad_sv_overlap_* fields keep them all """

    @staticmethod
    def _process(raw_db_data: dict, same_type: bool = False) -> dict:
        """ The version's snapshot, not settings, decides the filter and pick """
        cvf_list = vep_columns.filter_for(genome_build_name="GRCh37",
                                          pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT)
        variant_annotation_version = VariantAnnotationVersion(sv_overlap_same_type=same_type,
                                                              sv_overlap_single_value_method="lowest_af")
        SVOverlapProcessor(cvf_list, variant_annotation_version).process(raw_db_data)
        return raw_db_data

    def test_lowest_af_compares_numerically(self):
        """ gnomAD-SV mixes '4.6e-05' and '0.006085' notation, which sort the wrong way as strings """
        raw_db_data = self._process({
            "variant_class": "deletion",
            "gnomad_sv_overlap_name": "gnomAD-SV_v2.1_DEL_4_44027&gnomAD-SV_v2.1_DEL_4_44079",
            "gnomad_sv_overlap_af": "4.6e-05&0.006085",
            "gnomad_sv_overlap_percent": "94.5&12.1",
            "gnomad_sv_overlap_coords": "4:34558060-36509366&4:35146800-35152239",
            "gnomad_af": "4.6e-05&0.006085",
            "gnomad_ac": "1&132",
            "gnomad_popmax_af": "0.000131&0.015789",
        })

        self.assertEqual("4.6e-05", raw_db_data["gnomad_af"])
        self.assertEqual("1", raw_db_data["gnomad_ac"])
        self.assertEqual("0.000131", raw_db_data["gnomad_popmax_af"])
        # Every overlapping record is kept for display
        self.assertEqual("4.6e-05&0.006085", raw_db_data["gnomad_sv_overlap_af"])

    def test_same_type_filter_applies_before_pick(self):
        """ A DUP overlapping a rarer DEL takes the DUP's frequencies """
        raw_db_data = self._process({
            "variant_class": "duplication",
            "gnomad_sv_overlap_name": "gnomAD-SV_v2.1_DEL_17_160435&gnomAD-SV_v2.1_DUP_17_43398",
            "gnomad_sv_overlap_af": "4.6e-05&0.001383",
            "gnomad_sv_overlap_percent": "88.2&97.0",
            "gnomad_sv_overlap_coords": "17:6279501-6326805&17:6294033-6314646",
            "gnomad_af": "4.6e-05&0.001383",
        }, same_type=True)

        self.assertEqual("0.001383", raw_db_data["gnomad_af"])
        self.assertEqual("0.001383", raw_db_data["gnomad_sv_overlap_af"])

    def test_no_matching_records_blanks_columns(self):
        raw_db_data = self._process({
            "variant_class": "deletion",
            "gnomad_sv_overlap_name": "",
            "gnomad_sv_overlap_af": "",
            "gnomad_sv_overlap_percent": "",
            "gnomad_sv_overlap_coords": "",
            "gnomad_af": "",
        })

        self.assertIsNone(raw_db_data["gnomad_af"])
