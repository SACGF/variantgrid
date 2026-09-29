from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from annotation.gene_level_annotation import annotate_gene_level_run
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotationPipelineType,
    VariantAnnotationVersion,
)
from classification.autopopulate_evidence_keys.evidence_from_variant import (
    _get_gnomad_sv_overlap_note,
    _gnomad_sv_sourced_columns,
    get_evidence_fields_for_variant,
)
from classification.enums import SpecialEKeys
from classification.models import EvidenceKey
from genes.models import HGNC, GeneCopyNumberEventKind, GeneSymbol, HGNCImport, SpliceEvent
from genes.models_enums import HGNCStatus
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from genes.tests.gene_level_test_utils import (
    create_gene_copy_number_event,
    create_splice_event_variant,
)
from snpdb.models import GenomeBuild


class GnomADSVOverlapNoteTest(TestCase):
    """ gnomAD values for SVs come from an overlapping gnomAD-SV record - #1030 """

    def setUp(self):
        self.annotation_version = AnnotationVersion(
            variant_annotation_version=VariantAnnotationVersion(gnomad_sv="2.1"))
        self.variant_values = {
            "variantannotation__gnomad_sv_overlap_name": "gnomAD-SV_v2.1_DEL_1_13595&gnomAD-SV_v2.1_DEL_1_13596",
            "variantannotation__gnomad_sv_overlap_af": "0.013322&0.42",
            "variantannotation__gnomad_sv_overlap_percent": "100&55",
            "variantannotation__gnomad_sv_overlap_coords": "1:245636779-245648007&1:245600000-245700000",
            "variantannotation__gnomad_af": 0.013322,
        }

    def test_note_describes_record_used(self):
        note = _get_gnomad_sv_overlap_note(self.variant_values, self.annotation_version)
        self.assertIn("Used: gnomAD-SV_v2.1_DEL_1_13595 1:245636779-245648007 (11228bp), 100% overlap, AF: 0.013322",
                      note)
        self.assertIn("http://gnomad.broadinstitute.org/variant/DEL_1_13595?dataset=gnomad_sv_r2_1", note)
        self.assertIn("Other overlapping records (1): gnomAD-SV_v2.1_DEL_1_13596", note)

    def test_note_lists_all_records_when_none_matches_gnomad_af(self):
        variant_values = dict(self.variant_values, **{"variantannotation__gnomad_af": None})
        note = _get_gnomad_sv_overlap_note(variant_values, self.annotation_version)
        self.assertNotIn("Used:", note)
        self.assertIn("Overlapping records (2):", note)

    def test_no_note_without_overlap(self):
        self.assertIsNone(_get_gnomad_sv_overlap_note({}, self.annotation_version))

    def test_sourced_columns(self):
        sv_sourced_columns = _gnomad_sv_sourced_columns()
        self.assertIn("gnomad_af", sv_sourced_columns)
        self.assertIn("gnomad_popmax_af", sv_sourced_columns)
        self.assertNotIn("gnomad_sv_overlap_af", sv_sourced_columns)


class GeneLevelGeneSymbolTest(TestCase):
    """ A gene-level variant sits on no transcript, so gene_symbol comes off its resolved identity -
        the CNV caller's MYCL1 is MYCL """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)

        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(1014, "BCR"), (76, "ABL1")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
        # MYCL1 is what the panel's CNV VCF writes; HGNC renamed it MYCL
        GeneSymbol.objects.get_or_create(symbol="MYCL")
        HGNC.objects.create(pk=7553, gene_symbol_id="MYCL", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="MYCL proto-oncogene",
                            previous_symbols="MYCL1")

        # AR is seeded with the junctions the TSO 500 panel reports (genes/migrations/0093)
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=644, gene_symbol_id="AR", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="androgen receptor")

        cls.copy_number_variant = create_gene_copy_number_event(
            "MYCL1", GeneCopyNumberEventKind.GAIN).variant
        cls.fusion_variant = create_gene_fusion("BCR", "ABL1").variant
        cls.splice_variant = create_splice_event_variant("AR", "V7").variant
        cls.unnamed_splice_variant = create_splice_event_variant("AR", "grch37_x_1_2").variant
        cls._annotate_gene_level()

    @classmethod
    def _annotate_gene_level(cls):
        """ Run the GENE_LEVEL pipeline over them, so autopopulate sees what it sees in production """
        vav = cls.annotation_version.variant_annotation_version
        variants = sorted([cls.copy_number_variant, cls.fusion_variant, cls.splice_variant,
                           cls.unnamed_splice_variant], key=lambda v: v.pk)
        range_lock = AnnotationRangeLock.objects.create(version=vav, min_variant=variants[0],
                                                       max_variant=variants[-1], count=len(variants))
        annotate_gene_level_run(AnnotationRun.objects.create(
            annotation_range_lock=range_lock,
            pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL))

    def _autopopulated(self, variant, key: str):
        evidence_keys_list = list(EvidenceKey.objects.all().select_related("variantgrid_column"))
        data = get_evidence_fields_for_variant(self.genome_build, variant, None, None,
                                               evidence_keys_list, self.annotation_version)
        return EvidenceKey.get_value(data.data.get(key))

    def _autopopulated_gene_symbol(self, variant) -> str:
        return self._autopopulated(variant, SpecialEKeys.GENE_SYMBOL)

    def test_gain_autopopulates_the_approved_symbol(self):
        self.assertEqual("MYCL", self._autopopulated_gene_symbol(self.copy_number_variant))

    def test_fusion_autopopulates_the_anchor(self):
        """ The 5' partner - what the Variant is filed under and what the report sorts on """
        self.assertEqual("BCR", self._autopopulated_gene_symbol(self.fusion_variant))

    def test_splice_autopopulates_the_gene_and_the_junctions_name(self):
        """ The seeded SpliceEvent is what the report calls the junction """
        self.assertEqual("AR", self._autopopulated_gene_symbol(self.splice_variant))
        self.assertEqual("AR-V7 splice variant",
                         self._autopopulated(self.splice_variant, SpecialEKeys.SPLICE_LABEL))

    def test_a_junction_we_have_no_name_for_is_printed_as_its_breakpoints(self):
        """ Which reads as raw coordinates on a report, and is the prompt for the scientist to
            replace it with a name """
        self.assertEqual("AR GRCh37 X:1-2 splice",
                         self._autopopulated(self.unnamed_splice_variant, SpecialEKeys.SPLICE_LABEL))
        self.assertTrue(SpliceEvent.objects.exists(), "the seeded junctions are still there")
