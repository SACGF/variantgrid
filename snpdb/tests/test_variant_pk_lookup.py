import tempfile

from django.test import TestCase

from snpdb.models import GenomeBuild, Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.variant_pk_lookup import VariantPKLookup

LONG_ALTS = ["ACGTACGTACGTACGTAC", "TTTTTGGGGGCCCCCAAAAA", "GATTACAGATTACA"]


class VariantPKLookupTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        cls.long_variants = [slowly_create_test_variant("1", 1000 + i, "A", alt, cls.genome_build)
                             for i, alt in enumerate(LONG_ALTS)]

    def test_long_sequences_resolved_in_one_query(self):
        variant_pk_lookup = VariantPKLookup(self.genome_build)
        for i, alt in enumerate(LONG_ALTS):
            variant_pk_lookup.add(VariantCoordinate.from_explicit_no_svlen("1", 1000 + i, "A", alt))
        novel = VariantCoordinate.from_explicit_no_svlen("1", 2000, "A", "CCCCCCCCCCGGGGGGGGGG")
        variant_pk_lookup.add(novel)

        with self.assertNumQueries(2):  # Long sequences, then variants for the contig
            variant_pk_lookup.batch_check()
        self.assertEqual({v.pk for v in self.long_variants}, set(variant_pk_lookup.variant_pk_by_hash.values()))
        self.assertEqual([novel], variant_pk_lookup.unknown_variant_coordinates)

    def test_unknown_short_sequence_has_no_hash(self):
        variant_pk_lookup = VariantPKLookup(self.genome_build)
        variant_coordinate = VariantCoordinate.from_explicit_no_svlen("1", 3000, "A", "NANAN")
        self.assertIsNone(variant_pk_lookup.get_variant_coordinate_hash(variant_coordinate))
        self.assertIsNone(variant_pk_lookup.add(variant_coordinate))
        self.assertEqual([variant_coordinate], variant_pk_lookup.unknown_variant_coordinates)

    def test_insert_unknown_long_and_short_sequences(self):
        variant_coordinates = [
            VariantCoordinate.from_explicit_no_svlen("1", 4000, "A", "NANAN"),
            VariantCoordinate.from_explicit_no_svlen("1", 4001, "GATTACAGATTACA", "G"),  # Existing long seq as ref
            VariantCoordinate.from_explicit_no_svlen("1", 4002, "A", "CCCCCCCCCCGGGGGGGGGGT"),
        ]
        with tempfile.TemporaryDirectory() as working_dir:
            variant_pk_lookup = VariantPKLookup(self.genome_build, working_dir=working_dir)
            for variant_coordinate in variant_coordinates:
                variant_pk_lookup.add(variant_coordinate)
            variant_pk_lookup.batch_check(insert_unknown=True)

        variant_pk_lookup = VariantPKLookup(self.genome_build)
        variant_hashes = [variant_pk_lookup.get_variant_coordinate_hash(vc) for vc in variant_coordinates]
        variant_ids = variant_pk_lookup.get_variant_ids(variant_hashes)
        for variant_coordinate, variant_id in zip(variant_coordinates, variant_ids):
            self.assertEqual(variant_coordinate, Variant.objects.get(pk=variant_id).coordinate)
