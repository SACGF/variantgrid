"""
Backfilling VariantAnnotation columns from an annotated VCF (#1675).

@see claude/plans/1675_backfill_annotation_column_plan.md
"""
import os
import tempfile

from django.test import TestCase
from django.test.utils import override_settings

from annotation.backfill_columns import (
    BackfillColumnError,
    dump_annotated_variants,
    import_backfill_vcf,
    normalise_source_value,
    resolve_backfill_columns,
)
from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import AnnotationVersion, VariantAnnotation, VariantAnnotationVersion
from annotation.models.models import AnnotationRangeLock, AnnotationRun
from annotation.vep_field_formatters import format_pick_highest_int
from genes.models_enums import AnnotationConsortium
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

COSMIC_V99 = 99
COSMIC_V101_SAMPLE_COUNT = "GENOME_SCREEN_SAMPLE_COUNT"


@override_settings(**get_fake_annotation_settings_dict(columns_version=3))
class BackfillColumnResolutionTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        kwargs = get_fake_vep_version(cls.genome_build, AnnotationConsortium.ENSEMBL, 3)
        kwargs["cosmic"] = COSMIC_V99
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs,
                                                          status=VariantAnnotationVersion.Status.ACTIVE)

    def test_registry_supplies_source_field_and_formatter(self):
        target, = resolve_backfill_columns(self.vav, ["cosmic_count"])
        # COSMIC v99 writes the count as SAMPLE_COUNT, under the --custom label prefix
        self.assertEqual(target.source_field, "COSMIC_SAMPLE_COUNT")
        self.assertFalse(target.from_csq)
        self.assertEqual(target.formatter, format_pick_highest_int)
        # cosmic_count is variant-level, so only the representative transcript table holds it
        self.assertEqual(target.base_table_names,
                         (VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION,))

    def test_info_override_wins(self):
        target, = resolve_backfill_columns(self.vav, ["cosmic_count"],
                                           info_overrides={"cosmic_count": COSMIC_V101_SAMPLE_COUNT})
        self.assertEqual(target.source_field, COSMIC_V101_SAMPLE_COUNT)
        self.assertEqual(target.formatter, format_pick_highest_int)

    def test_csq_override_targets_representative_table_only(self):
        target, = resolve_backfill_columns(self.vav, ["impact"], csq_overrides={"impact": "IMPACT"})
        self.assertTrue(target.from_csq)
        self.assertEqual(target.base_table_names,
                         (VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION,))

    def test_column_outside_columns_version_reported(self):
        # cadd_raw arrived in columns_version 4
        with self.assertRaises(BackfillColumnError) as cm:
            resolve_backfill_columns(self.vav, ["cadd_raw"])
        self.assertIn("cadd_raw", str(cm.exception))

    def test_unknown_column_reported(self):
        with self.assertRaises(BackfillColumnError) as cm:
            resolve_backfill_columns(self.vav, ["not_a_column"])
        self.assertIn("not_a_column", str(cm.exception))

    def test_override_for_column_not_being_backfilled_reported(self):
        with self.assertRaises(BackfillColumnError):
            resolve_backfill_columns(self.vav, ["cosmic_count"],
                                     info_overrides={"cosmic_legacy_id": "LEGACY_ID"})


class NormaliseSourceValueTests(TestCase):
    def test_multi_valued_info_formats_as_vep_joined_string(self):
        """ cyvcf2 hands back a tuple where VEP would have written 12&7 - the registry formatter has to
            see the same thing either way """
        self.assertEqual(normalise_source_value((12, 7)), "12&7")
        self.assertEqual(format_pick_highest_int(normalise_source_value((12, 7))), 12)

    def test_single_value(self):
        self.assertEqual(normalise_source_value(12), "12")
        self.assertEqual(format_pick_highest_int(normalise_source_value(12)), 12)

    def test_empty_values(self):
        for empty in (None, "", ".", " "):
            self.assertIsNone(normalise_source_value(empty), empty)

    def test_info_flag(self):
        self.assertEqual(normalise_source_value(True), "1")


@override_settings(**get_fake_annotation_settings_dict(columns_version=3))
class BackfillColumnRoundTripTests(TestCase):
    """ Dump an annotated version, hand back a VCF carrying a COSMIC count, import it. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        kwargs = get_fake_vep_version(cls.genome_build, AnnotationConsortium.ENSEMBL, 3)
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs,
                                                          status=VariantAnnotationVersion.Status.ACTIVE)
        AnnotationVersion.latest(cls.genome_build, validate=False)

        cls.variants = [
            slowly_create_test_variant("1", 100 + i, "A", "G", cls.genome_build) for i in range(3)
        ]
        annotation_range_lock = AnnotationRangeLock.objects.create(version=cls.vav,
                                                                   min_variant=cls.variants[0],
                                                                   max_variant=cls.variants[-1],
                                                                   count=len(cls.variants))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock)
        # The pipeline COPYs annotation straight into the version's partition (the base table is
        # partitioned by INHERITS, with no routing trigger), and that's the table the backfill updates
        partition = cls.vav.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, partition):
            for variant in cls.variants:
                VariantAnnotation.objects.create(version=cls.vav, variant=variant,
                                                 annotation_run=annotation_run,
                                                 symbol="RUNX1", cosmic_legacy_id="COSM1")
        # An existing value the file has nothing to say about - what --clear-missing decides the fate of
        VariantAnnotation.objects.filter(version=cls.vav, variant=cls.variants[2]).update(cosmic_count=3)

    def setUp(self):
        self.targets = resolve_backfill_columns(self.vav, ["cosmic_count"],
                                                info_overrides={"cosmic_count": COSMIC_V101_SAMPLE_COUNT})

    def _cosmic_counts(self) -> list:
        qs = VariantAnnotation.objects.filter(version=self.vav).order_by("variant_id")
        return list(qs.values_list("cosmic_count", flat=True))

    @staticmethod
    def _annotate_dump(dump_filename: str, annotated_filename: str, counts: dict[int, str]):
        """ Stand in for the operator's `bcftools annotate`: add a COSMIC count to the dumped records. """
        info_header = (f'##INFO=<ID={COSMIC_V101_SAMPLE_COUNT},Number=1,Type=Integer,'
                       f'Description="How many genome screens samples have this mutation">\n')
        with open(dump_filename) as f_in, open(annotated_filename, "w") as f_out:
            for line in f_in:
                if line.startswith("#CHROM"):
                    f_out.write(info_header)
                if line.startswith("#"):
                    f_out.write(line)
                    continue
                columns = line.rstrip("\n").split("\t")
                variant_id = int(columns[7].split("variant_id=")[1].split(";")[0])
                if count := counts.get(variant_id):
                    columns[7] += f";{COSMIC_V101_SAMPLE_COUNT}={count}"
                f_out.write("\t".join(columns) + "\n")

    def _dump_and_annotate(self, tmp_dir: str, counts: dict[int, str], only_missing: bool = False) -> str:
        dump_filename = os.path.join(tmp_dir, "dump.vcf")
        count = dump_annotated_variants(self.vav, dump_filename, targets=self.targets,
                                        only_missing=only_missing)
        self.assertEqual(count, 2 if only_missing else 3)
        annotated_filename = os.path.join(tmp_dir, "annotated.vcf")
        self._annotate_dump(dump_filename, annotated_filename, counts)
        return annotated_filename

    def test_round_trip_writes_formatted_values(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            # The '&'-joined form COSMIC produces where several records overlap a variant
            counts = {self.variants[0].pk: "12&7", self.variants[1].pk: "5"}
            annotated_filename = self._dump_and_annotate(tmp_dir, counts)
            report = import_backfill_vcf(self.vav, annotated_filename, self.targets)

        self.assertEqual(report["records_read"], 3)
        # Third variant keeps its existing value - the file said nothing about it
        self.assertEqual(self._cosmic_counts(), [12, 5, 3])
        # Neighbouring columns untouched
        annotation = VariantAnnotation.objects.get(version=self.vav, variant=self.variants[0])
        self.assertEqual(annotation.symbol, "RUNX1")
        self.assertEqual(annotation.cosmic_legacy_id, "COSM1")

    def test_clear_missing_nulls_values_the_file_is_silent_on(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            annotated_filename = self._dump_and_annotate(tmp_dir, {self.variants[0].pk: "12"})
            import_backfill_vcf(self.vav, annotated_filename, self.targets, clear_missing=True)

        self.assertEqual(self._cosmic_counts(), [12, None, None])

    def test_dry_run_changes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            annotated_filename = self._dump_and_annotate(tmp_dir, {self.variants[0].pk: "12"})
            report = import_backfill_vcf(self.vav, annotated_filename, self.targets, dry_run=True)

        self.assertEqual(report["rows_changed"], 1)
        self.assertEqual(self._cosmic_counts(), [None, None, 3])

    def test_only_missing_skips_populated_rows(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            # variants[2] already has a count, so it isn't dumped (asserted in _dump_and_annotate)
            annotated_filename = self._dump_and_annotate(tmp_dir, {self.variants[1].pk: "5"},
                                                         only_missing=True)
            import_backfill_vcf(self.vav, annotated_filename, self.targets)

        self.assertEqual(self._cosmic_counts(), [None, 5, 3])

    def test_not_null_skips_variants_the_source_never_matched(self):
        """ cosmic_legacy_id stands in for a sibling column the same source writes - only variants
            holding one can gain a count, so the rest stay out of the dump """
        VariantAnnotation.objects.filter(version=self.vav, variant=self.variants[1]) \
                                 .update(cosmic_legacy_id=None)
        with tempfile.TemporaryDirectory() as tmp_dir:
            dump_filename = os.path.join(tmp_dir, "dump.vcf")
            count = dump_annotated_variants(self.vav, dump_filename, targets=self.targets,
                                            only_missing=True, not_null_columns=["cosmic_legacy_id"])
        self.assertEqual(count, 1)

    def test_not_null_rejects_a_column_that_is_not_a_field(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            dump_filename = os.path.join(tmp_dir, "dump.vcf")
            with self.assertRaises(BackfillColumnError):
                dump_annotated_variants(self.vav, dump_filename, targets=self.targets,
                                        not_null_columns=["not_a_column"])

    def test_batches_smaller_than_the_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            counts = {v.pk: str(i + 1) for i, v in enumerate(self.variants)}
            annotated_filename = self._dump_and_annotate(tmp_dir, counts)
            report = import_backfill_vcf(self.vav, annotated_filename, self.targets, batch_size=2)

        self.assertEqual(report["rows_copied"], 3)
        self.assertEqual(self._cosmic_counts(), [1, 2, 3])

    def test_record_without_variant_id_reported(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            annotated_filename = self._dump_and_annotate(tmp_dir, {})
            stripped_filename = os.path.join(tmp_dir, "no_variant_id.vcf")
            with open(annotated_filename) as f_in, open(stripped_filename, "w") as f_out:
                for line in f_in:
                    if not line.startswith("#"):
                        columns = line.rstrip("\n").split("\t")
                        columns[7] = "."
                        line = "\t".join(columns) + "\n"
                    f_out.write(line)

            with self.assertRaises(BackfillColumnError):
                import_backfill_vcf(self.vav, stripped_filename, self.targets)
