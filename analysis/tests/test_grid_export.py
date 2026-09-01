"""
Tests for the analysis node CSV/VCF export path (issue #1257).

Builds a small but real node export - a SampleNode over a cohort whose VCF has FT filters, with
variants spread over several contigs - and runs it through node_grid_get_export_iterator, the same
entry point the view and the Celery export tasks use.
"""
import csv
import json
import tempfile
import zipfile
from datetime import datetime, timezone
from unittest import mock
from urllib.parse import urlencode

from django.contrib.auth.models import User
from django.db import connection
from django.test import SimpleTestCase, TestCase, override_settings
from django.test.client import Client
from django.test.utils import CaptureQueriesContext
from django.urls.base import resolve, reverse

from analysis.grid_export import format_items_iterator, node_grid_get_export_iterator
from analysis.grids import ExportVariantGrid
from analysis.models import Analysis
from analysis.models.enums import NodeStatus
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.models.models_variant_tag import VariantTag
from analysis.tasks.analysis_grid_export_tasks import (
    NODE_EXPORT_GENERATOR,
    export_node_to_downloadable_file,
)
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils import FakeRequest
from library.django_utils.grid_export import EXPORT_ROWS_PER_CHUNK, grid_export_csv
from snpdb.models import CachedGeneratedFile, CohortGenotype, GenomeBuild, Tag
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.models.models_enums import CohortGenotypeCollectionType
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

# (contig name, position) - deliberately out of both PK order and contig-name string order, so an
# export that leaks either shows up as a difference
VARIANT_COORDINATES = [
    ("3", 3000),
    ("1", 2000),
    ("10", 1000),
    ("1", 1000),
    ("2", 5000),
    ("10", 500),
]


class TestCsvChunking(SimpleTestCase):
    """ (A5) Rows are buffered into chunks rather than yielded one WSGI chunk at a time """

    def test_rows_are_buffered_into_chunks(self):
        csv_columns = [{"name": "id", "label": "id"}, {"name": "value", "label": "value"}]
        num_rows = EXPORT_ROWS_PER_CHUNK * 2 + 5
        items = ({"id": i, "value": f"row_{i}"} for i in range(num_rows))

        chunks = list(grid_export_csv(csv_columns, items))
        # header, 2 full chunks, then the remainder
        self.assertEqual(len(chunks), 4)
        lines = "".join(chunks).splitlines()
        self.assertEqual(len(lines), num_rows + 1)
        self.assertEqual(lines[-1], f"{num_rows - 1},row_{num_rows - 1}")


class GridExportTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_1257_export')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)

        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.vcf = cls.cohort.vcf
        # FT is only offered as a grid column when the VCF declares a sample filters field
        cls.vcf.sample_filters_field = "FT"
        cls.vcf.save()
        cls.sample = cls.cohort.cohortsample_set.get(sample__name='proband').sample

        cgc = CohortGenotypeCollection.objects.get(cohort=cls.cohort, cohort_version=cls.cohort.version,
                                                   collection_type=CohortGenotypeCollectionType.UNCOMMON)
        cls.variants = []
        for chrom, position in VARIANT_COORDINATES:
            variant = slowly_create_test_variant(chrom, position, "A", "T", cls.genome_build)
            cls.variants.append(variant)
            cls._add_genotype(cgc, variant)

        # One variant carries the missing sentinel in every numeric sample column, so exports can be
        # checked to render it as '.'
        cls.missing_variant = slowly_create_test_variant("2", 9000, "A", "T", cls.genome_build)
        cls._add_genotype(cgc, cls.missing_variant, read_depth=CohortGenotype.MISSING_NUMBER_VALUE)

        cls.analysis = Analysis(genome_build=cls.genome_build)
        cls.analysis.set_defaults_and_save(cls.user)

    @classmethod
    def _add_genotype(cls, cgc, variant, read_depth=40):
        """ Insert a CohortGenotype row into the collection's partition table (proband at index 0) """
        old_db_table = CohortGenotype._meta.db_table
        try:
            CohortGenotype._meta.db_table = cgc.get_partition_table()
            CohortGenotype.objects.create(
                collection=cgc, variant=variant,
                ref_count=0, het_count=1, hom_count=0, unk_count=0,
                filters="X",
                samples_zygosity="E..",
                samples_allele_depth=[20, 0, 0],
                samples_allele_frequency=[0.5, 0.0, 0.0],
                samples_read_depth=[read_depth, 0, 0],
                samples_genotype_quality=[30, 30, 30],
                samples_phred_likelihood=[0, 0, 0],
            )
        finally:
            CohortGenotype._meta.db_table = old_db_table

    def _sample_node(self) -> SampleNode:
        node = SampleNode.objects.create(analysis=self.analysis, sample=self.sample)
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()
        return node

    def _request(self, position_less_than=None, rules=None) -> FakeRequest:
        request = FakeRequest(user=self.user)
        if position_less_than is not None:
            rules = [{"op": "lt", "field": "locus__position", "data": str(position_less_than)}]
        if rules:
            request.GET = {"filters": json.dumps({"groupOp": "AND", "rules": rules})}
        return request

    def _export_lines(self, node, export_type="csv", request=None) -> list[str]:
        if request is None:
            request = self._request()
        _filename, file_iterator = node_grid_get_export_iterator(request, node, export_type)
        return "".join(file_iterator).splitlines()

    def _export_csv(self, node, request=None) -> tuple[list[str], list[list[str]]]:
        """ (header, rows) of a CSV export """
        rows = list(csv.reader(self._export_lines(node, request=request)))
        return rows[0], rows[1:]


class TestExportPerRowQueries(GridExportTestCase):
    """ (A1) The VCF's filter code lookup is per-VCF, so it must not scale with the number of rows """

    def _export_and_count_filter_lookups(self, node, request) -> tuple[int, int]:
        with CaptureQueriesContext(connection) as queries:
            rows = self._export_lines(node, request=request)
        # get_code_lookup selects the VCF's filter rows - other snpdb_vcffilter queries (the node's own
        # filter settings, existence checks) join different tables
        lookups = [q["sql"] for q in queries.captured_queries
                   if 'FROM "snpdb_vcffilter"' in q["sql"] and "filter_code" in q["sql"]]
        return len(rows), len(lookups)

    def test_vcf_filter_lookup_is_not_per_row(self):
        node = self._sample_node()
        all_rows, all_lookups = self._export_and_count_filter_lookups(node, self._request())
        some_rows, some_lookups = self._export_and_count_filter_lookups(node, self._request(1500))

        self.assertEqual(all_rows, node.count + 1)  # header
        self.assertGreater(all_rows, some_rows)  # the two exports really do differ in row count
        self.assertEqual(all_lookups, some_lookups,
                         "VCFFilter code lookups should not scale with the number of exported rows")


class TestFormatItemsIterator(GridExportTestCase):
    """ (A2) format_items_iterator only summarises tags - the missing value substitution is done
        upstream by the sample columns' packed_data_replace """

    def test_tags_are_summarised_and_analysis_tags_applied(self):
        items = [
            {"id": 1, "tags_global": "foo:2024-03-01|bar:2019-07-12|foo:2019-07-12", "tags": None},
            {"id": 2, "tags_global": None, "tags": None},
        ]
        formatted = list(format_items_iterator(iter(items), {2: "in_analysis"}))
        self.assertEqual(formatted[0]["tags_global"], "foo x 2, bar")
        self.assertIsNone(formatted[1]["tags_global"])
        self.assertIsNone(formatted[0]["tags"])
        self.assertEqual(formatted[1]["tags"], "in_analysis")

    def test_tags_fresh_counts_with_stale_cutoff(self):
        payload = "foo:2024-03-01|foo:2019-07-12|foo:2019-06-01|bar:2019-01-01|baz:2024-05-06"
        items = [{"id": 1, "tags_global": payload, "tags": None}]
        stale_date = datetime(2020, 1, 1, tzinfo=timezone.utc)
        formatted = list(format_items_iterator(iter(items), tag_stale_date=stale_date))
        # Multi-event tags show fresh of total; a single entirely-stale tag shows (0 fresh);
        # a single fresh tag stays as today
        self.assertEqual(formatted[0]["tags_global"], "foo x 3 (1 fresh), bar (0 fresh), baz")

    def test_missing_sample_values_export_as_dot(self):
        node = self._sample_node()
        header, rows = self._export_csv(node)
        variant_id_index = header.index("variant_id")
        read_depth_index = header.index("DP proband")

        by_variant_id = {int(row[variant_id_index]): row for row in rows}
        self.assertEqual(by_variant_id[self.missing_variant.pk][read_depth_index], ".")
        other_variant = self.variants[0]
        self.assertEqual(by_variant_id[other_variant.pk][read_depth_index], "40")


class TestFilteredExportRowCount(GridExportTestCase):
    """ (A3) Grid filters must be applied exactly once. A second .filter() over a multi-valued
        relation adds a second JOIN, multiplying the rows that relation already fans out """

    def test_filtered_export_applies_filter_once(self):
        tagged_variant = self.variants[0]
        for tag_id in ["tag_a", "tag_b"]:
            tag = Tag.objects.get_or_create(id=tag_id)[0]
            VariantTag.objects.create(variant=tagged_variant, genome_build=self.genome_build, tag=tag,
                                       analysis=self.analysis, user=self.user)

        node = self._sample_node()
        expected = node.get_queryset().filter(varianttag__tag__isnull=False).count()
        self.assertGreater(expected, 1)  # the relation fans out, so a second JOIN would show up

        request = self._request(rules=[{"op": "nn", "field": "varianttag__tag", "data": ""}])
        header, rows = self._export_csv(node, request=request)
        self.assertEqual(len(rows), expected)
        variant_id_index = header.index("variant_id")
        self.assertEqual({int(row[variant_id_index]) for row in rows}, {tagged_variant.pk})


class TestPkBatchedExport(GridExportTestCase):
    """ (B) Export walks the node's PKs a contig at a time, in genome build contig order, running the
        annotated queryset against batches of them """

    def _expected_pk_order(self, node) -> list[int]:
        pks = []
        for contig in self.genome_build.standard_contigs:
            contig_qs = node.get_queryset().filter(locus__contig=contig).order_by("locus__position", "pk")
            pks.extend(contig_qs.values_list("pk", flat=True))
        return pks

    def _exported_pks(self, node) -> list[int]:
        header, rows = self._export_csv(node)
        variant_id_index = header.index("variant_id")
        return [int(row[variant_id_index]) for row in rows]

    def test_export_is_in_contig_then_position_order(self):
        node = self._sample_node()
        expected = self._expected_pk_order(node)
        self.assertEqual(len(expected), node.count)
        # the fixture spans several contigs, so contig ordering is actually exercised
        contigs = {v.locus.contig_id for v in self.variants}
        self.assertGreater(len(contigs), 1)

        self.assertEqual(self._exported_pks(node), expected)

    def test_batch_size_does_not_change_output(self):
        node = self._sample_node()
        expected = self._exported_pks(node)
        # small enough that every contig's PKs span multiple batches
        with mock.patch.object(ExportVariantGrid, "EXPORT_PK_BATCH_SIZE", 1):
            self.assertEqual(self._exported_pks(node), expected)


class TestNodeExportLaunch(GridExportTestCase):
    """ (C) The view launches the export under Celery and redirects to the CachedGeneratedFile poll URL.
        get_or_create_and_launch de-dupes on the params hash, so a double-click can't start two of them,
        while a differently filtered view of the same node is a separate download """

    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.client = Client()
        self.client.force_login(self.user)

    def _grid_params(self, **extra_params) -> dict:
        """ The whitelisted params the view passes on to the export task """
        params = {
            "node_id": str(self.node.pk),
            "version_id": str(self.node.version),
            "extra_filters": "",
        }
        params.update(extra_params)
        return params

    def _launch_export(self, **extra_params) -> CachedGeneratedFile:
        params = dict(self._grid_params(**extra_params), export_type="csv")
        url = reverse("node_grid_export", kwargs={"analysis_id": self.analysis.pk}) + "?" + urlencode(params)
        with mock.patch.object(export_node_to_downloadable_file, "apply_async") as mock_apply_async:
            mock_apply_async.return_value = mock.Mock(id="fake-celery-task-id", result=None)
            response = self.client.get(url)
        self.assertEqual(response.status_code, 302)
        cgf_id = resolve(response.url).kwargs["cgf_id"]
        return CachedGeneratedFile.objects.get(pk=cgf_id)

    def _num_cached_files(self) -> int:
        return CachedGeneratedFile.objects.filter(generator=NODE_EXPORT_GENERATOR).count()

    def test_identical_launch_reuses_cached_generated_file(self):
        first = self._launch_export()
        second = self._launch_export()
        self.assertEqual(first.pk, second.pk)
        self.assertEqual(self._num_cached_files(), 1)

    def test_stale_version_id_exports_latest_version(self):
        """ A node update (eg a tag delete) bumps the version while the grid is open - the export
            uses the latest version rather than rejecting the client's stale one (#789) """
        stale_version = self.node.version
        SampleNode.objects.filter(pk=self.node.pk).update(version=stale_version + 1)
        stale = self._launch_export(version_id=str(stale_version))
        self.node.refresh_from_db()
        current = self._launch_export()
        self.assertEqual(stale.pk, current.pk)

    def test_different_grid_filters_are_separate_downloads(self):
        unfiltered = self._launch_export()
        filters = json.dumps({"groupOp": "AND",
                              "rules": [{"op": "lt", "field": "locus__position", "data": "1500"}]})
        filtered = self._launch_export(filters=filters)
        self.assertNotEqual(unfiltered.pk, filtered.pk)
        self.assertEqual(self._num_cached_files(), 2)

    def test_export_task_writes_zipped_csv(self):
        cgf = self._launch_export()
        # The view strips version_id from the grid params it hands the task (node.version is authoritative)
        task_grid_params = self._grid_params()
        task_grid_params.pop("version_id")
        with tempfile.TemporaryDirectory() as generated_dir:
            with override_settings(GENERATED_DIR=generated_dir):
                export_node_to_downloadable_file.apply(args=(self.node.pk, self.node.version, self.user.pk,
                                                             "csv", None, task_grid_params))
                cgf.refresh_from_db()
                self.assertEqual(cgf.task_status, "SUCCESS")
                self.assertEqual(cgf.progress, 1)
                self.assertTrue(cgf.filename.endswith(".csv.zip"), cgf.filename)
                with zipfile.ZipFile(cgf.filename) as zipf:
                    csv_name = zipf.namelist()[0]
                    lines = zipf.read(csv_name).decode().splitlines()
                self.assertEqual(len(lines), self.node.count + 1)  # header
