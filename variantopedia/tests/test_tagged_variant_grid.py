from django.conf import settings
from django.contrib.auth.models import User
from django.db import connection
from django.test import RequestFactory, TestCase
from django.test.utils import CaptureQueriesContext
from django.urls import resolve, reverse
from guardian.shortcuts import assign_perm

from analysis.models import Analysis, VariantTag
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version, create_gata2_transcript_version
from library.django_utils import FakeRequest
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import (
    Allele,
    AlleleOrigin,
    GenomeBuild,
    Tag,
    TagColor,
    TagColorsCollection,
    UserGridConfig,
    UserSettingsOverride,
    Variant,
    VariantAllele,
)
from variantopedia.grids import TaggedVariantGrid, VariantTagCountsColumns, VariantTagsColumns


class TaggedVariantGridTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='tagged_variant_grid_user')[0]
        cls.other_user = User.objects.get_or_create(username='tagged_variant_grid_other_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)

        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="SomaticReportable")

        # both_variant carries both tags, artefact_variant only one
        cls.both_variant, cls.artefact_variant, cls.other_user_variant = list(Variant.objects.order_by("pk")[:3])
        cls._tag(cls.both_variant, cls.artefact)
        cls._tag(cls.both_variant, cls.reportable)
        cls._tag(cls.artefact_variant, cls.artefact)

        # Initial group read permissions depend on deployment settings, so grant explicitly
        other_user_tag = cls._tag(cls.other_user_variant, cls.artefact, user=cls.other_user)
        assign_perm(VariantTag.get_read_perm(), cls.user, other_user_tag)

    @classmethod
    def _tag(cls, variant: Variant, tag: Tag, user: User = None) -> VariantTag:
        allele, _ = VariantAllele.objects.get_or_create(
            variant=variant, genome_build=cls.genome_build, origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            defaults={"allele": Allele.objects.create()})
        return VariantTag.objects.create(variant=variant, allele=allele.allele, tag=tag,
                                         genome_build=cls.genome_build, user=user or cls.user)

    def _grid_rows(self, extra_filters=None) -> dict[int, dict]:
        request = FakeRequest(user=self.user)
        grid = TaggedVariantGrid(request, self.genome_build.name, extra_filters=extra_filters)
        qs = grid.get_initial_queryset().values(*grid.value_columns())
        return {row["id"]: row for row in qs}

    def _grid_variant_ids(self, extra_filters) -> set[int]:
        return set(self._grid_rows(extra_filters))

    def _tags_grid_variant_ids(self, filters) -> set[int]:
        url = reverse('variant_tags_datatable', kwargs={"genome_build_name": self.genome_build.name})
        request = RequestFactory().get(url, filters)
        request.resolver_match = resolve(url)
        request.user = self.user
        config = VariantTagsColumns(request)
        qs = config.filter_queryset(config.get_initial_queryset())
        return set(qs.values_list("variant__id", flat=True))

    def test_single_tag_filter(self):
        self.assertEqual(self._grid_variant_ids({"tag": "Artefact"}),
                         {self.both_variant.pk, self.artefact_variant.pk, self.other_user_variant.pk})

    def test_multiple_tags_require_all(self):
        """ The co-occurrence card links here expecting variants carrying every tag, not any of them """
        self.assertEqual(self._grid_variant_ids({"tags": ["Artefact", "SomaticReportable"]}),
                         {self.both_variant.pk})

    def test_tag_count_column(self):
        rows = self._grid_rows()
        self.assertEqual(rows[self.both_variant.pk]["tag_count"], 2)
        self.assertEqual(rows[self.artefact_variant.pk]["tag_count"], 1)

    def test_user_filter(self):
        """ The user page links here to show a single user's tags """
        self.assertEqual(self._grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})
        self.assertEqual(self._tags_grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})

    def test_variant_tags_datatable_row(self):
        """ Rows carry what the client renderers need - tag pill, variant link, delete URL, constant build """
        self.client.force_login(self.user)
        url = reverse('variant_tags_datatable', kwargs={"genome_build_name": self.genome_build.name})
        response = self.client.get(url, {"tag": self.reportable.pk})
        self.assertEqual(response.status_code, 200)

        rows = response.json()["data"]
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["genome_build"], self.genome_build.name)
        self.assertEqual(row["tag"], {"tag": self.reportable.pk, "variant_id": self.both_variant.pk})
        self.assertEqual(row["variant_string"]["url"],
                         reverse("view_variant", kwargs={"variant_id": self.both_variant.pk}))
        self.assertTrue(row["delete"], "Tag owner can delete their own tag")

    def test_variant_tags_datatable_classify_button(self):
        """ A RequiresClassification tag is a to-do item - offer to complete it from the analysis it came from """
        annotation_version = AnnotationVersion.latest(self.genome_build)
        analysis = Analysis.objects.create(genome_build=self.genome_build, annotation_version=annotation_version,
                                           user=self.user)
        requires_classification, _ = Tag.objects.get_or_create(pk=settings.TAG_REQUIRES_CLASSIFICATION)
        variant_tag = self._tag(self.artefact_variant, requires_classification)
        variant_tag.analysis = analysis
        variant_tag.save()

        self.client.force_login(self.user)
        url = reverse('variant_tags_datatable', kwargs={"genome_build_name": self.genome_build.name})
        response = self.client.get(url, {"tag": requires_classification.pk})
        row = response.json()["data"][0]
        self.assertEqual(row["variant_string"]["classify_url"],
                         reverse("create_classification_for_variant_tag",
                                 kwargs={"analysis_id": analysis.pk, "variant_tag_id": variant_tag.pk}))

    def test_tag_awaiting_liftover_keeps_its_own_coordinate(self):
        """ A tag gets its allele assigned asynchronously (@see _liftover_variant_tag), so a freshly
            made one has nothing to reach a variant through. It was made on a variant in this build,
            so show that rather than dropping the row - the grid disagreed with the tag counts. """
        variant = self.artefact_variant
        VariantTag.objects.create(variant=variant, allele=None, tag=self.reportable,
                                  genome_build=self.genome_build, user=self.user)

        url = reverse('variant_tags_datatable', kwargs={"genome_build_name": self.genome_build.name})
        request = RequestFactory().get(url, {"tag": self.reportable.pk})
        request.resolver_match = resolve(url)
        request.user = self.user
        config = VariantTagsColumns(request)
        rows = config.filter_queryset(config.get_initial_queryset()).values("variant__id", "variant_string")

        by_variant = {r["variant__id"]: r["variant_string"] for r in rows}
        self.assertIn(variant.pk, by_variant, "Tag with no allele yet should still be listed")
        self.assertEqual(f"{variant.locus.contig.name}:{variant.locus.position} "
                         f"{variant.locus.ref}>{variant.alt}",
                         by_variant[variant.pk])

    def test_grids_and_the_no_coordinate_warning_account_for_every_tag(self):
        """ The page warns how many tags the grids can't place in this build. Counting those by a
            different rule than the grids use left tags counted as both shown and not shown. """
        other_build = GenomeBuild.get_name_or_alias("GRCh38")
        for_build = VariantTag.get_for_build(self.genome_build)
        without_coordinate = VariantTag.objects.exclude(pk__in=for_build.values_list("pk", flat=True))

        self.assertEqual(VariantTag.objects.count(),
                         for_build.count() + without_coordinate.count())
        # A tag made in the other build with no allele yet has no coordinate here, so it is warned about
        stray = VariantTag.objects.create(variant=self.both_variant, allele=None, tag=self.artefact,
                                          genome_build=other_build, user=self.user)
        self.assertIn(stray.pk, set(without_coordinate.values_list("pk", flat=True)))
        self.assertNotIn(stray.pk, set(for_build.values_list("pk", flat=True)))

    def test_variant_tags_export(self):
        """ CSV comes off the queryset - the DataTable itself only ever pages 100 rows at a time """
        self.client.force_login(self.user)
        url = reverse('variant_tags_export', kwargs={"genome_build_name": self.genome_build.name})
        response = self.client.get(url, {"tag": self.reportable.pk})
        self.assertEqual(response.status_code, 200)

        lines = b"".join(response.streaming_content).decode().strip().splitlines()
        self.assertEqual(len(lines), 2, lines)
        self.assertTrue(lines[1].startswith(f"{self.genome_build.name},"), lines[1])
        self.assertIn(self.reportable.pk, lines[1])

    def _variant_tag_counts_tags(self) -> list[str]:
        url = reverse('variant_tag_counts_datatable', kwargs={"variant_id": self.both_variant.pk})
        request = RequestFactory().get(url)
        request.resolver_match = resolve(url)
        request.user = self.user
        config = VariantTagCountsColumns(request)
        tag_column = config.rich_columns[0]
        qs = config.get_initial_queryset().order_by(*tag_column.sort_string(False))
        return list(qs.values_list("tag", flat=True))

    def test_variant_tag_counts_custom_sort_order(self):
        """ The variant page tag table follows the sort order from the user's tag colours collection """
        self.assertEqual(self._variant_tag_counts_tags(), ["Artefact", "SomaticReportable"])

        collection = TagColorsCollection.objects.create(name="sort test colors", user=self.user)
        TagColor.objects.create(collection=collection, tag=self.artefact, rgb="", sort_order=10)
        user_settings_override, _ = UserSettingsOverride.objects.get_or_create(user=self.user)
        user_settings_override.tag_colors = collection
        user_settings_override.save()

        self.assertEqual(self._variant_tag_counts_tags(), ["SomaticReportable", "Artefact"])

    def test_user_filter_overrides_show_group_data(self):
        """ An explicit user filter must still show another user's (permission-visible) tags
            even when the grid config is set to only show your own data """
        for caption in [TaggedVariantGrid.grid_name, VariantTagsColumns.GRID_NAME]:
            config = UserGridConfig.get(self.user, caption)
            config.show_group_data = False
            config.save()

        self.assertEqual(self._grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})
        self.assertEqual(self._tags_grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})


class VariantTagsGridQueryTest(TestCase):
    """ The tags grid pages 100 rows out of a query over every tag in the build, so what it costs per
        row and per annotation version is the whole cost of the page (@see issue #1794) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='variant_tags_grid_user')[0]
        cls.other_user = User.objects.get_or_create(username='variant_tags_grid_other_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.old_annotation_version = get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)
        cls.tag = Tag.objects.create(pk="GridQuery")
        cls.variant, cls.other_variant = list(Variant.objects.order_by("pk")[:2])

    @classmethod
    def _tag(cls, variant: Variant, user: User = None, analysis: Analysis = None) -> VariantTag:
        allele, _ = VariantAllele.objects.get_or_create(
            variant=variant, genome_build=cls.genome_build, origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            defaults={"allele": Allele.objects.create()})
        return VariantTag.objects.create(variant=variant, allele=allele.allele, tag=cls.tag, analysis=analysis,
                                         genome_build=cls.genome_build, user=user or cls.user)

    def _response(self, user: User = None):
        url = reverse('variant_tags_datatable', kwargs={"genome_build_name": self.genome_build.name})
        self.client.force_login(user or self.user)
        response = self.client.post(url, {"length": 100, "tag": self.tag.pk})
        self.assertEqual(response.status_code, 200)
        return response.json()

    def _rows(self, user: User = None) -> list[dict]:
        return self._response(user)["data"]

    def _annotate(self, variant_annotation_version: VariantAnnotationVersion, transcript_version):
        """ Annotation lives in the version's own partition, so write it there not the parent table """
        annotation_range_lock, _ = AnnotationRangeLock.objects.get_or_create(
            version=variant_annotation_version, min_variant=self.variant, max_variant=self.variant, count=1)
        annotation_run, _ = AnnotationRun.objects.get_or_create(annotation_range_lock=annotation_range_lock)
        partition_table = variant_annotation_version.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, partition_table):
            VariantAnnotation.objects.create(version=variant_annotation_version, variant=self.variant,
                                             transcript_version=transcript_version,
                                             gene=transcript_version.gene_version.gene,
                                             annotation_run=annotation_run,
                                             predictions_num_pathogenic=0, predictions_num_benign=0)

    def _new_variant_annotation_version(self) -> VariantAnnotationVersion:
        """ Retire the fake version and stand up a copy of it as the build's current annotation """
        old_pk = self.old_annotation_version.variant_annotation_version_id
        variant_annotation_version = VariantAnnotationVersion.objects.get(pk=old_pk)
        variant_annotation_version.pk = None
        variant_annotation_version.status = VariantAnnotationVersion.Status.NEW
        variant_annotation_version.save()  # Creates the partition to write annotation into

        VariantAnnotationVersion.objects.filter(pk=old_pk).update(
            status=VariantAnnotationVersion.Status.HISTORICAL)
        variant_annotation_version.status = VariantAnnotationVersion.Status.ACTIVE
        variant_annotation_version.save()
        AnnotationVersion.new_sub_version(self.genome_build)
        return variant_annotation_version

    def test_gene_symbol_comes_from_the_build_annotation_version_only(self):
        """ Joining every version of the annotation multiplied the rows the DISTINCT had to sort """
        self._tag(self.variant)
        self._annotate(self.old_annotation_version.variant_annotation_version,
                       create_fake_transcript_version(self.genome_build))  # RUNX1
        self._annotate(self._new_variant_annotation_version(),
                       create_gata2_transcript_version(self.genome_build))  # GATA2

        rows = self._rows()
        self.assertEqual(len(rows), 1, "One row per tag, however many versions the variant is annotated in")
        self.assertEqual(rows[0]["gene_symbol"], "GATA2")

    def test_delete_permissions_cost_the_same_however_many_rows(self):
        """ can_write delegates to the analysis, so checking it per row was 2 Guardian queries a row """
        analysis = Analysis.objects.create(genome_build=self.genome_build, user=self.user)
        self._tag(self.variant, analysis=analysis)
        self._rows()  # Warm the caches the first request of a session fills
        with CaptureQueriesContext(connection) as one_row:
            self._rows()

        other_analysis = Analysis.objects.create(genome_build=self.genome_build, user=self.user)
        self._tag(self.other_variant, analysis=analysis)
        self._tag(self.other_variant, analysis=other_analysis)
        with CaptureQueriesContext(connection) as three_rows:
            rows = self._rows()

        self.assertEqual(len(rows), 3)
        self.assertEqual(len(three_rows.captured_queries), len(one_row.captured_queries))

    def test_only_writable_tags_offer_delete(self):
        other_analysis = Analysis.objects.create(genome_build=self.genome_build, user=self.other_user)
        mine = self._tag(self.variant)
        theirs = self._tag(self.other_variant, user=self.other_user, analysis=other_analysis)
        assign_perm(VariantTag.get_read_perm(), self.user, theirs)
        assign_perm(Analysis.get_read_perm(), self.user, other_analysis)

        delete_by_tag_id = {row["id"]: row["delete"] for row in self._rows()}
        self.assertTrue(delete_by_tag_id[mine.pk], "Own tag is deletable")
        self.assertIsNone(delete_by_tag_id[theirs.pk], "Tag on someone else's analysis is not deletable")

    def test_records_total_skips_counting_the_unfiltered_queryset(self):
        """ The initial queryset is every tag in the build - too expensive to count for a footer """
        self._tag(self.variant)
        VariantTag.objects.create(variant=self.other_variant, tag=Tag.objects.create(pk="OtherGridQuery"),
                                  genome_build=self.genome_build, user=self.user)
        data = self._response()
        self.assertEqual(1, data["recordsFiltered"])
        self.assertEqual(data["recordsFiltered"], data["recordsTotal"])
