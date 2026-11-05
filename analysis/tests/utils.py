"""
Shared fixtures for analysis node tests: AnalysisSetupMixin (a user, GRCh37, the fake annotation version and an
analysis), and the row-membership helpers - create_annotated_variants / annotate_variants write named variants'
VariantAnnotation rows into a version's partition, node_keeps reports which of those names a node's Q keeps.
"""
from typing import Optional

from django.contrib.auth.models import User

from analysis.models import Analysis
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.fake_data import get_fake_annotation_version
from annotation.models.models import VariantAnnotation
from annotation.models.models_annotation_run import AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from library.django_utils.django_partition import bulk_create_in_partition
from snpdb.models import GenomeBuild
from snpdb.models.models_variant import Variant
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class AnalysisSetupMixin:
    """setUpTestData mixin that provides cls.analysis and cls.grch37.

    Shared by multiple analysis test classes to avoid repeating the same
    user + genome build + annotation version + analysis setup boilerplate.
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)


def annotate_variants(vav: VariantAnnotationVersion,
                      variants_and_columns: dict[str, tuple[Variant, dict]]) -> dict[str, Variant]:
    """ A VariantAnnotation row per name holding those columns, in vav's partition - the table node
        querysets join through (create_fake_variant_annotation writes the base table, which they never read) """
    annotation_run = AnnotationRun.objects.create()
    records = [VariantAnnotation(version=vav, variant=variant, annotation_run=annotation_run, **columns)
               for variant, columns in variants_and_columns.values()]
    bulk_create_in_partition(vav, VariantAnnotation, VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION,
                             records)
    return {name: variant for name, (variant, _) in variants_and_columns.items()}


def create_annotated_variants(vav: VariantAnnotationVersion, rows: dict[str, dict]) -> dict[str, Variant]:
    """ An SNV per name, 100bp apart on chr1, annotated with that name's columns """
    return annotate_variants(vav, {
        name: (slowly_create_test_variant("1", 1000 + i * 100, "A", "T", vav.genome_build), columns)
        for i, (name, columns) in enumerate(rows.items())
    })


def node_keeps(node, variants_by_name: dict[str, Variant], among: Optional[set[str]] = None) -> set[str]:
    """ The names (of among, default all) whose variants pass the node's Q - a node with no Q keeps them all """
    if among is None:
        among = set(variants_by_name)
    pk_names = {variants_by_name[name].pk: name for name in among}
    qs = get_variant_queryset_for_annotation_version(node.analysis.annotation_version).filter(pk__in=pk_names)
    if (q := node._get_node_q()) is not None:
        qs = qs.filter(q)
    return {pk_names[pk] for pk in qs.values_list("pk", flat=True)}
