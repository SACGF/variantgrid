from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.tests.utils import AnalysisSetupMixin
from snpdb.models.models_enums import ImportStatus
from snpdb.models.models_genomic_interval import (
    GenomicIntervalsCategory,
    GenomicIntervalsCollection,
)


class IntersectionNodeBedFileTest(AnalysisSetupMixin, TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        category = GenomicIntervalsCategory.objects.get_or_create(name="Uploaded")[0]
        cls.gic = GenomicIntervalsCollection.objects.create(name="exome", category=category,
                                                            user=User.objects.get_or_create(username="bed_owner")[0],
                                                            genome_build=cls.grch37,
                                                            processed_file="/no/such/dir/exome.processed.bed",
                                                            import_status=ImportStatus.SUCCESS)

    def test_missing_processed_file_is_a_node_configuration_error(self):
        """ eg a database cloned without its data directory - the node says why, the collection is marked """
        node = IntersectionNode(analysis=self.analysis, accordion_panel=IntersectionNode.SELECTED_INTERVALS,
                                genomic_intervals_collection=self.gic)
        errors = node._get_configuration_errors()
        self.assertTrue(any("exome.processed.bed" in e for e in errors), errors)

        self.gic.refresh_from_db()
        self.assertEqual(self.gic.import_status, ImportStatus.ERROR)
