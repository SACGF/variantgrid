"""
The VEP pool rule deployment_check warns on (#2477): annotation_workers concurrency x ANNOTATION_VEP_FORK
CPU-bound VEP processes must fit the box's cores.
"""
from django.test import SimpleTestCase

from variantgrid.deployment_validation.vep_check import vep_pool_oversubscription


class VepPoolOversubscriptionTest(SimpleTestCase):

    def test_pool_filling_the_cores_is_fine(self):
        self.assertIsNone(vep_pool_oversubscription(worker_slots=4, fork=1, cpu_count=4))

    def test_fork_multiplies_the_pool(self):
        warning = vep_pool_oversubscription(worker_slots=4, fork=2, cpu_count=4)
        self.assertIn("= 8 VEP processes for 4 cores", warning)

    def test_unset_fork_counts_as_one(self):
        self.assertIsNone(vep_pool_oversubscription(worker_slots=4, fork=0, cpu_count=4))
