"""
    https://github.com/SACGF/variantgrid/issues/1601

    Need to trigger reloads of bad metrics, so we die properly

"""

import logging

from django.core.management.base import BaseCommand

from seqauto.models import IlluminaFlowcellQC
from snpdb.models import DataState


class Command(BaseCommand):

    def handle(self, *args, **options):

        qs = IlluminaFlowcellQC.objects.exclude(data_state=DataState.ERROR)
        qs = qs.filter(mean_cluster_density__isnull=True)

        if not qs.exists():
            logging.info("No potentially bad IlluminaFlowcellQC records")

        for iqc in qs:
            logging.info(f"Reloading: {iqc}")
            iqc.load_from_file(None)

            logging.info(f"{iqc}: {iqc.get_data_state_display()}")
