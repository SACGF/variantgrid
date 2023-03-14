#!/usr/bin/env python3

from django.core.management.base import BaseCommand

from genes.models import GeneList, GeneListCategory
from library.guardian_utils import add_public_group_read_permission


class Command(BaseCommand):

    def handle(self, *args, **options):
        for gene_list in GeneList.objects.filter(category__name=GeneListCategory.PANEL_APP_CACHE):
            add_public_group_read_permission(gene_list)
