#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from annotation.models.models import VariantAnnotationVersion
from annotation.models.models_gene_counts import CohortGeneCounts
from annotation.tasks.cohort_sample_gene_damage_counts import get_or_create_gene_count_type_and_values
from snpdb.models import VCF, ImportStatus


def launch_task_for_vcf(gene_count_type, vcf):
    variant_annotation_version = VariantAnnotationVersion.latest(vcf.genome_build)
    cohort = vcf.cohort
    cgc, created = CohortGeneCounts.objects.get_or_create(variant_annotation_version=variant_annotation_version,
                                                          gene_count_type=gene_count_type,
                                                          cohort=cohort,
                                                          cohort_version=cohort.version)
    if created:
        cgc.launch_task()
    else:
        logging.warning("Warning: cohort %s already had CohortGeneCounts for gene count type %s", cohort, gene_count_type)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', default=False)
        parser.add_argument('--clear', action='store_true', default=False)
        parser.add_argument('--gene-count-type', required=True)
        parser.add_argument('--vcf-id', type=int)

    def handle(self, *args, **options):
        gene_count_type_name = options["gene_count_type"]
        vcf_id = options.get("vcf_id")
        all_vcfs = options.get("all")
        clear = options.get("clear")

        if not (all_vcfs or vcf_id):
            msg = "Need to supply at least one of '--vcf-id' or '--all'"
            raise ValueError(msg)

        print(f"GeneCountType: {gene_count_type_name}")
        # Call this now, so that it's done before celery tasks each call it
        gene_count_type, _ = get_or_create_gene_count_type_and_values(gene_count_type_name)

        if clear:
            print("Clearing existing CohortGeneCounts")
            gene_count_type.cohortgenecounts_set.all().delete()
            gene_count_type.genevaluecountcollection_set.all().delete()

        if vcf_id:
            vcf = VCF.objects.get(pk=vcf_id)
            launch_task_for_vcf(gene_count_type, vcf)
        if all_vcfs:
            for vcf in VCF.objects.filter(import_status=ImportStatus.SUCCESS):
                launch_task_for_vcf(gene_count_type, vcf)
