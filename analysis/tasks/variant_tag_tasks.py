import logging
from time import time

import celery
from django.db.models import Q

from analysis.models import VariantTag, Analysis, TagNode
from analysis.models.nodes.node_utils import update_nodes
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import ImportSource, Tag, VariantAlleleSource, VariantAllele


@celery.task
def analysis_tag_created_task(variant_tag_id):
    variant_tag = VariantTag.objects.get(pk=variant_tag_id)
    _update_analysis_on_variant_tag_change(variant_tag.analysis, variant_tag.tag)
    _liftover_variant_tag(variant_tag)


@celery.task
def analysis_tag_deleted_task(analysis_id, tag_id):
    analysis = Analysis.objects.get(pk=analysis_id)
    tag = Tag.objects.get(pk=tag_id)
    _update_analysis_on_variant_tag_change(analysis, tag)


def _update_analysis_on_variant_tag_change(analysis: Analysis, tag: Tag):
    """ TagNodes relying on VariantTag need to be reloaded """

    start = time()
    tag_filter = Q(tag__isnull=True) | Q(tag=tag)
    # TagNode doesn't have any output, so no need to check children
    need_to_reload = False
    for node in TagNode.objects.filter(analysis=analysis).filter(tag_filter):
        node.queryset_dirty = True
        node.save()
        need_to_reload = True

    if need_to_reload:
        logging.info("Reloading %s (%d) as tag %s changed", analysis, analysis.pk, tag)
        update_nodes(analysis.pk)

    end = time()
    time_taken = end - start
    logging.info("Time taken %.3f", time_taken)


def _liftover_variant_tag(variant_tag: VariantTag):
    genome_build = variant_tag.analysis.genome_build
    populate_clingen_alleles_for_variants(genome_build, [variant_tag.variant])
    variant_allele = VariantAllele.objects.get(variant=variant_tag.variant)
    allele_source = VariantAlleleSource.objects.create(variant_allele=variant_allele)
    create_liftover_pipelines(admin_bot(), allele_source, ImportSource.WEB, genome_build)
