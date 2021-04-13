import celery

from analysis.models import VariantTag, Analysis
from analysis.models.nodes.node_utils import update_nodes
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import ImportSource, VariantAlleleSource, VariantAllele


@celery.task
def analysis_tag_created_task(variant_tag_id):
    """ Do this async to save a few miliseconds when adding/removing tags """
    try:
        variant_tag = VariantTag.objects.get(pk=variant_tag_id)
    except VariantTag.DoesNotExist:
        return  # Deleted before this got run, doesn't matter...
    update_nodes(variant_tag.analysis.pk)
    _liftover_variant_tag(variant_tag)


@celery.task
def analysis_tag_deleted_task(analysis_id, _tag_id):
    """ Do this async to save a few miliseconds when adding/removing tags """
    analysis = Analysis.objects.get(pk=analysis_id)
    update_nodes(analysis.pk)


def _liftover_variant_tag(variant_tag: VariantTag):
    genome_build = variant_tag.analysis.genome_build
    populate_clingen_alleles_for_variants(genome_build, [variant_tag.variant])
    variant_allele = VariantAllele.objects.get(variant=variant_tag.variant)
    allele_source = VariantAlleleSource.objects.create(variant_allele=variant_allele)
    create_liftover_pipelines(admin_bot(), allele_source, ImportSource.WEB, genome_build)
