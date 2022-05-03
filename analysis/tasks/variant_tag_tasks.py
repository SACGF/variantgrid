import celery

from analysis.models import VariantTag, Analysis
from analysis.models.nodes.node_utils import update_analysis
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import ImportSource, VariantAlleleSource, VariantAllele


@celery.shared_task
def variant_tag_created_task(variant_tag_id):
    """ Do this async to save a few miliseconds when adding/removing tags """
    try:
        variant_tag = VariantTag.objects.get(pk=variant_tag_id)
    except VariantTag.DoesNotExist:
        return  # Deleted before this got run, doesn't matter...
    if variant_tag.analysis:
        update_analysis(variant_tag.analysis.pk)
    _liftover_variant_tag(variant_tag)


@celery.shared_task
def variant_tag_deleted_in_analysis_task(analysis_id, _tag_id):
    """ Do this async to save a few miliseconds when adding/removing tags """
    analysis = Analysis.objects.get(pk=analysis_id)
    update_analysis(analysis.pk)


def _liftover_variant_tag(variant_tag: VariantTag):
    populate_clingen_alleles_for_variants(variant_tag.genome_build, [variant_tag.variant])
    variant_allele = VariantAllele.objects.get(variant=variant_tag.variant, genome_build=variant_tag.genome_build)

    # Assign allele to VariantTag
    if variant_tag.allele is None:
        variant_tag.allele = variant_allele.allele
        variant_tag.save()

    allele_source = VariantAlleleSource.objects.create(variant_allele=variant_allele)
    create_liftover_pipelines(admin_bot(), allele_source, ImportSource.WEB, variant_tag.genome_build)
