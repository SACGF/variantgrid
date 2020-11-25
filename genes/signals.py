from django.conf import settings
import logging

from annotation.models.models import CachedWebResource
from genes.models import GeneListCategory
from genes.tasks.cached_web_resource_tasks import PanelAppEnglandPanelsWebResourceTask, \
    PanelAppAustraliaPanelsWebResourceTask, GnomADGeneConstraintWebResourceTask, PfamWebResourceTask

# For some reason this doesn't work as a variable, has to be stored here...
panel_app_england_panels_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_PANEL_APP_ENGLAND_PANELS,
                                                                                     PanelAppEnglandPanelsWebResourceTask)

panel_app_australia_panels_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_PANEL_APP_AUSTRALIA_PANELS,
                                                                                       PanelAppAustraliaPanelsWebResourceTask)

gnomad_gene_constraint_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_GNOMAD_GENE_CONSTRAINT,
                                                                                   GnomADGeneConstraintWebResourceTask)

pfam_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_PFAM, PfamWebResourceTask)


def cached_third_part_gene_list_pre_delete_handler(sender, instance, **kwargs):
    logging.info("Deleting all associated gene lists.")
    GeneListCategory.objects.filter(company=instance.company).delete()
