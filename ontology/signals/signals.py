from django.conf import settings

from ontology.tasks.cached_web_resource_tasks import ClinGenCCWebResourceTask
from snpdb.models.models_cached_web_resource import CachedWebResource

gencc_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_GENCC,
                                                                  ClinGenCCWebResourceTask)
