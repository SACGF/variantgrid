from django.conf import settings

from annotation.models.models import CachedWebResource
from ontology.tasks.cached_web_resource_tasks import ClinGenCCWebResourceTask

gencc_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_GENCC,
                                                                  ClinGenCCWebResourceTask)
