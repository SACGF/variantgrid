"""
CachedWebResource: a named data set that populates itself from the web (HGNC, MANE, PanelApp, GenCC,
ClinVar citations etc). Creating one with a name from settings (CACHED_WEB_RESOURCE_*) fires a post_save
handler built by CachedWebResource.named_handler_factory, which each app connects in its apps.py to start
its loader task (a snpdb/tasks/cached_web_resource_tasks.py:CachedWebResourceTask subclass).
"""
from django.db import models
from django_extensions.db.models import TimeStampedModel

from snpdb.models.models_enums import ImportStatus

__all__ = ["CachedWebResource"]


class CachedWebResource(TimeStampedModel):
    """ These are annotations that can populate themselves via the web

        * Can be manually updated via button kicking off annotation.views.load_cached_web_resource
        * List of names each system uses is in settings
        * Creation of one with name triggers signals which apps can configure to listen to / run populators
    """
    name = models.TextField(primary_key=True)
    description = models.TextField(blank=True)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)

    class Meta(TimeStampedModel.Meta):
        # Moved from annotation (#2278) as a state-only migration, keeping its table
        db_table = "annotation_cachedwebresource"

    @staticmethod
    def named_handler_factory(name, celery_task_class):
        def handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
            created = kwargs.get("created")

            if created and instance.name == name:
                task_signature = celery_task_class.si(instance.pk)
                task_signature.apply_async()

        return handler
