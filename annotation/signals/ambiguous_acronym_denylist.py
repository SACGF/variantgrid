""" A new OntologyVersion is a new denylist cache key - build it in the background before a page asks for it """
from django.db import transaction
from django.db.models.signals import post_save
from django.dispatch import receiver

from annotation.tasks.ambiguous_acronym_denylist_task import build_ambiguous_acronym_denylist_task
from ontology.models import OntologyVersion


@receiver(post_save, sender=OntologyVersion)
def ontology_version_post_save_handler(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    if created:
        transaction.on_commit(build_ambiguous_acronym_denylist_task.delay)
