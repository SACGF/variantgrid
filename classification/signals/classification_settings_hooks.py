"""
Gives every Lab and Organization its classification settings row as it is created, so `lab.classification_settings`
and `organization.classification_settings` always resolve.
"""
from django.db.models.signals import post_save
from django.dispatch import receiver

from classification.models.classification_settings_models import (
    LabClassificationSettings,
    OrganizationClassificationSettings,
)
from snpdb.models.models import Lab, Organization


@receiver(post_save, sender=Organization)
def create_organization_classification_settings(sender, instance: Organization, created: bool, **kwargs):  # pylint: disable=unused-argument
    if created:
        OrganizationClassificationSettings.objects.get_or_create(organization=instance)


@receiver(post_save, sender=Lab)
def create_lab_classification_settings(sender, instance: Lab, created: bool, **kwargs):  # pylint: disable=unused-argument
    if created:
        LabClassificationSettings.objects.get_or_create(lab=instance)
