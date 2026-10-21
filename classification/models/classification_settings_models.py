"""
The classification side of a Lab and an Organization: LabClassificationSettings and OrganizationClassificationSettings,
one row each, reached as `lab.classification_settings` / `organization.classification_settings`. Every Lab and
Organization has one - classification/signals/classification_settings_hooks.py creates it on save. Like Lab, lookups
are cached per request in production (ObjectManagerCachingRequest as the base manager, which the reverse one-to-one uses).
"""
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL

from classification.models.clinvar_key_models import ClinVarKey
from library.django_utils.django_object_managers import ObjectManagerCachingRequest
from snpdb.models.models import Lab, Organization


class OrganizationClassificationSettings(models.Model):
    organization = models.OneToOneField(Organization, primary_key=True, on_delete=CASCADE,
                                        related_name="classification_settings")
    # EvidenceKeyOverrides dict, the lab's classification_config is applied over it
    classification_config = models.JSONField(null=True, blank=True)

    objects = ObjectManagerCachingRequest()

    class Meta:
        base_manager_name = 'objects'
        verbose_name_plural = "organization classification settings"

    def __str__(self):
        return f"{self.organization} classification settings"


class LabClassificationSettings(models.Model):
    lab = models.OneToOneField(Lab, primary_key=True, on_delete=CASCADE, related_name="classification_settings")
    clinvar_key = models.ForeignKey(ClinVarKey, null=True, blank=True, on_delete=SET_NULL,
                                    related_name="lab_classification_settings")
    # EvidenceKeyOverrides dict, applied over the organization's classification_config
    classification_config = models.JSONField(null=True, blank=True)
    # Badged when settings.CLASSIFICATION_DISTINGUISH_RESEARCH, left out of discordance unless settings.DISCORDANCE_RESEARCH_ENABLED
    research = models.BooleanField(default=False, blank=True)
    # One classification per variant rather than per patient, so patient counts including the lab are a lower bound
    consolidates_variant_classifications = models.BooleanField(default=False)
    # Where the lab uploads classification files (may be an S3 directory), enables the classifications/upload page
    upload_location = models.TextField(null=True, blank=True)
    upload_automatic = models.BooleanField(default=False, blank=True)
    upload_instructions = models.TextField(default="", blank=True)

    objects = ObjectManagerCachingRequest()

    class Meta:
        base_manager_name = 'objects'
        verbose_name_plural = "lab classification settings"

    def __str__(self):
        return f"{self.lab} classification settings"
