"""
Manual variant entry: a user's search / "create variant" text becomes a ManualVariantEntryCollection of one
ManualVariantEntry per line, and CreatedManualVariant records the Variants made from each.
annotation/manual_variant_entry.py does the parsing and creation.
"""

from functools import cached_property
from typing import Optional

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models.deletion import CASCADE

from annotation.enums import ManualVariantEntryType
from annotation.models.models_annotation_run import AnnotationRun
from snpdb.enums import ImportStatus
from snpdb.models import (
    DBSNP_PATTERN,
    HGVS_UNCLEANED_PATTERN,
    VARIANT_PATTERN,
    VARIANT_SYMBOLIC_PATTERN,
    GenomeBuild,
    Variant,
)


class ManualVariantEntryCollection(models.Model):
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)
    created = models.DateTimeField(auto_now_add=True)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    celery_task = models.CharField(max_length=36, null=True)

    @cached_property
    def first_entry(self) -> Optional['ManualVariantEntry']:
        """ Often a user will enter a single variant into search and then wait for it to be created/annotated
            We have the "first" one to deal with this """
        return self.manualvariantentry_set.all().order_by("pk").first()

    @cached_property
    def first_variant(self) -> Optional['Variant']:
        variant = None
        if mve := self.first_entry:
            if cve := mve.unique_created_variants.filter(variant__isnull=False).order_by("variant_id").first():
                variant = cve.variant
        return variant

    def first_variant_annotation_run(self) -> Optional['AnnotationRun']:
        annotation_run = None
        if mve := self.first_entry:
            if variant := self.first_variant:
                annotation_run = AnnotationRun.get_for_variant(variant, mve.genome_build)
        return annotation_run

    @classmethod
    def filter_for_user(cls, user):
        """ Superusers see all collections, other users only their own """
        qs = cls.objects.all()
        if not user.is_superuser:
            qs = qs.filter(user=user)
        return qs

    @classmethod
    def get_for_user(cls, user, mvec_id):
        try:
            return cls.filter_for_user(user).get(pk=mvec_id)
        except cls.DoesNotExist:
            msg = f"User {user} doesn't have permission to access mvec_id={mvec_id}"
            raise PermissionDenied(msg)

    def __str__(self):
        status = self.get_import_status_display()
        if self.import_status == ImportStatus.IMPORTING:
            status += f", task: {self.celery_task}"

        count = self.manualvariantentry_set.count()
        return f"{self.pk}: {count} manual entries ({status})"


class ManualVariantEntry(models.Model):
    manual_variant_entry_collection = models.ForeignKey(ManualVariantEntryCollection, on_delete=CASCADE)
    line_number = models.IntegerField()
    entry_text = models.TextField()
    entry_type = models.CharField(max_length=1, choices=ManualVariantEntryType.choices, default=ManualVariantEntryType.UNKNOWN)
    warning_message = models.TextField(blank=True, null=True)
    error_message = models.TextField(blank=True, null=True)  # Set if any error...

    @property
    def unique_created_variants(self):
        return self.createdmanualvariant_set.distinct('variant_id').all()

    @property
    def genome_build(self) -> GenomeBuild:
        return self.manual_variant_entry_collection.genome_build

    @staticmethod
    def get_entry_type(line: str) -> ManualVariantEntryType:
        MATCHERS = {
            DBSNP_PATTERN: ManualVariantEntryType.DBSNP,
            HGVS_UNCLEANED_PATTERN: ManualVariantEntryType.HGVS,
            VARIANT_PATTERN: ManualVariantEntryType.VARIANT,
            VARIANT_SYMBOLIC_PATTERN: ManualVariantEntryType.VARIANT,
        }
        for k, v in MATCHERS.items():
            if k.match(line):
                return v
        return ManualVariantEntryType.UNKNOWN


class CreatedManualVariant(models.Model):
    manual_variant_entry = models.ForeignKey(ManualVariantEntry, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
