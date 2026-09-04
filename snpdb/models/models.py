"""

This is the main app for VariantGrid, and this file contains the main models for Variants / Samples
etc and things that don't fit anywhere else.

'snpdb' was the highly unoriginal name I used before 'VariantGrid'
"""
import json
import logging
import os
import re
import uuid
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property, total_ordering
from html import escape
from re import RegexFlag
from typing import Optional, TypedDict

from celery import signature
from celery.result import AsyncResult
from django.conf import settings
from django.contrib.auth.models import Group, User
from django.core.cache import cache
from django.core.exceptions import FieldDoesNotExist, PermissionDenied, ValidationError
from django.db import models, transaction
from django.db.models import QuerySet, TextChoices
from django.db.models.deletion import CASCADE, PROTECT, SET_NULL
from django.db.models.signals import pre_delete
from django.dispatch import receiver
from django.urls import reverse
from django.utils import timezone
from django.utils.safestring import SafeString
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager
from more_itertools import first

from classification.enums.classification_enums import AlleleOriginBucket, ShareLevel
from eventlog.models import create_event
from library.django_utils import get_url_from_media_root_filename
from library.django_utils.django_object_managers import ObjectManagerCachingRequest
from library.enums.log_level import LogLevel
from library.guardian_utils import admin_bot
from library.log_utils import AdminNotificationBuilder
from library.preview_request import PreviewModelMixin
from library.utils import JsonObjType, import_class
from snpdb.models.models_enums import AwardPeriod, UserAwardKind, UserAwardLevel
from snpdb.user_awards import get_award_definition
from user_messages.models import Message


# A tag used for both origins is stored as UNKNOWN so AlleleOriginFilterDefault.buckets works verbatim.
# On a tag that means "applies to both", which is what people see
TAG_ALLELE_ORIGIN_CHOICES = [
    (AlleleOriginBucket.UNKNOWN, "Both"),
    (AlleleOriginBucket.GERMLINE, "Germline only"),
    (AlleleOriginBucket.SOMATIC, "Somatic only"),
]


class Tag(models.Model):
    id = models.CharField(max_length=50, primary_key=True)
    # Retired tags are kept rather than deleted so old names still resolve, and so merges leave a trail.
    # They're hidden from anywhere you pick a tag, but rows already pointing at them keep working
    retired = models.DateTimeField(null=True, blank=True)
    merged_into = models.ForeignKey('self', null=True, blank=True, on_delete=SET_NULL)
    allele_origin_bucket = models.CharField(max_length=1, choices=TAG_ALLELE_ORIGIN_CHOICES,
                                            default=AlleleOriginBucket.UNKNOWN)

    @classmethod
    def live_qs(cls) -> QuerySet['Tag']:
        return cls.objects.filter(retired__isnull=True)

    @property
    def active(self) -> bool:
        return self.retired is None

    def __str__(self):
        return self.id


class CachedGeneratedFile(TimeStampedModel):
    # We use a UUID for these so people can't guess the IDs and see other people's graphs/downloads etc
    # The generator and params hash are unique - we look these up in a view where we do the appropriate
    # security checks on any objects - then can just pass the hash around (big enough to be secret)
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    filename = models.TextField(null=True)
    exception = models.TextField(null=True)
    generator = models.TextField()
    params_hash = models.TextField()  # sha256 of params
    task_id = models.CharField(max_length=36, null=True)
    task_status = models.TextField(null=True)
    generate_start = models.DateTimeField(null=True)
    generate_end = models.DateTimeField(null=True)
    progress = models.FloatField(null=True)

    class Meta:
        unique_together = ("generator", "params_hash")

    def __str__(self):
        if self.filename:
            description = f"file: {self.filename} made: {self.generate_end}"
        else:
            description = f"task: {self.task_id} sent: {self.generate_start}"
        return f"{self.generator}({self.params_hash}): {description}"

    def get_absolute_url(self):
        return reverse("cached_generated_file_check", kwargs={"cgf_id": self.pk})

    def get_media_url(self):
        if self.filename is None:
            raise ValueError(f"{self}.filename is None")
        return get_url_from_media_root_filename(self.filename)

    @property
    def file_missing(self) -> bool:
        """ We generated a file but it's not there anymore - eg media_root cleaned up, or the database
            was copied from a deployment with a different MEDIA_ROOT """
        if not self.filename:
            return False  # Still generating - nothing promised yet
        if not self.filename.startswith(os.path.join(settings.MEDIA_ROOT, "")):
            return True
        return not os.path.exists(self.filename)

    @staticmethod
    def get_or_create_and_launch(generator, params_hash, task: signature) -> 'CachedGeneratedFile':
        cgf, created = CachedGeneratedFile.objects.get_or_create(generator=generator,
                                                                 params_hash=params_hash)
        if cgf.file_missing:
            # Drop the row so the get_or_create below regenerates it lazily
            logging.info("Discarding CachedGeneratedFile %s - %s is gone", cgf.pk, cgf.filename)
            cgf.delete()
            cgf, created = CachedGeneratedFile.objects.get_or_create(generator=generator,
                                                                     params_hash=params_hash)
        if created or not cgf.task_id:
            logging.debug("Launching Celery Job for CachedGeneratedFile(generator=%s, params_hash=%s)",
                          generator, params_hash)
            async_result = task.apply_async()
            cgf.task_id = async_result.id
            cgf.generate_start = timezone.now()
            cgf.progress = 0.0
            cgf.save()
        else:
            async_result = AsyncResult(cgf.task_id)

        if async_result.result:
            cgf.save_from_async_result(async_result)

        return cgf

    def save_from_async_result(self, async_result):
        self.task_status = async_result.status
        self.generate_end = timezone.now()
        if async_result.successful():
            self.filename = async_result.result
        else:
            self.exception = async_result.result
        self.save()

    def check_ready(self):
        if not self.generate_end:
            async_result = AsyncResult(self.task_id)
            if async_result.result:
                self.save_from_async_result(async_result)

    def save(self, *args, **kwargs):
        if settings.CACHE_GENERATED_FILES:
            logging.debug("Save cachedgeneratedfile!")
            super().save(*args, **kwargs)
        else:
            logging.debug("Not caching generated files - skipped save!")


@receiver(pre_delete, sender=CachedGeneratedFile)
def cgf_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    if instance.filename:
        if os.path.exists(instance.filename):
            os.unlink(instance.filename)

class Company(models.Model):
    name = models.TextField(primary_key=True)

    @staticmethod
    def get_our_company():
        from genes.models import GeneListCategory

        company = None
        company_name = getattr(settings, "COMPANY", None)
        if company_name:
            company, created = Company.objects.get_or_create(name=company_name)
            if created:
                pathology_test_category = GeneListCategory.get_or_create_category(GeneListCategory.PATHOLOGY_TEST)
                pathology_test_category.company = company
                pathology_test_category.save()

        return company

    def __str__(self):
        return self.name


class Manufacturer(models.Model):
    name = models.TextField()

    def __str__(self):
        return self.name


class SoftwareVersion(models.Model):
    name = models.TextField()
    version = models.TextField()

    def __str__(self):
        return f"{self.name} ({self.version})"

    class Meta:
        abstract = True


class Wiki(TimeStampedModel):
    # Inherit and add a 1-to-1 field for whatever object you want
    objects = InheritanceManager()
    last_edited_by = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    markdown = models.TextField()

    def _get_restricted_object(self):
        return None

    def can_write(self, user) -> bool:
        obj = self._get_restricted_object()
        if obj and not obj.can_write(user):
            return False
        return True

    def check_user_edit_permission(self, user):
        """ Overwrite to restrict users """
        if not self.can_write(user):
            msg = f"You do not have permission to write to wiki {self.pk}"
            raise PermissionDenied(msg)

    @staticmethod
    def get_subclass_by_name(class_name):
        klass = import_class(class_name)
        if not issubclass(klass, Wiki):
            msg = f"Wiki.get_or_create passed class_name='{class_name}' which is non-subclass of Wiki"
            raise PermissionDenied(msg)
        return klass

    @staticmethod
    def get_or_create(class_name, unique_keyword, unique_value, user=None):
        klass = Wiki.get_subclass_by_name(class_name)
        try:
            klass._meta.get_field(unique_keyword)
        except FieldDoesNotExist:
            raise PermissionDenied(f"'{unique_keyword}' is not a valid field for {class_name}")
        if user is not None:
            # Check permission on an unsaved instance before creating the row (throws 403)
            klass(**{unique_keyword: unique_value}).check_user_edit_permission(user)
        wiki, _ = klass.objects.get_or_create(**{unique_keyword: unique_value})
        return wiki

    def __str__(self):
        return self.markdown


class ImportedWikiCollection(models.Model):
    """ VCF imports only support 1 file at a time, so we can't import a multi build file easily """
    match_column_name = models.TextField()  # e.g. gene/variant
    genome_build = models.ForeignKey('GenomeBuild', null=True, on_delete=CASCADE)


class ImportedWiki(models.Model):
    collection = models.ForeignKey(ImportedWikiCollection, on_delete=CASCADE)
    match_column_value = models.TextField()  # From collection.match_column_name
    markdown = models.TextField()
    matched_wiki = models.ForeignKey(Wiki, null=True, on_delete=SET_NULL)
    username = models.TextField()
    created = models.DateTimeField()  # Time on original server
    modified = models.DateTimeField()  # Time on original server


class Organization(models.Model, PreviewModelMixin):
    # If you add fields @see OrganizationAdmin
    name = models.TextField()
    short_name = models.TextField(blank=False, null=True)  # Don't use for anything other than human labels
    group_name = models.TextField(blank=True, null=True, unique=True)
    classification_config = models.JSONField(null=True, blank=True)
    active = models.BooleanField(default=True, blank=True)

    objects = ObjectManagerCachingRequest()

    class Meta:
        ordering = ['name']
        verbose_name = 'Organisation'
        verbose_name_plural = 'Organisations'
        default_manager_name = 'objects'

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-building"

    def __lt__(self, other):
        return self.name < other.name

    def get_absolute_url(self):
        return reverse('view_organization', kwargs={"organization_id": self.pk})

    def __str__(self):
        return self.name

    @property
    def shortest_name(self) -> str:
        return self.short_name or self.name

    @cached_property
    def classifying_labs(self) -> list['Lab']:
        return [lab for lab in self.lab_set.all().order_by('name') if lab.total_classifications > 0]

    @cached_property
    def sharing_labs(self) -> list['Lab']:
        return [lab for lab in self.lab_set.all().order_by('name') if lab.total_shared_classifications > 0]

    @property
    def group(self) -> Optional[Group]:
        if self.group_name:
            group, _ = Group.objects.get_or_create(name=self.group_name)
            return group
        return None

    def is_member(self, user: User) -> bool:
        return Lab.valid_labs_qs(user).filter(organization=self).exists()

    def add_member(self, user: User):
        """ Org group membership only - it grants no lab. A user already in the group is a no-op,
            so this is safe to call on every login """
        group = self.group
        if group is None or group.user_set.filter(pk=user.pk).exists():
            return

        group.user_set.add(user)
        org_message = getattr(settings, "USER_CREATE_ORG_MESSAGE", {})
        if message := org_message.get(self.group_name):
            Message.objects.create(subject="Your initial lab / organisation",
                                   body=message,
                                   sender=admin_bot(),
                                   recipient=user)

    def can_write(self, user: User) -> bool:
        return user.is_superuser or self.lab_set.filter(labhead__user=user).exists()

    def check_can_write(self, user):
        if not self.can_write(user):
            msg = f"You do not have WRITE permission for {self.pk}"
            raise PermissionDenied(msg)


class LabUserRole(TextChoices):
    HEAD = "head", "Lab Head"
    USER = "user", "User"


@total_ordering
class LabUser:

    def __init__(self, user: User, role: LabUserRole):
        self.user = user
        self.role = role

    @property
    def is_head(self) -> bool:
        return self.role == LabUserRole.HEAD

    @cached_property
    def preferred_label(self) -> str:
        from snpdb.models import AvatarDetails
        return AvatarDetails.avatar_for(self.user).preferred_label

    def __lt__(self, other):
        return self.preferred_label < other.preferred_label

    def __eq__(self, other):
        return self.user == other.user


class ClinVarAssertionMethods(TextChoices):
    # "yes", "no", "unknown", "not provided", "not applicable"
    yes = "yes"
    no = "no"
    unknown = "unknown"
    not_provided = "not provided"
    not_applicable = "not applicable"


class ClinVarCitationsModes(TextChoices):
    all = "all"
    interpretation_summary_only = "interpret"


class ClinVarKey(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar key"

    id = models.TextField(primary_key=True)
    api_key = models.TextField(null=True, blank=True)
    org_id = models.TextField(null=False, blank=True, default='')  # maybe this should have been the id?
    name = models.TextField(blank=True, default='')
    last_full_run = models.DateTimeField(null=True)

    default_affected_status = models.TextField(choices=ClinVarAssertionMethods.choices, null=True, blank=True)
    inject_acmg_description = models.BooleanField(blank=True, default=False)
    include_interpretation_summary = models.BooleanField(blank=True, default=True)
    assertion_method_lookup = models.JSONField(null=False, default=dict)
    citations_mode = models.TextField(choices=ClinVarCitationsModes.choices, default=ClinVarCitationsModes.all)

    def assertion_criteria_vg_to_code(self, vg_value: str) -> Optional[JsonObjType]:
        """
        :param vg_value: Value as stored in the evidence key
        :return: The code we map to, None indicates we don't have a mapping
        """
        def map_value(raw_value: str):
            if not raw_value:
                raw_value = ""

            if self.assertion_method_lookup:
                if lookups := self.assertion_method_lookup.get("lookups"):
                    for lookup in lookups:
                        if match_text := lookup.get("match"):
                            if re.compile(match_text, RegexFlag.IGNORECASE).match(raw_value):
                                return lookup.get("citation")
                        else:
                            raise ValueError("Assertion Method Lookup, lacking a 'match' field")
                else:
                    raise ValueError("Assertion Method Lookup, needs a key 'lookups', with dict entries of 'match' and 'citation'")
            elif raw_value == "acmg":
                # if no method lookups have been setup, only accept "acmg"
                return raw_value
            else:
                raise ValueError(f"Unable to map {raw_value=}")

        mapped_value = map_value(vg_value)
        if mapped_value == "acmg":
            return {
                "db": "PubMed",
                "id": "PMID:25741868"
            }
        return mapped_value

    @property
    def label(self) -> str:
        return self.name if self.name else self.id

    def __str__(self):
        return f"ClinVarKey ({self.label})"

    def __lt__(self, other):
        return self.id < other.id

    @staticmethod
    def clinvar_keys_for_user(user: User) -> QuerySet:
        """
        Ideally this would be on ClinVarKey but can't be due to ordering
        """
        if user.is_superuser:
            return ClinVarKey.objects.all().order_by('pk')

        labs = Lab.valid_labs_qs(user).filter(clinvar_key__isnull=False)
        if labs:
            return ClinVarKey.objects.filter(pk__in=labs.values_list('clinvar_key', flat=True)).order_by('pk')
        else:
            return ClinVarKey.objects.none()

    def check_user_can_access(self, user: User) -> None:
        """
        :throw PermissionDenied if user is not associated to ClinVarKey
        """
        if not user.is_superuser:
            allowed_clinvar_keys = ClinVarKey.clinvar_keys_for_user(user)
            if not allowed_clinvar_keys.filter(pk=self.pk).exists():
                raise PermissionDenied("User does not belong to a lab that uses the submission key")


class ClinVarKeyExcludePatternMode(TextChoices):
    EXCLUDE_IF_MATCH = 'X'
    EXCLUDE_IF_NO_MATCH = 'N'


class ClinVarKeyExcludePattern(TimeStampedModel):
    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=CASCADE)
    evidence_key = models.TextField()  # Actually EvidenceKey but can't link to that from snpdb
    pattern = models.TextField()
    name = models.TextField(blank=True)
    case_insensitive = models.BooleanField(default=True, blank=True)
    mode = models.TextField(choices=ClinVarKeyExcludePatternMode.choices, default=ClinVarKeyExcludePatternMode.EXCLUDE_IF_MATCH)

    @cached_property
    def regex(self):
        flags = 0
        if self.case_insensitive:
            flags = RegexFlag.IGNORECASE
        return re.compile(self.pattern, flags=flags)

    def should_exclude(self, text: str) -> bool:
        matches = bool(self.regex.search(text))
        if self.mode == ClinVarKeyExcludePatternMode.EXCLUDE_IF_MATCH:
            return matches
        else:
            return not matches

    def clean(self):
        try:
            re.compile(self.pattern)
        except Exception as exc:
            raise ValidationError({'pattern': ValidationError(f'{self.pattern} is not a valid regular expression')}) from exc

        from classification.models import EvidenceKeyMap
        if EvidenceKeyMap.cached_key(self.evidence_key).is_dummy:
            raise ValidationError({'evidence_key': ValidationError(f'{self.evidence_key} is not a valid EvidenceKey')})

    def __str__(self):
        from classification.models import EvidenceKeyMap
        return EvidenceKeyMap.cached_key(self.evidence_key).pretty_label + " : " + (self.name or self.pattern)


class Country(models.Model):
    name = models.TextField(primary_key=True)
    short_name = models.TextField(unique=True, null=True)
    population = models.IntegerField(null=True)

    def __str__(self):
        return self.name


class State(models.Model):
    name = models.TextField(primary_key=True)
    short_name = models.TextField(unique=True, null=True)
    country = models.ForeignKey(Country, on_delete=CASCADE)
    population = models.IntegerField(null=True)

    def __str__(self):
        return self.name


@dataclass
class ContactDetails:
    website: Optional[str] = ""
    name: Optional[str] = ""
    phone: Optional[str] = ""
    email: Optional[str] = ""

    def __bool__(self):
        # contact name by its self isn't useful
        return bool(self.website) or bool(self.phone) or bool(self.email)


class Lab(models.Model, PreviewModelMixin):
    name = models.TextField()
    external = models.BooleanField(default=False, blank=True)  # From somewhere else, e.g. Shariant
    research = models.BooleanField(default=False, blank=True)
    city = models.TextField()
    state = models.ForeignKey(State, null=True, on_delete=PROTECT)
    country = models.ForeignKey(Country, null=True, on_delete=PROTECT)
    url = models.TextField(blank=True)
    css_class = models.TextField(blank=True)
    lat = models.FloatField(null=True, blank=True)
    long = models.FloatField(null=True, blank=True)

    contact_name = models.TextField(blank=True)
    contact_phone = models.TextField(blank=True)
    contact_email = models.TextField(blank=True)

    # Opt-in to participate in MatchMaker Exchange (patient matchmaking). Distinct from
    # a record's share_level=PUBLIC (which only means "publicly shareable"); enabling MME
    # is a deliberate lab decision and requires a resolvable MME contact (see clean()).
    mme_enabled = models.BooleanField(default=False, blank=True)

    clinvar_key = models.ForeignKey(ClinVarKey, null=True, blank=True, on_delete=SET_NULL)

    group_name = models.TextField(blank=True, null=True, unique=True)
    classification_config = models.JSONField(null=True, blank=True)

    # want every lab to have an organization, but not going to have them
    # at point of migration
    organization = models.ForeignKey(Organization, null=False, blank=False, on_delete=CASCADE)
    # location where the lab can upload files to, (in some environments may refer to s3 directory)
    upload_location = models.TextField(null=True, blank=True)
    upload_automatic = models.BooleanField(default=False, blank=True)
    upload_instructions = models.TextField(default="", blank=True)

    consolidates_variant_classifications = models.BooleanField(default=False)
    # does this lab produce one classification

    """
    If provided, and filename matches, file upload will be automatically set to auto_processed
    """

    email = models.TextField(blank=True)
    slack_webhook = models.TextField(blank=True)

    def __lt__(self, other):
        if self.organization != other.organization:
            return self.organization < other.organization
        return self.name < other.name

    objects = ObjectManagerCachingRequest()

    class Meta:
        ordering = ['name']
        base_manager_name = 'objects'

    @classmethod
    def preview_category(cls) -> str:
        return "Organisation Lab"

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-flask"

    def clean(self):
        super().clean()
        # An MME match is an invitation to a person-to-person clinical conversation, so a
        # lab can't opt in without an address a matching lab can reply to.
        if self.mme_enabled and not (self.contact_email or "").strip():
            raise ValidationError({
                "mme_enabled": "Set an MME contact email before enabling MatchMaker Exchange."
            })

    @property
    def contact_details(self) -> ContactDetails:
        return ContactDetails(
            name=self.contact_name,
            phone=self.contact_phone,
            email=self.contact_email,
            website=self.url
        )

    @property
    def group(self) -> Optional[Group]:
        if self.group_name:
            group, _ = Group.objects.get_or_create(name=self.group_name)
            return group
        return None

    @property
    def group_institution(self) -> Optional[Group]:
        if self.group_name:
            parts = self.group_name.split('/')
            if len(parts) >= 2:
                inst_group_name = '/'.join(parts[:-1])
                group, _ = Group.objects.get_or_create(name=inst_group_name)
                return group

    @cached_property
    def active_users(self) -> QuerySet[User]:
        return self.group.user_set.filter(is_active=True)

    @cached_property
    def lab_users(self) -> list[LabUser]:
        users = list(self.group.user_set.filter(is_active=True))
        heads = set(self.labhead_set.values_list('user_id', flat=True))
        lab_users: list[LabUser] = []
        for user in users:
            role = LabUserRole.USER
            if user.id in heads:
                role = LabUserRole.HEAD
            lab_users.append(LabUser(user=user, role=role))
        lab_users.sort()
        return lab_users

    @cached_property
    def shared_classifications(self) -> QuerySet['Classification']:
        return self.classification_set.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).exclude(withdrawn=True)

    @cached_property
    def classifications(self) -> QuerySet['Classification']:
        """ Shared or all classifications based on settings.CLASSIFICATION_STATS_USE_SHARED """
        if settings.CLASSIFICATION_STATS_USE_SHARED:
            qs = self.shared_classifications
        else:
            qs = self.classification_set.all()
        return qs

    @cached_property
    def total_classifications(self) -> int:
        """ Total 'classifications' as per above method """
        return self.classifications.count()

    @cached_property
    def total_shared_classifications(self) -> int:
        return self.shared_classifications.count()

    @cached_property
    def total_unshared_classifications(self) -> int:
        all_classifications = self.classification_set.exclude(withdrawn=True).count()
        return all_classifications - self.total_shared_classifications

    @cached_property
    def classifications_by_created(self) -> QuerySet['Classification']:
        return self.classifications.order_by("created")

    @cached_property
    def classifications_by_modified(self) -> QuerySet['Classification']:
        return self.classifications.order_by("-modified")

    @property
    def first_classification_ever_shared_date(self) -> Optional[datetime]:
        return self.classification_set.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).values_list('created', flat=True).order_by('created').first()

    @cached_property
    def classifications_per_day(self) -> float:
        try:
            latest = now()
            days = (latest - self.first_classification_ever_shared_date).days
            if days == 0:
                days = 1
            cpd = self.total_classifications / days
        except Exception:
            cpd = 0
        return cpd

    @cached_property
    def classifications_per_week(self) -> float:
        return self.classifications_per_day * 7

    @staticmethod
    def valid_labs_qs(user: User, admin_check=False) -> QuerySet['Lab']:
        # as organization is used for sorting, it's generally always a good idea to select related it
        if admin_check and user.is_superuser:
            return Lab.objects.select_related('organization')

        group_names = list(user.groups.values_list('name', flat=True))
        return Lab.objects.select_related('organization').filter(group_name__in=group_names).order_by('name')

    @staticmethod
    def has_active_lab(user: User, admin_check=False) -> bool:
        """ Whether the user has access to any lab belonging to an active organization """
        return Lab.valid_labs_qs(user=user, admin_check=admin_check).filter(organization__active=True).exists()

    """
    # these methods have been superseeded by having full classification activity by lab
    def classifications_activity(self, time_period: timedelta):
        trunc_func = TimePeriod.truncate_func(time_period)
        qs = self.classifications.annotate(time_period=trunc_func("created")).values("time_period") \
            .annotate(num_classifications=Count("id")).order_by("-time_period")
        return qs.values_list("time_period", "num_classifications")

    def classifications_activity_by_day(self):
        return self.classifications_activity(TimePeriod.DAY)
    """

    def is_member(self, user: User, admin_check=False) -> bool:
        return self.valid_labs_qs(user=user, admin_check=admin_check).filter(pk=self.pk).exists()

    def is_head(self, user: User) -> bool:
        return self.labhead_set.filter(user=user).exists()

    def can_write(self, user: User) -> bool:
        return user.is_superuser or self.is_head(user)

    def check_can_write(self, user):
        if not self.can_write(user):
            msg = f"You do not have WRITE permission for {self.pk}"
            raise PermissionDenied(msg)

    def can_manage_members(self, user: User) -> bool:
        return settings.LAB_HEAD_MANAGE_MEMBERS and (user.is_superuser or self.is_head(user))

    def check_can_manage_members(self, user: User):
        if not self.can_manage_members(user):
            msg = f"You do not have permission to manage members of {self.pk}"
            raise PermissionDenied(msg)

    @property
    def external_membership(self) -> Optional[str]:
        """ Why membership of this lab is decided outside VariantGrid (eg an Active Directory group) """
        return settings.LAB_EXTERNAL_MEMBERSHIP.get(self.group_name)

    def candidate_members_qs(self, for_user: User) -> QuerySet[User]:
        """ Active users that for_user may add to this lab - their organisation, or everyone for a superuser """
        group = self.group
        if group is None:
            return User.objects.none()

        qs = User.objects.filter(is_active=True)
        if not for_user.is_superuser:
            group_institution = self.group_institution
            if group_institution is None:
                return User.objects.none()
            qs = qs.filter(groups=group_institution)
        return qs.exclude(groups=group).order_by('last_name', 'first_name', 'username')

    def _clear_member_caches(self):
        for cached in ("active_users", "lab_users"):
            self.__dict__.pop(cached, None)

    @transaction.atomic
    def add_member(self, user: User, added_by: User):
        """ Lab group membership, plus the organisation group that goes with it. Already a member is a no-op """
        group = self.group
        if group is None or group.user_set.filter(pk=user.pk).exists():
            return

        group.user_set.add(user)
        if group_institution := self.group_institution:
            group_institution.user_set.add(user)

        from snpdb.models.models_user_settings import UserSettingsOverride  # Circular import
        user_settings_override, _ = UserSettingsOverride.objects.get_or_create(user=user)
        user_settings_override.auto_set_default_lab()
        user_settings_override.save()

        self._clear_member_caches()
        self._log_membership_change(user, added_by, added=True)

        Message.objects.create(subject=f"You have been added to {self}",
                               body=f"You are now a member of the lab [{self}]({self.get_absolute_url()}).",
                               sender=admin_bot(),
                               recipient=user)

    def remove_member_blocked_reason(self, user: User, removed_by: User) -> Optional[str]:
        """ Why removed_by may not remove user from this lab, or None if they may. Superusers bypass """
        if removed_by.is_superuser:
            return None
        if reason := self.external_membership:
            return f"Membership of {self} is managed outside {settings.SITE_NAME}. {reason}"
        if user == removed_by:
            return "You cannot remove yourself from a lab"
        if self.is_head(user) and self.labhead_set.count() == 1:
            return f"{user} is the only lab head of {self}"
        return None

    @transaction.atomic
    def remove_member(self, user: User, removed_by: User):
        """ Removes the lab group only - organisation membership is the user's standing in the org, not the lab """
        group = self.group
        if group is None:
            return

        if reason := self.remove_member_blocked_reason(user, removed_by):
            raise PermissionDenied(reason)

        group.user_set.remove(user)
        self.labhead_set.filter(user=user).delete()

        self._clear_member_caches()
        self._log_membership_change(user, removed_by, added=False)

    def _log_membership_change(self, user: User, changed_by: User, added: bool):
        """ Membership is a permission change, so it's logged for the lab and notified to admins """
        verb = "added to" if added else "removed from"
        details = f"{user} ({user.username}) {verb} lab {self} by {changed_by} ({changed_by.username})"
        event_name = "lab_member_added" if added else "lab_member_removed"
        create_event(user=changed_by, name=event_name, details=details, app_name="snpdb")

        nb = AdminNotificationBuilder(message="Lab Membership Changed")
        nb.add_markdown(details)
        nb.add_field("Lab", str(self))
        nb.add_field("User", user.username)
        nb.add_field("Changed by", changed_by.username)
        nb.send()

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        if self.group_name:
            # pre-create the groups
            _ = self.group
            _ = self.group_institution
            # TODO assign the lab permissions to the groups
            # but also make sure this doesn't break RunX1

    def get_absolute_url(self):
        return reverse('view_lab', kwargs={"lab_id": self.pk})

    def __str__(self):
        return f"{self.organization.shortest_name} / {self.name}"


class LabHead(models.Model):
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)

    class Meta:
        unique_together = ('lab', 'user')

    def __str__(self):
        return f"{self.lab}: {self.user}"


class UserAward(TimeStampedModel):
    """ Three kinds (see UserAwardKind). Computed rows (titles/badges) are keyed on
        (definition_key, subject, period) and reused: 'created' is first earned, 'modified' the last
        change of holder/tier. See snpdb.user_awards for the definitions and the computation """
    user = models.ForeignKey(User, on_delete=CASCADE)
    award_text = models.TextField(null=False, blank=False)
    award_level = models.TextField(max_length=1, choices=UserAwardLevel.choices, default=UserAwardLevel.GOLD)
    active = models.BooleanField(null=False, blank=True, default=True)
    kind = models.CharField(max_length=1, choices=UserAwardKind.choices, default=UserAwardKind.KUDOS)
    definition_key = models.TextField(null=True, blank=True)  # AwardDefinition.key, null for kudos
    subject = models.TextField(null=True, blank=True)  # e.g. tag name for per-tag titles
    period = models.CharField(max_length=1, choices=AwardPeriod.choices, null=True, blank=True)  # titles only
    count = models.IntegerField(null=True, blank=True)  # score (titles) / raw progress (badges)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["user", "definition_key", "subject", "period"],
                                    condition=models.Q(definition_key__isnull=False),
                                    nulls_distinct=False,
                                    name="user_award_unique_computed"),
        ]

    @property
    def level_css_class(self) -> str:
        return f"user-award-{self.get_award_level_display().lower()}"

    @property
    def icon_name(self) -> str:
        """ font-awesome icon (without the fa-solid prefix) """
        if self.kind == UserAwardKind.TITLE and self.period:
            return AwardPeriod(self.period).icon
        if self.kind == UserAwardKind.BADGE:
            if definition := self.definition:
                return definition.icon
            return "fa-award"
        return "fa-trophy"

    @property
    def definition(self):
        """ The AwardDefinition this was computed from (None for kudos or a retired definition) """
        if self.definition_key:
            return get_award_definition(self.definition_key)
        return None

    @property
    def icon_class(self):
        return SafeString(f"fa-solid {self.icon_name} user-award {self.level_css_class}")

    @property
    def icon(self):
        return SafeString(f"<i class='{self.icon_class}'></i>")

    @property
    def period_rank(self) -> int:
        return AwardPeriod(self.period).rank if self.period else 0

    def __str__(self):
        return f"{self.get_award_level_display()} {self.user}: {self.award_text}"


class UserAwards:
    """ A user's awards, partitioned by kind. Only titles decorate the user elsewhere on the site
        (see AvatarDetails) - badges and kudos live on the profile pages """

    def __init__(self, user: User):
        award_qs = UserAward.objects.filter(user=user).all()
        award_list: list[UserAward] = sorted(award_qs, key=lambda x: (not x.active, -x.period_rank, 100 - UserAwardLevel(x.award_level).int_value, x.award_text))

        self.all_awards = award_list
        self.awards = [award for award in award_list if award.active]

    def __bool__(self):
        return bool(self.awards)

    @cached_property
    def titles(self) -> list[UserAward]:
        """ Active titles, ALL_TIME -> MONTH -> DAY """
        return [a for a in self.awards if a.kind == UserAwardKind.TITLE]

    @cached_property
    def badges(self) -> list[UserAward]:
        """ Earned badges only - see badge_row() for the progress of unearned ones """
        return [a for a in self.awards if a.kind == UserAwardKind.BADGE]

    @cached_property
    def kudos(self) -> list[UserAward]:
        return [a for a in self.awards if a.kind == UserAwardKind.KUDOS]

    def badge_row(self, definition_key: str) -> Optional[UserAward]:
        """ The (possibly not yet earned, hence inactive) badge row for a definition """
        return first((a for a in self.all_awards if a.kind == UserAwardKind.BADGE and a.definition_key == definition_key), None)

    def progress(self, definition) -> "BadgeProgress":
        return BadgeProgress(definition=definition, award=self.badge_row(definition.key))

    @cached_property
    def highest_award(self) -> Optional[UserAward]:
        return first(self.awards, None)

    @property
    def award_text_html(self):
        return "<br/>".join([f"<i class='{award.icon_class}'></i>" + escape(award.award_text) for award in self.awards])

    @property
    def html(self):
        if not self.awards:
            return ""
        return SafeString(f'<i class="{self.highest_award.icon_class}" title="{ self.award_text_html }"></i>')

    def __iter__(self):
        return iter(self.awards)

    def __len__(self):
        return len(self.awards)

    def __getitem__(self, item):
        return self.awards[item]


@dataclass(frozen=True)
class BadgeProgress:
    """ Where a user is on a badge's bronze/silver/gold ladder - for the profile page """
    definition: object  # AwardDefinition
    award: Optional[UserAward]

    @property
    def count(self) -> int:
        return (self.award.count if self.award else None) or 0

    @property
    def earned(self) -> bool:
        return bool(self.award and self.award.active)

    @property
    def next_threshold(self) -> Optional[int]:
        """ The next tier's threshold, None once gold """
        return first((t for t in self.definition.tiers if t > self.count), None)

    @property
    def next_tier_name(self) -> Optional[str]:
        if (threshold := self.next_threshold) is None:
            return None
        return ["bronze", "silver", "gold"][self.definition.tiers.index(threshold)]

    @property
    def percent(self) -> int:
        if (threshold := self.next_threshold) is None:
            return 100
        return min(100, int(100 * self.count / threshold))

    @property
    def visible(self) -> bool:
        """ Hidden definitions (Night Owl etc) only show once earned """
        return self.earned or not self.definition.hidden


class LabProject(models.Model):
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    leader = models.TextField()
    members = models.TextField(blank=True)
    families = models.IntegerField(default=0)
    involved = models.BooleanField(default=False)
    date = models.DateField()


class SiteMessageDict(TypedDict, total=False):
    severity: str
    message: str
    timestamp: int


class SiteMessage(models.Model):
    """ Used to display a message at the top of the site (eg warnings/messages) """

    # log_level not currently used.
    log_level = models.CharField(max_length=1, choices=LogLevel.CHOICES, default=LogLevel.INFO)
    message = models.TextField()
    date_time = models.DateTimeField(null=True, blank=True)

    @staticmethod
    def _get_site_messages() -> list[SiteMessageDict]:
        site_messages: list[SiteMessageDict] = []
        if settings.SITE_MESSAGE:
            site_messages.append({
                "message": settings.SITE_MESSAGE
            })

        sm: SiteMessage
        for sm in SiteMessage.objects.order_by("pk"):
            site_messages.append({
                "message": sm.message,
                "severity": sm.log_level,
                "timestamp": sm.date_time.timestamp() if sm.date_time else None
            })

        return site_messages

    @staticmethod
    def get_site_messages() -> list[SiteMessageDict]:
        """ This is cached as otherwise we'd be introducing a DB query for every page """
        SITE_MESSAGE_KEY = "variantgrid_site_message"
        site_message_dict = None
        site_message_str = cache.get(SITE_MESSAGE_KEY)
        reload = True
        if site_message_str is not None:
            try:
                site_message_dict = json.loads(site_message_str)
                reload = False
            except Exception:
                pass

        if reload:
            site_message_dict = SiteMessage._get_site_messages()
            site_message_str = json.dumps(site_message_dict)
            cache.set(SITE_MESSAGE_KEY, site_message_str, timeout=30)
        return site_message_dict

    def __str__(self):
        return f"{self.get_log_level_display()}: {self.message}"
