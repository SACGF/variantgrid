"""

This is the main app for VariantGrid, and this file contains the main models for Variants / Samples
etc and things that don't fit anywhere else.

'snpdb' was the highly unoriginal name I used before 'VariantGrid'
"""
import json
from functools import total_ordering

from celery.result import AsyncResult
from django.conf import settings
from django.contrib.auth.models import User, Group
from django.core.cache import cache
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models import QuerySet
from django.db.models.aggregates import Min, Count
from django.db.models.deletion import SET_NULL, CASCADE
from django.urls import reverse
from django.utils import timezone
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from lazy import lazy
import logging
from model_utils.managers import InheritanceManager
from typing import List, TypedDict, Optional, Dict
from library.enums.log_level import LogLevel
from library.enums.time_enums import TimePeriod
from library.log_utils import send_notification
from library.utils import import_class
from classification.enums.classification_enums import ShareLevel


class Tag(models.Model):
    id = models.CharField(max_length=50, primary_key=True)

    def __str__(self):
        return self.id


class CachedGeneratedFile(models.Model):
    """ Matplotlib graph generated via async Celery task
        @see snpdb.graphs.graphcache.CacheableGraph """
    filename = models.TextField(null=True)
    exception = models.TextField(null=True)
    generator = models.TextField()
    params_hash = models.TextField()
    task_id = models.CharField(max_length=36, null=True)
    task_status = models.TextField(null=True)  # TODO: what's the actual size?
    generate_start = models.DateTimeField(null=True)
    generate_end = models.DateTimeField(null=True)

    class Meta:
        unique_together = ("generator", "params_hash")

    def __str__(self):
        if self.filename:
            description = f"file: {self.filename} made: {self.generate_end}"
        else:
            description = f"task: {self.task_id} sent: {self.generate_start}"
        return f"{self.generator}({self.params_hash}): {description}"

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

    def save(self, **kwargs):
        if settings.CACHE_GENERATED_FILES:
            logging.debug("Save cachedgeneratedfile!")
            super().save(**kwargs)
        else:
            logging.debug("Not caching generated files - skipped save!")


class Company(models.Model):
    name = models.TextField(primary_key=True)

    @staticmethod
    def get_our_company():
        company = None
        company_name = getattr(settings, "COMPANY", None)
        if company_name:
            company, _ = Company.objects.get_or_create(name=company_name)
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

    class Meta:
        abstract = True


class Wiki(TimeStampedModel):
    # Inherit and add a 1-to-1 field for whatever object you want
    objects = InheritanceManager()
    last_edited_by = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    markdown = models.TextField()

    def _get_restricted_object(self):
        return None

    def can_write(self, user):
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
    def get_or_create(class_name, unique_keyword, unique_value):
        klass = Wiki.get_subclass_by_name(class_name)
        wiki, _ = klass.objects.get_or_create(**{unique_keyword: unique_value})
        return wiki

    def __str__(self):
        return self.markdown


class Organization(models.Model):
    # If you add fields @see OrganizationAdmin
    name = models.TextField()
    short_name = models.TextField(blank=False, null=True)  # Don't use for anything other than human labels
    group_name = models.TextField(blank=True, null=True, unique=True)
    classification_config = models.JSONField(null=True, blank=True)
    active = models.BooleanField(default=True, blank=True)

    class Meta:
        ordering = ['name']
        verbose_name = 'Organisation'
        verbose_name_plural = 'Organisations'

    def get_absolute_url(self):
        return reverse('view_organization', kwargs={"pk": self.pk})

    def __str__(self):
        return self.name

    def shortest_name(self):
        return self.short_name or self.name

    @lazy
    def classifying_labs(self) -> List['Lab']:
        return [lab for lab in self.lab_set.all().order_by('name') if lab.total_classifications > 0]

    @lazy
    def sharing_labs(self) -> List['Lab']:
        return [lab for lab in self.lab_set.all().order_by('name') if lab.total_shared_classifications > 0]

    def is_member(self, user: User):
        return Lab.valid_labs_qs(user).filter(organization=self).exists()

    def can_write(self, user: User):
        return user.is_superuser or self.lab_set.filter(labhead__user=user).exists()

    def check_can_write(self, user):
        if not self.can_write(user):
            msg = f"You do not have WRITE permission for {self.pk}"
            raise PermissionDenied(msg)


@total_ordering
class LabUser:

    def __init__(self, user: User, role: str):
        self.user = user
        self.role = role

    @lazy
    def preferred_label(self) -> str:
        from snpdb.models import UserSettings
        return UserSettings.preferred_label_for(self.user)

    def __lt__(self, other):
        return self.preferred_label < other.preferred_label

    def __eq__(self, other):
        return self.user == other.user


class Lab(models.Model):
    name = models.TextField()
    external = models.BooleanField(default=False, blank=True)  # From somewhere else, eg Shariant
    city = models.TextField()
    state = models.TextField(null=True)
    country = models.TextField()
    url = models.TextField(blank=True)
    css_class = models.TextField(blank=True)
    lat = models.FloatField(null=True, blank=True)
    long = models.FloatField(null=True, blank=True)

    group_name = models.TextField(blank=True, null=True, unique=True)
    classification_config = models.JSONField(null=True, blank=True)

    # want every lab to have an organization, but not going to have them
    # at point of migration
    organization = models.ForeignKey(Organization, null=False, blank=False, on_delete=CASCADE)
    # location where the lab can upload files to, (in some environments may refer to s3 directory)
    upload_location = models.TextField(null=True, blank=True)

    email = models.TextField(blank=True)
    slack_webhook = models.TextField(blank=True)

    def send_notification(self,
            message: str,
            blocks: Optional[Dict] = None,
            username: Optional[str] = None,
            emoji: str = ":dna:"):
        if slack_url := self.slack_webhook:
            send_notification(message=message, blocks=blocks, username=username, emoji=emoji, slack_webhook_url=slack_url)
        else:
            raise ValueError("Lab is not configured to send Slack notifications")

    class Meta:
        ordering = ['name']

    @property
    def group(self):
        if self.group_name:
            group, _ = Group.objects.get_or_create(name=self.group_name)
            return group
        return None

    @property
    def group_institution(self):
        if self.group_name:
            parts = self.group_name.split('/')
            if len(parts) >= 2:
                inst_group_name = '/'.join(parts[:-1])
                group, _ = Group.objects.get_or_create(name=inst_group_name)
                return group

    @lazy
    def active_users(self):
        return self.group.user_set.filter(is_active=True)

    @lazy
    def lab_users(self):
        users = list(self.group.user_set.filter(is_active=True))
        heads = set(self.labhead_set.values_list('user_id', flat=True))
        lab_users: List[LabUser] = list()
        for user in users:
            role = 'user'
            if user.id in heads:
                role = 'head'
            lab_users.append(LabUser(user=user, role=role))
        lab_users.sort()
        return lab_users

    @lazy
    def shared_classifications(self):
        return self.classification_set.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).exclude(withdrawn=True)

    @lazy
    def classifications(self):
        """ Shared or all classifications based on settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED """
        if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
            qs = self.shared_classifications
        else:
            qs = self.classification_set.all()
        return qs

    @lazy
    def total_classifications(self):
        """ Total 'classifications' as per above method """
        return self.classifications.count()

    @lazy
    def total_shared_classifications(self) -> int:
        return self.shared_classifications.count()

    @lazy
    def total_unshared_classifications(self) -> int:
        all_classifications = self.classification_set.exclude(withdrawn=True).count()
        return all_classifications - self.total_shared_classifications

    @lazy
    def classifications_by_created(self) -> QuerySet:
        return self.classifications.order_by("created")

    @lazy
    def classifications_by_modified(self) -> QuerySet:
        return self.classifications.order_by("-modified")

    @property
    def first_classification_ever_shared(self) -> 'Classification':
        return self.classification_set.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).order_by('-created').first()

    @lazy
    def classifications_per_day(self):
        try:
            data = self.classifications.aggregate(earliest=Min("created"))
            earliest = data["earliest"]
            latest = now()
            days = (latest - earliest).days
            if days == 0:
                days = 1
            cpd = self.total_classifications / days
        except:
            cpd = 0
        return cpd

    @lazy
    def classifications_per_week(self):
        return self.classifications_per_day * 7

    @staticmethod
    def valid_labs_qs(user: User, admin_check=False) -> QuerySet:
        if admin_check and user.is_superuser:
            return Lab.objects.all()

        group_names = list(user.groups.values_list('name', flat=True))
        return Lab.objects.filter(group_name__in=group_names)

    def classifications_activity(self, time_period):
        trunc_func = TimePeriod.truncate_func(time_period)
        qs = self.classifications.annotate(time_period=trunc_func("created")).values("time_period") \
            .annotate(num_classifications=Count("id")).order_by("-time_period")
        return qs.values_list("time_period", "num_classifications")

    def classifications_activity_by_day(self):
        return self.classifications_activity(TimePeriod.DAY)

    def is_member(self, user: User, admin_check=False):
        return self.valid_labs_qs(user=user, admin_check=admin_check).filter(pk=self.pk).exists()

    def can_write(self, user: User):
        return user.is_superuser or self.labhead_set.filter(user=user).exists()

    def check_can_write(self, user):
        if not self.can_write(user):
            msg = f"You do not have WRITE permission for {self.pk}"
            raise PermissionDenied(msg)

    def save(self, **kwargs):
        super().save(**kwargs)
        if self.group_name:
            # pre-create the groups
            self.group
            self.group_institution
            # TODO assign the lab permissions to the groups
            # but also make sure this doesn't break RunX1

    def get_absolute_url(self):
        return reverse('view_lab', kwargs={"pk": self.pk})

    def __str__(self):
        name = self.name
        if self.organization:
            name += f" ({self.organization})"
        return name


class LabHead(models.Model):
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)

    class Meta:
        unique_together = ('lab', 'user')

    def __str__(self):
        return f"{self.lab}: {self.user}"


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
    def _get_site_messages() -> List[SiteMessageDict]:
        site_messages: List[SiteMessageDict] = []
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
    def get_site_messages() -> List[SiteMessageDict]:
        """ This is cached as otherwise we'd be introducing a DB query for every page """
        SITE_MESSAGE_KEY = "variantgrid_site_message"
        site_message_dict = None
        site_message_str = cache.get(SITE_MESSAGE_KEY)
        reload = True
        if site_message_str is not None:
            try:
                site_message_dict = json.loads(site_message_str)
                reload = False
            except:
                pass

        if reload:
            site_message_dict = SiteMessage._get_site_messages()
            site_message_str = json.dumps(site_message_dict)
            cache.set(SITE_MESSAGE_KEY, site_message_str, timeout=30)
        return site_message_dict
