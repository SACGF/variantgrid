"""
Records what a block of code wrote, per model: rows created, updated, unchanged and deleted, and the
names of the fields an update changed (#1413), plus free-text notes about what it did (#2349).
IntegrationApiMiddleware wraps each tracked API call in collect_changes() and stores the result on its
ApiRequestLog row.

Entry points:
- collect_changes(): context manager that installs a fresh ChangeCollector for the block and yields it
- active_collector(): that collector, or None outside collect_changes()
- update_or_create_if_changed(model, defaults, **lookup): update_or_create that only saves when a default differs
- get_or_create_recorded(model, defaults, **lookup): get_or_create that records a row it found as unchanged
- track_deletes(*models): count deletes of these models

Creates and updates are counted by a post_save receiver for every model. Deletes are opt-in per model,
because a post_delete receiver for every model would stop Django fast-deleting any of them.
bulk_create, QuerySet.update and fast deletes send no signals, so code a tracked call runs that writes
that way records its rows with ChangeCollector.record.
"""
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, Optional

from django.core.exceptions import FieldDoesNotExist
from django.db import models, router, transaction
from django.db.models.signals import post_delete, post_save

IGNORED_APP_LABELS = {"eventlog"}  # the API log never counts its own rows

class WriteOutcome(StrEnum):
    CREATED = "created"
    UPDATED = "updated"
    UNCHANGED = "unchanged"


@dataclass
class ModelChanges:
    created: int = 0
    updated: int = 0
    unchanged: int = 0
    deleted: int = 0
    fields: set[str] = field(default_factory=set)

    @property
    def num_changes(self) -> int:
        return self.created + self.updated + self.deleted

    def as_dict(self) -> dict[str, Any]:
        return {
            "created": self.created,
            "updated": self.updated,
            "unchanged": self.unchanged,
            "deleted": self.deleted,
            "fields": sorted(self.fields),
        }


@dataclass
class ChangeCollector:
    by_model: dict[str, ModelChanges] = field(default_factory=dict)  # key "app_label.model"
    notes: list[str] = field(default_factory=list)

    def record(self, model: type[models.Model], created: int = 0, updated: int = 0, unchanged: int = 0,
               deleted: int = 0, fields: Iterable[str] = ()):
        model_changes = self.by_model.setdefault(model._meta.label_lower, ModelChanges())
        model_changes.created += created
        model_changes.updated += updated
        model_changes.unchanged += unchanged
        model_changes.deleted += deleted
        model_changes.fields.update(fields)

    @property
    def num_changes(self) -> int:
        return sum(mc.num_changes for mc in self.by_model.values())

    def as_dict(self) -> dict[str, dict[str, Any]]:
        """ {"app_label.model": {"created": n, "updated": n, "unchanged": n, "deleted": n, "fields": [sorted]}},
            leaving out models nothing was recorded for """
        return {key: mc.as_dict() for key, mc in sorted(self.by_model.items())
                if mc.num_changes or mc.unchanged}


_active_collector: ContextVar[Optional[ChangeCollector]] = ContextVar("change_collector", default=None)


def active_collector() -> Optional[ChangeCollector]:
    return _active_collector.get()


@contextmanager
def collect_changes() -> Iterator[ChangeCollector]:
    """ A context variable rather than an attribute on the request, as nested serializers have no request """
    collector = ChangeCollector()
    token = _active_collector.set(collector)
    try:
        yield collector
    finally:
        _active_collector.reset(token)


def _tracked_collector(model: type[models.Model]) -> Optional[ChangeCollector]:
    collector = _active_collector.get()
    if collector is None or model._meta.app_label in IGNORED_APP_LABELS:
        return None
    return collector


def _auto_now_field_names(model: type[models.Model]) -> set[str]:
    return {f.name for f in model._meta.concrete_fields if getattr(f, "auto_now", False)}


def _record_save(sender, instance, created, update_fields=None, raw=False, **kwargs):  # pylint: disable=unused-argument
    if raw or not (collector := _tracked_collector(sender)):
        return
    if created:
        collector.record(sender, created=1)
    else:
        fields = set(update_fields or ()) - _auto_now_field_names(sender)
        collector.record(sender, updated=1, fields=fields)


def _record_delete(sender, instance, **kwargs):  # pylint: disable=unused-argument
    if collector := _tracked_collector(sender):
        collector.record(sender, deleted=1)


post_save.connect(_record_save, dispatch_uid="change_tracking_record_save")


def track_deletes(*model_classes: type[models.Model]):
    """ Connected per model - a receiver on that model costs it Django's fast delete """
    for model in model_classes:
        post_delete.connect(_record_delete, sender=model,
                            dispatch_uid=f"change_tracking_record_delete_{model._meta.label_lower}")


def _differs(instance: models.Model, name: str, value) -> Optional[str]:
    """ The field name if setting it to value would change instance (and sets it), else None. A foreign key
        compares by its column, so a model instance and its pk are equal """
    try:
        model_field = instance._meta.get_field(name)
    except FieldDoesNotExist:
        model_field = None

    if model_field is not None and model_field.concrete and model_field.is_relation:
        if isinstance(value, models.Model):
            new_value = getattr(value, model_field.target_field.attname)
            attribute = model_field.name
        else:
            new_value = value
            attribute = model_field.attname
        if getattr(instance, model_field.attname) == new_value:
            return None
        setattr(instance, attribute, value)
        return model_field.name

    if getattr(instance, name) == value:
        return None
    setattr(instance, name, value)
    return name


def update_or_create_if_changed[ModelT: models.Model](
        model: type[ModelT], /, defaults: Optional[dict[str, Any]] = None, **lookup) -> tuple[ModelT, WriteOutcome]:
    """ update_or_create that saves an existing row only when a default differs, and then only those fields
        (plus auto_now ones, so `modified` means changed). Callables in defaults are called, as Django does.
        model is positional-only so a field called "model" can be a lookup """
    defaults = defaults or {}
    with transaction.atomic(using=router.db_for_write(model)):
        instance, created = model.objects.select_for_update().get_or_create(defaults=defaults, **lookup)
        if created:
            return instance, WriteOutcome.CREATED

        changed_fields = []
        for name, value in defaults.items():
            if callable(value):
                value = value()
            if changed_field := _differs(instance, name, value):
                changed_fields.append(changed_field)

        if changed_fields:
            instance.save(update_fields=set(changed_fields) | _auto_now_field_names(model))
            return instance, WriteOutcome.UPDATED

    if collector := _tracked_collector(model):
        collector.record(model, unchanged=1)
    return instance, WriteOutcome.UNCHANGED


def get_or_create_recorded[ModelT: models.Model](
        model: type[ModelT], /, defaults: Optional[dict[str, Any]] = None, **lookup) -> tuple[ModelT, WriteOutcome]:
    """ get_or_create - defaults only apply to a new row - that records a row it found as unchanged """
    instance, created = model.objects.get_or_create(defaults=defaults, **lookup)
    if created:
        return instance, WriteOutcome.CREATED
    if collector := _tracked_collector(model):
        collector.record(model, unchanged=1)
    return instance, WriteOutcome.UNCHANGED
