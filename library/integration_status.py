from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Any, Callable, Optional

import django.dispatch
from django.db import models
from django.utils.text import slugify

from library.utils import flatten_nested_lists

"""
Integration statuses answer "when did each external system last check in?" - they're collected
live for the Server Status page, where health checks are a nightly Slack digest.

Apps register a receiver on integration_status_signal from their AppConfig.ready(), so a
deployment's integrations arrive with its apps and core never imports them.
"""


class IntegrationDirection(models.TextChoices):
    INBOUND = "I", "Inbound"  # something outside pushes, or we pull into VariantGrid
    OUTBOUND = "O", "Outbound"  # VariantGrid calls out


@dataclass
class IntegrationDetail:
    """ One "Last Run: 4 minutes ago" item inside an integration's box """
    label: str
    timestamp: Optional[datetime] = None
    status: str = "info"  # bootstrap contextual class - info/success/warning/danger
    extra: Optional[str] = None  # e.g. the error message against "Last Failure"

    @property
    def is_bad(self) -> bool:
        return self.status == "danger"


@dataclass
class IntegrationTrigger:
    """ A "Sync now" button. The callable is held here rather than named in the POST, so the
        page can only ever run something an app registered """
    action_id: str  # unique across the deployment, e.g. "sapath-helix-sync"
    run: Callable[[], Any]  # usually some_task.delay
    label: str = "Sync now"


@dataclass
class IntegrationStatus:
    name: str  # "Helix NGS Database"
    details: list[IntegrationDetail] = field(default_factory=list)
    key: Optional[str] = None  # identifies the same integration across providers - see fallback
    direction: IntegrationDirection = IntegrationDirection.INBOUND
    description: Optional[str] = None  # one line under the title
    url: Optional[str] = None  # link to the records this integration writes
    record_count: Optional[int] = None
    trigger: Optional[IntegrationTrigger] = None
    warning_age: Optional[timedelta] = None  # feeds the nightly Slack digest
    enabled: bool = True  # renders greyed with a "disabled" badge
    sort_order: int = 0
    fallback: bool = False
    """ A shipped default row, superseded when an app registers a richer status for the same key -
        that's how IntegrationActivity rows appear unaided but still let their app dress them up """

    @property
    def slug(self) -> str:
        return slugify(self.name)

    @property
    def last_activity(self) -> Optional[datetime]:
        """ The newest timestamp across all details, None if this integration has never run """
        timestamps = [detail.timestamp for detail in self.details if detail.timestamp]
        if timestamps:
            return max(timestamps)
        return None

    @property
    def has_run(self) -> bool:
        return self.last_activity is not None

    @property
    def ok(self) -> bool:
        return not any(detail.is_bad for detail in self.details)


integration_status_signal = django.dispatch.Signal()


def get_integration_statuses() -> tuple[list[IntegrationStatus], list[str]]:
    """
    Collected with send_robust: returns the statuses plus one message per receiver that raised,
    so a broken provider shows as a message on the page instead of a 500
    """
    results = []
    messages: list[str] = []
    for caller, result in integration_status_signal.send_robust(sender=None):
        if isinstance(result, Exception):
            messages.append(f"Exception generating integration status by {caller}: {result}")
        else:
            results.append(result)

    statuses: list[IntegrationStatus] = flatten_nested_lists(results)
    claimed_keys = {status.key for status in statuses if status.key and not status.fallback}
    statuses = [status for status in statuses if not (status.fallback and status.key in claimed_keys)]
    statuses.sort(key=lambda status: (status.sort_order, status.name))
    return statuses, messages


def run_integration_trigger(action_id: str) -> Optional[IntegrationStatus]:
    """ Re-collects the statuses, finds the one whose trigger matches, calls run() """
    statuses, _ = get_integration_statuses()
    for status in statuses:
        if (trigger := status.trigger) and trigger.action_id == action_id:
            trigger.run()
            return status
    return None
