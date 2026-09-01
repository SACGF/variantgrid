"""
Recompute of titles and badges from the registered AwardDefinitions (#1819) - run by the beat
task in snpdb/tasks/user_award_tasks.py and manage.py update_user_awards.
"""
import logging
from datetime import datetime
from typing import Optional

from django.conf import settings
from django.utils import timezone

from library.guardian_utils import admin_bot
from snpdb.models import UserAward
from snpdb.models.models_enums import AwardPeriod, UserAwardKind, UserAwardLevel
from snpdb.user_awards import AwardDefinition, get_award_definitions

BADGE_LEVELS = (UserAwardLevel.BRONZE, UserAwardLevel.SILVER, UserAwardLevel.GOLD)


def period_start(period: AwardPeriod) -> Optional[datetime]:
    """ None = all time. Calendar month / day in the server timezone """
    if period == AwardPeriod.ALL_TIME:
        return None
    now = timezone.localtime()
    day_start = now.replace(hour=0, minute=0, second=0, microsecond=0)
    if period == AwardPeriod.DAY:
        return day_start
    return day_start.replace(day=1)


def update_user_awards(periods: list[AwardPeriod], badges: bool = False,
                       definitions: Optional[list[AwardDefinition]] = None) -> bool:
    """ Returns False when awards are off for the deployment. definitions defaults to the registry """
    if not settings.USER_AWARDS_ENABLED:
        return False
    if definitions is None:
        definitions = get_award_definitions()
    titles = [d for d in definitions if d.kind == UserAwardKind.TITLE]
    badge_definitions = [d for d in definitions if d.kind == UserAwardKind.BADGE]

    for definition in titles:
        for period in periods:
            period = AwardPeriod(period)
            if period in definition.periods:
                update_title(definition, period)
    if badges:
        for definition in badge_definitions:
            update_badge(definition)
    deactivate_retired_definitions()
    return True


def update_title(definition: AwardDefinition, period: AwardPeriod):
    """ The user(s) with the max count per subject hold the title (ties all hold it), everyone else
        holding it for that subject/period loses it. Nobody holds a title with count 0 """
    if not definition.enabled:
        return
    counts = definition.compute(period_start(period))
    bot_id = admin_bot().pk
    holder_pks = []
    for subject, user_counts in counts.items():
        user_counts = {user_id: c for user_id, c in user_counts.items() if user_id != bot_id and c >= 1}
        if not user_counts:
            continue
        top = max(user_counts.values())
        for user_id, count in user_counts.items():
            if count == top:
                award, _ = UserAward.objects.update_or_create(
                    user_id=user_id, definition_key=definition.key, subject=subject, period=period,
                    defaults={
                        "kind": UserAwardKind.TITLE,
                        "active": True,
                        "count": count,
                        "award_level": UserAwardLevel.GOLD,
                        "award_text": definition.award_text(subject=subject, period=period),
                    })
                holder_pks.append(award.pk)

    lost = UserAward.objects.filter(kind=UserAwardKind.TITLE, definition_key=definition.key, period=period,
                                    active=True).exclude(pk__in=holder_pks)
    if num_lost := lost.update(active=False, modified=timezone.now()):
        logging.info("%s (%s): %d user(s) lost the title", definition.key, period.label, num_lost)


def update_badge(definition: AwardDefinition):
    """ Tier = highest threshold <= count. Below bronze the row is kept inactive so the profile can
        show progress. Badges are never revoked and the tier only moves up """
    if not definition.enabled:
        return
    counts = definition.compute(None).get(None, {})
    bot_id = admin_bot().pk
    existing = {a.user_id: a for a in UserAward.objects.filter(kind=UserAwardKind.BADGE, definition_key=definition.key)}
    for user_id, count in counts.items():
        if user_id == bot_id or count < 1:
            continue
        level = badge_level(definition, count)
        award = existing.get(user_id)
        if award is None:
            UserAward.objects.create(user_id=user_id, kind=UserAwardKind.BADGE, definition_key=definition.key,
                                     active=level is not None, count=count,
                                     award_level=level or UserAwardLevel.BRONZE,
                                     award_text=definition.award_text())
            continue
        changed = award.count != count
        award.count = count
        if level is not None and (not award.active or UserAwardLevel(award.award_level) < level):
            award.active = True
            award.award_level = level
            award.award_text = definition.award_text()
            changed = True
        if changed:
            award.save()


def badge_level(definition: AwardDefinition, count: int) -> Optional[UserAwardLevel]:
    level = None
    for threshold, tier_level in zip(definition.tiers, BADGE_LEVELS):
        if count >= threshold:
            level = tier_level
    return level


def deactivate_retired_definitions():
    """ Awards for definitions that are disabled (settings.USER_AWARDS_DISABLED_KEYS) or no longer
        registered stop showing. They're kept so a re-enabled definition picks up where it left off """
    live_keys = [d.key for d in get_award_definitions() if d.enabled]
    retired = UserAward.objects.filter(definition_key__isnull=False, active=True).exclude(definition_key__in=live_keys)
    if num_retired := retired.update(active=False, modified=timezone.now()):
        logging.info("Deactivated %d award(s) for disabled/retired definitions", num_retired)
