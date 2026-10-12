"""
Recompute of badges from the registered AwardDefinitions (#1819) - run by the nightly beat task in
snpdb/tasks/user_award_tasks.py and manage.py update_user_awards.
"""
import logging
from typing import Optional

from django.conf import settings
from django.utils import timezone

from library.guardian_utils import admin_bot
from snpdb.models import UserAward
from snpdb.models.models_enums import UserAwardKind, UserAwardLevel
from snpdb.user_awards import AwardDefinition, get_award_definitions

BADGE_LEVELS = (UserAwardLevel.BRONZE, UserAwardLevel.SILVER, UserAwardLevel.GOLD)


def update_user_awards(definitions: Optional[list[AwardDefinition]] = None) -> bool:
    """ Returns False when awards are off for the deployment. definitions defaults to the registry """
    if not settings.USER_AWARDS_ENABLED:
        return False
    if definitions is None:
        definitions = get_award_definitions()
    for definition in definitions:
        update_badge(definition)
    deactivate_retired_definitions()
    return True


def update_badge(definition: AwardDefinition):
    """ Tier = highest threshold <= count. Below bronze the row is kept inactive so the profile can
        show progress. Badges are never revoked and the tier only moves up """
    if not definition.enabled:
        return
    counts = definition.compute()
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
                                     award_text=definition.title)
            continue
        changed = award.count != count
        award.count = count
        if level is not None and (not award.active or UserAwardLevel(award.award_level) < level):
            award.active = True
            award.award_level = level
            award.award_text = definition.title
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
