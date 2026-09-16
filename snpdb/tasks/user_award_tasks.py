import celery
from django.conf import settings

from snpdb.models.models_enums import AwardPeriod
from snpdb.user_award_updates import update_user_awards as update_user_awards_now


@celery.shared_task
def update_user_awards(periods: list[str], badges: bool = False):
    """ Beat: hourly with [DAY], nightly with [ALL_TIME, MONTH] + badges - see variantgrid/celery.py """
    if not settings.USER_AWARDS_ENABLED:
        return
    update_user_awards_now([AwardPeriod(p) for p in periods], badges=badges)
