import celery
from django.conf import settings

from snpdb.user_award_updates import update_user_awards as update_user_awards_now


@celery.shared_task
def update_user_awards():
    """ Beat: nightly - see variantgrid/celery.py """
    if not settings.USER_AWARDS_ENABLED:
        return
    update_user_awards_now()
