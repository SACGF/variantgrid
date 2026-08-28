import json
import random
from typing import Optional

from django import template
from django.conf import settings

from snpdb.models.models_user_settings import UserSettings
from variantgrid.tips import get_tips

register = template.Library()


@register.inclusion_tag("tips/tip_box.html", takes_context=True)
def tip_box(context, app: Optional[str] = None):
    """ A "Tip: ..." box, optionally restricted to one app's tips - see variantgrid/data/tips.csv

        The server picks the tip shown initially (so it works without JavaScript), and hands the
        whole applicable list to tips.js for the next-tip button. """

    if not settings.TIPS_ENABLED:
        return {}

    user = context.get("user")
    if not (user and user.is_authenticated):
        return {}

    if not UserSettings.get_for_user(user).show_tips:
        return {}

    tips = get_tips(app=app)
    if not tips:
        return {}

    return {
        "tip": random.choice(tips),
        "tips_json": json.dumps([t.tip for t in tips]),
    }
