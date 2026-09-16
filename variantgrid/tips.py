""" Feature tips shown on loading screens and blank grids - see variantgrid/data/tips.csv """
import csv
import logging
import os
from dataclasses import dataclass
from typing import Any, Optional

from django.conf import settings

from library.cache import timed_cache
from variantgrid.perm_path import get_visible_url_names

TIPS_CSV = os.path.join(settings.VARIANTGRID_APP_DIR, "data", "tips.csv")


@dataclass(frozen=True)
class Tip:
    tip: str
    app: str
    url: str  # Registered URL name for the page the tip is about - used as the visibility key
    settings_enabled: str  # Dotted settings path that must be truthy (blank = always applicable)


@timed_cache()
def _load_tips() -> list[Tip]:
    tips = []
    with open(TIPS_CSV, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                tips.append(Tip(tip=row["tip"].strip(), app=row["app"].strip(),
                                url=row["url"].strip(), settings_enabled=row["settings_enabled"].strip()))
            except (KeyError, AttributeError):
                # A typo in the data file skips that row rather than killing every tip
                logging.warning("Skipping unparseable row in %s: %s", TIPS_CSV, row)
    return tips


def _setting_enabled(path: str) -> bool:
    """ Walks a dotted path (eg 'SOMALIER.enabled') through settings, stepping into dicts and objects """
    if not path:
        return True

    value: Any = settings
    for step in path.split("."):
        if isinstance(value, dict):
            value = value.get(step)
        else:
            value = getattr(value, step, None)
        if value is None:
            return False
    return bool(value)


@timed_cache()
def get_tips(app: Optional[str] = None) -> list[Tip]:
    """ Tips applicable to this deployment - hidden URLs and switched off features drop out """
    if not settings.TIPS_ENABLED:
        return []

    url_name_visible = get_visible_url_names()
    tips = []
    for tip in _load_tips():
        if app and tip.app != app:
            continue
        if not url_name_visible[tip.url]:
            continue
        if not _setting_enabled(tip.settings_enabled):
            continue
        tips.append(tip)
    return tips
