"""
/classification/reclassification_analytics: the permission check, bringing the timelines up to date and
rendering. The numbers come from classification/reclassification_analytics.
"""
from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render

from classification.enums.classification_enums import LabExternalFilter
from classification.models.classification_reclassification_models import (
    ReclassificationBuildResult,
    ReclassificationEventBuilder,
)
from classification.reclassification_analytics.analytics import ReclassificationAnalytics
from classification.reclassification_analytics.charts import chart_data
from classification.reclassification_analytics.results import (
    SIGNIFICANCE_ORDER,
    significance_label,
)
from classification.tasks.classification_reclassification_tasks import (
    reclassification_events_update,
)

RECLASSIFICATION_PAGE_BUILD_LIMIT = 5000


def _timelines_up_to_date() -> ReclassificationBuildResult:
    """ Opening the page is what keeps the timelines current, unless there's more waiting than a
        web worker should take on - then the nightly task's queue picks the batch up. """
    build_result = ReclassificationEventBuilder.bring_up_to_date(
        max_classifications=RECLASSIFICATION_PAGE_BUILD_LIMIT)
    if build_result.outstanding:
        reclassification_events_update.delay()
    return build_result


def view_reclassification_analytics(request: HttpRequest) -> HttpResponseBase:
    if not (request.user.is_superuser or settings.CLASSIFICATION_RECLASSIFICATION_ANALYTICS_NON_ADMIN):
        raise PermissionDenied("You must be a super user to view this page")
    build_result = _timelines_up_to_date()
    analytics = ReclassificationAnalytics.from_params(request.GET)
    return render(request, "classification/classification_reclassification_analytics.html", {
        "analytics": analytics,
        "chart_data": chart_data(analytics),
        "build_result": build_result,
        "lab_external_choices": LabExternalFilter.choices,
        "lab_external_default": analytics.lab_external,
        "significance_choices": [(cs, significance_label(cs)) for cs in SIGNIFICANCE_ORDER],
    })
