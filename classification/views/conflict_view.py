from typing import Optional, Union

from django.shortcuts import render

from classification.models import Conflict
from classification.views.classification_dashboard_view import ClassificationDashboard
from library.utils.django_utils import render_ajax_view
from snpdb.lab_picker import LabPickerData


def conflict_detail(request, conflict_id: int):
    conflict = Conflict.objects.get(id=conflict_id)
    return render_ajax_view(request, 'classification/conflict_detail.html', context={'conflict': conflict})


# def allele_groupings(request, lab_id: Optional[Union[str, int]] = None):

def conflicts_view(request, lab_id: Optional[Union[str, int]] = None):
    lab_picker = LabPickerData.from_request(request, lab_id, 'conflicts')
    return render(request, 'classification/conflicts.html', context={
        "dlab": ClassificationDashboard(lab_picker=lab_picker)
    })