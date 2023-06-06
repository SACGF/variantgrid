from typing import Optional, Union

from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render

from classification.models.allele_overlaps import OverlapsCalculator
from snpdb.lab_picker import LabPickerData


def view_overlaps_vus(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'vus')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/vus.html", {"lab_picker_data": lab_picker})


def view_overlaps_vus_detail(request: HttpRequest, lab_id: Optional[Union[str, int]] = None) -> HttpResponseBase:
    return render(request, "classification/vus_detail.html", {
        "overlaps": OverlapsCalculator(perspective=LabPickerData.from_request(request, lab_id), shared_only=True)
    })
