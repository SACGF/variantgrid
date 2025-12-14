from typing import Optional, Union

from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render

from classification.classification_stats import get_classification_counts, get_lab_gene_counts, get_vus_lab_gene_counts
from classification.models.allele_overlaps import OverlapsCalculator
from classification.templatetags.classification_tags import clinical_significance
from snpdb.lab_picker import LabPickerData

from rest_framework.generics import get_object_or_404
from snpdb.models import Lab

def view_my_lab(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'my_lab_lab')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/my_lab.html", {"lab_picker_data": lab_picker})


def view_my_lab_detail(request: HttpRequest, lab_id: Optional[Union[str, int]] = None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id)
    labs = lab_picker.lab_selection.selected_labs

    if len(labs) == 1:
        selected_lab = next(iter(labs))

    gene_vus_count = get_vus_lab_gene_counts(user=request.user, lab=selected_lab)

    return render(request, "classification/my_lab_detail.html", {
        "user": request.user,
        "selected_lab": selected_lab,
        "gene_vus_count": gene_vus_count})