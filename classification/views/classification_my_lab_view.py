from typing import Optional, Union

from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render
from library.utils import ExportRow, export_column

from classification.classification_stats import get_vus_lab_gene_counts
from snpdb.lab_picker import LabPickerData

def view_my_lab(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'my_lab_lab')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/my_lab.html", {"lab_picker_data": lab_picker})


def view_my_lab_detail(request: HttpRequest, lab_id: Optional[Union[str, int]] = None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id)
    labs = lab_picker.lab_selection.selected_labs

    gene_vus_count = get_vus_lab_gene_counts(user=request.user, labs=labs, allele_level=False)

    vus_present = any(
        len(d.get("x", [])) > 0
        for d in gene_vus_count
    )

    gene_tuples = zip(gene_vus_count[0]['x'], gene_vus_count[0]['y'])

    return render(request,
                  "classification/my_lab_detail.html",
                  {"user": request.user,
                   "lab_picker_data": lab_picker,
                   "vus_present": vus_present,
                   "gene_vus_count": gene_vus_count,
                   "gene_tuples": gene_tuples})


class GeneClinSigCountDownload(ExportRow):

    def __init__(self, gene_counts: tuple):
        self.gene_counts: tuple = gene_counts

    @export_column(label="gene")
    def gene(self):
        return self.gene_counts[0]

    @export_column(label="vus")
    def vus(self):
        return self.gene_counts[1]

    @export_column(label="benign")
    def ben(self):
        return self.gene_counts[2]

    @export_column(label="likely_benign")
    def likely_ben(self):
        return self.gene_counts[3]

    @export_column(label="likely_pathogenic")
    def likely_pat(self):
        return self.gene_counts[4]

    @export_column(label="pathogenic")
    def pat(self):
        return self.gene_counts[5]


def my_lab_download(request: HttpRequest, lab_id: Optional[Union[str, int]]= None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id)
    labs = lab_picker.lab_selection.selected_labs

    gene_vus_count = get_vus_lab_gene_counts(user=request.user, labs=labs, allele_level=False)
    all_counts = [g['y'] for g in gene_vus_count]
    gene_tuples = [ (x, *ys) for x, *ys in zip(gene_vus_count[0]['x'], *all_counts)]

    return GeneClinSigCountDownload.streaming_csv(iter(gene_tuples),
                                                  filename="trial_gene_vus")