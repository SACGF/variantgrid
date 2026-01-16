from typing import Optional, Union

from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render
from library.utils import ExportRow, export_column

from classification.classification_stats import get_lab_clinsig_gene_counts
from snpdb.lab_picker import LabPickerData

def view_my_lab(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'my_lab_lab')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/my_lab.html", {"lab_picker_data": lab_picker})


def view_my_lab_detail(request: HttpRequest, lab_id: Optional[Union[str, int]] = None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id)
    labs = lab_picker.lab_selection.selected_labs


    # if allele_level := bool(request.GET.get('allele_level')):
    #     gene_vus_count = get_lab_clinsig_gene_counts(user=request.user, labs=labs, allele_level=allele_level)
    # else:
    #     gene_vus_count = get_lab_clinsig_gene_counts(user=request.user, labs=labs)

    gene_clinsig_count = get_lab_clinsig_gene_counts(user=request.user, labs=labs)
    unique_vus_count = get_lab_clinsig_gene_counts(user=request.user, labs=labs, allele_level=True)


    vus_present = any(
        len(d.get("x", [])) > 0
        for d in gene_clinsig_count
    )

    return render(request,
                  "classification/my_lab_detail.html",
                  {"lab_picker_data": lab_picker,
                   "vus_present": vus_present,
                   "gene_vus_count": gene_clinsig_count,
                   "unique_vus_count": unique_vus_count})


class GeneClinSigCountDownload(ExportRow):

    def __init__(self, gene_counts: tuple):
        self.gene_counts: tuple = gene_counts

    @export_column(label="Gene")
    def gene(self):
        return self.gene_counts[0]

    @export_column(label="Total VUS")
    def total_vus(self):
        return self.gene_counts[1]

    @export_column(label="VUS")
    def vus(self):
        return self.gene_counts[2]

    @export_column(label="VUS A")
    def vus_a(self):
        return self.gene_counts[3]

    @export_column(label="VUS B")
    def vus_b(self):
        return self.gene_counts[4]

    @export_column(label="VUS C")
    def vus_c(self):
        return self.gene_counts[5]

    @export_column(label="Benign")
    def ben(self):
        return self.gene_counts[6]

    @export_column(label="Likely Benign")
    def likely_ben(self):
        return self.gene_counts[7]

    @export_column(label="Likely Pathogenic")
    def likely_pat(self):
        return self.gene_counts[8]

    @export_column(label="Pathogenic")
    def pat(self):
        return self.gene_counts[9]

class GeneUniqueVUSCountDownload(ExportRow):
    def __init__(self, gene_counts: tuple):
        self.gene_counts: tuple = gene_counts

    @export_column(label="Gene")
    def gene(self):
        return self.gene_counts[0]

    @export_column(label="Unique VUS")
    def vus(self):
        return self.gene_counts[1]

def my_lab_download(request: HttpRequest, lab_id: Optional[Union[str, int]]= None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id)
    labs = lab_picker.lab_selection.selected_labs

    if allele_level := bool(request.GET.get('allele_level')):
        print(allele_level)
        gene_vus_count = get_lab_clinsig_gene_counts(user=request.user, labs=labs, allele_level=allele_level)
    else:
        gene_vus_count = get_lab_clinsig_gene_counts(user=request.user, labs=labs)

    all_counts = [g['y'] for g in gene_vus_count]
    gene_tuples = [ (x, *ys) for x, *ys in zip(gene_vus_count[0]['x'], *all_counts)]
    print(gene_tuples)

    if allele_level:
        return GeneUniqueVUSCountDownload.streaming_csv(iter(gene_tuples),
                                                        filename="unique_VUS")
    else:
        return GeneClinSigCountDownload.streaming_csv(iter(gene_tuples),
                                                  filename="trial_gene_vus")