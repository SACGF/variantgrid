from django.db.models import QuerySet
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from assay_detailed.models import AssayDetailedRNA
from library.utils.django_utils import render_ajax_view
from snpdb.views.datatable_view import DatatableConfig, DC, RichColumn


class AssayDetailedRNAColumns(DatatableConfig[AssayDetailedRNA]):

    def get_initial_queryset(self) -> QuerySet[DC]:
        return AssayDetailedRNA.objects.all()

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('view_assay_detailed_rna_detail',
                                                                       expected_height=108)

        self.rich_columns = [
            RichColumn("lab"),
            RichColumn("date", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn("id", visible=False)
        ]


def view_assay_detailed_rna_detail(request: HttpRequest, assay_detailed_rna_pk: int) -> HttpResponse:
    full_detail = request.GET.get("detail") == "full"
    context = {
        "assay": AssayDetailedRNA.objects.get(pk=assay_detailed_rna_pk),
        "detail": "full" if full_detail else "partial"
    }
    return render_ajax_view(
        request,
        template_name="assay_detailed/assay_detailed_rna_detail.html",
        menubar="classification",
        context=context)