from datetime import datetime
from datetime import timedelta, date
from dataclasses import dataclass, field
from django.http import StreamingHttpResponse, HttpRequest, HttpResponseBase, JsonResponse
from library.utils import ExportRow, export_column
from typing import Iterator


@dataclass
class InternalLabRow(ExportRow):

    @export_column(label="Share Level")
    def sharelevel_column(self) -> str:
        return 'self.share_level'

    @export_column(label="Genome Build")
    def build_column(self) -> str:
        return 'build'

    @export_column(label="Allele Origin")
    def allele_origin_column(self) -> str:
        return 'allele_origin'

    @export_column(label="Genomic Location")
    def genomic_locations_column(self) -> str:
        return 'genomic_locations'


def stream_internal_lab_rows(share_level, build, allele_origin, genomic_locations) -> Iterator[InternalLabRow]:
    print('fields data', genomic_locations, share_level, build, allele_origin)
    report_data = InternalLabRow()
    yield report_data


def internal_lab_download(request):
    if request.method == 'POST':
        share_level = request.POST.get('share_level')
        build = request.POST.get('genome_build')
        allele_origin = request.POST.get('allele_origin')
        genomic_locations = request.POST.get('genomic_location')

        return InternalLabRow.streaming_csv(
            stream_internal_lab_rows(share_level, build, allele_origin, genomic_locations),
            filename="Internal_lab_report.csv"
        )
    else:
        return JsonResponse({'status': 'error'})
