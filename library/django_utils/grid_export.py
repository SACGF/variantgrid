import csv
from collections.abc import Iterator

from django.http.response import StreamingHttpResponse
from django.utils.text import slugify

from library.utils import StashFile
from library.utils.date_utils import local_date_string

# Rows written per yielded chunk - one chunk per row costs far more in WSGI/nginx (and file write)
# overhead than it saves in memory
EXPORT_ROWS_PER_CHUNK = 1000

MAX_FILE_NAME_LENGTH = 100  # Shortened as someone got a DDE error opening it in Windows


def csv_streaming_response(basename: str, csv_iterator: Iterator[str]) -> StreamingHttpResponse:
    basename = f"{basename}_{local_date_string()}"[:MAX_FILE_NAME_LENGTH]
    response = StreamingHttpResponse(csv_iterator, content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{slugify(basename)}.csv"'
    return response


def grid_export_csv(columns, items) -> Iterator[str]:
    """ columns: [{name, label}] - @see DatatableConfig.csv_columns """
    pseudo_buffer = StashFile()
    header, labels = _header_labels(columns)
    # Don't use dictwriter as some sample names may be the same
    writer = csv.writer(pseudo_buffer, dialect='excel', quoting=csv.QUOTE_MINIMAL)
    if header and str(header[0]).upper().startswith("ID"):
        # Excel reads a file whose first letters are "ID" as a SYLK file - quoting stops that
        csv.writer(pseudo_buffer, dialect='excel', quoting=csv.QUOTE_ALL).writerow(header)
    else:
        writer.writerow(header)

    yield pseudo_buffer.value
    for i, obj in enumerate(items, start=1):
        # labels dict is in same sorted order as header
        row = [obj.get(k) for k in labels]
        writer.writerow(row)
        if i % EXPORT_ROWS_PER_CHUNK == 0:
            yield pseudo_buffer.value
    if remaining := pseudo_buffer.value:
        yield remaining


def _header_labels(columns):
    labels = {}

    header = []
    for c in columns:
        name = c['name']
        column_label = c.get("label", name)
        labels[name] = column_label
        header.append(column_label)

    return header, labels
