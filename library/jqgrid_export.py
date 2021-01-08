import csv
from django.http.response import StreamingHttpResponse
from django.utils.text import slugify


class StashFile:
    """ File-like object that holds a value """

    def __init__(self):
        self.data = ''

    def write(self, value):
        self.data += value

    @property
    def value(self):
        data = self.data
        self.data = ''
        return data


def colmodel_header_labels(colmodels, label_overrides=None):
    labels = {}

    header = []
    for c in colmodels:
        name = c['name']
        column_label = c.get("label", name)
        if label_overrides:
            column_label = label_overrides.get(name, column_label)

        labels[name] = column_label
        header.append(column_label)

    return header, labels


def grid_export_request(request, grid, basename):
    request.GET = request.GET.copy()  # Immutable
    request.GET['rows'] = 0  # No pagination
    items = grid.get_items(request)[2]
    colmodels = grid.get_colmodels()
    return grid_export_csv(basename, colmodels, items)


def grid_export_csv(basename, colmodels, items):
    MAX_FILE_NAME_LENGTH = 100  # Shorted as someone got a DDE error opening it in Windows
    # If you make the first letters “ID” of a text file
    # Excel incorrectly assumes you are trying to open an SYLK file.
    label_overrides = {"id": "variant_id"}

    pseudo_buffer = StashFile()
    header, labels = colmodel_header_labels(colmodels, label_overrides=label_overrides)
    writer = csv.DictWriter(pseudo_buffer, header, dialect='excel')

    def iter_row_writer():
        writer.writeheader()
        yield pseudo_buffer.value
        for obj in items:
            row = {}
            for k, value in obj.items():
                index = labels.get(k)
                if index:
                    row[index] = value
            writer.writerow(row)
            yield pseudo_buffer.value

    basename = basename[:MAX_FILE_NAME_LENGTH]

    response = StreamingHttpResponse(iter_row_writer(), content_type="text/csv")
    response['Content-Disposition'] = 'attachment; filename="%s.csv"' % slugify(basename)
    return response
