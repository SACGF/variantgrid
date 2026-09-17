import json

import pandas as pd
from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase

from library.django_utils.datatable_dataframe import (
    INDEX_COLUMN_KEY,
    DataFrameDatatableConfig,
    DataFrameTableView,
)
from snpdb.views.datatable_view import SortOrder


class FakeDataFrameConfig(DataFrameDatatableConfig):
    index_label = "Gene Symbol"
    csv_name = "fake_dataframe"

    def get_dataframe(self) -> pd.DataFrame:
        df = pd.DataFrame({"count": [3, 1, 2], "other.column": ["c", "a", "b"]},
                          index=["GATA2", "BRCA1", "RUNX1"])
        return df.sort_index()

    def get_column_label(self, column_name):
        return f"Label: {column_name}"


class DataFrameTableViewTests(TestCase):
    def setUp(self):
        self.user = User.objects.create(username='dataframe_datatable_user')

    def _view(self, **params) -> DataFrameTableView:
        request = RequestFactory().get("/fake/dataframe/", params)
        request.user = self.user
        view = DataFrameTableView(column_class=FakeDataFrameConfig)
        view.request = request
        view.config = FakeDataFrameConfig(request)
        return view

    def _data(self, **params):
        return self._view(**params).get_context_data()

    def test_definition_columns(self):
        definition = self._view().json_definition()
        columns = definition["columns"]
        # Column names go up positionally, so a name holding a dot isn't read as a nested path
        self.assertEqual([c["data"] for c in columns], [INDEX_COLUMN_KEY, "c1", "c2"])
        self.assertEqual([c["label"] for c in columns],
                         ["Gene Symbol", "Label: count", "Label: other.column"])
        self.assertEqual(definition["order"], [[0, SortOrder.ASC.value]])
        self.assertIn("dataTableCsv=1", definition["downloadUrl"])

    def test_rows_and_index(self):
        data = self._data()
        self.assertEqual(data["recordsTotal"], 3)
        self.assertEqual([row[INDEX_COLUMN_KEY] for row in data["data"]], ["BRCA1", "GATA2", "RUNX1"])
        self.assertEqual([row["c1"] for row in data["data"]], [1, 3, 2])

    def test_sort_by_column(self):
        data = self._data(**{"order[0][column]": "1", "order[0][dir]": "desc"})
        self.assertEqual([row[INDEX_COLUMN_KEY] for row in data["data"]], ["GATA2", "RUNX1", "BRCA1"])

    def test_sort_by_index_descending(self):
        data = self._data(**{"order[0][column]": "0", "order[0][dir]": "desc"})
        self.assertEqual([row[INDEX_COLUMN_KEY] for row in data["data"]], ["RUNX1", "GATA2", "BRCA1"])

    def test_paging(self):
        data = self._data(start=1, length=1)
        self.assertEqual(data["recordsTotal"], 3)
        self.assertEqual([row[INDEX_COLUMN_KEY] for row in data["data"]], ["GATA2"])

    def test_csv_download(self):
        response = self._view(dataTableCsv=1).download_csv()
        csv_text = b"".join(response.streaming_content).decode()
        lines = csv_text.strip().splitlines()
        self.assertEqual(lines[0], "Gene Symbol,Label: count,Label: other.column")
        self.assertEqual(lines[1], "BRCA1,1,a")
        self.assertEqual(len(lines), 4)  # header + every row, not just the first page

    def test_json_serializable(self):
        json.dumps(self._data()["data"])
