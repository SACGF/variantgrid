"""
The case report templates this repo ships, under classification/test_data/.

The generic one is what a deployment gets before it writes its own: an HTML document with the
header, results summary, tier sections and method block. A deployment's own case template (SA Path's
TSO 500 one) replaces the content of the row, not this file.
"""
import os

TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data")
GENERIC_CASE_TEMPLATE_FILENAME = "generic_case_report.html"


def _read(filename: str) -> str:
    with open(os.path.join(TEST_DATA_DIR, filename), encoding="utf-8") as f:
        return f.read()


def generic_case_template() -> str:
    return _read(GENERIC_CASE_TEMPLATE_FILENAME)


def generic_json_template() -> str:
    """ Nothing has a json_template any more (@see report/renderers.py:render_json) - this survives
        only because the pushed migration 0179_default_case_report_template imports it """
    return ""
