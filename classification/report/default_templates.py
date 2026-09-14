"""
The case report templates this repo ships, under classification/test_data/.

The generic pair is what a deployment gets before it writes its own: an HTML document with the
header, results summary, tier sections and method block, and a blank JSON template, which means the
canonical context dump. A deployment's own templates (SA Path's TSO 500 pair) replace the content of
the row, not these files.
"""
import os

TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data")
GENERIC_CASE_TEMPLATE_FILENAME = "generic_case_report.html"
GENERIC_JSON_TEMPLATE_FILENAME = "generic_case_report.json"


def _read(filename: str) -> str:
    with open(os.path.join(TEST_DATA_DIR, filename), encoding="utf-8") as f:
        return f.read()


def generic_case_template() -> str:
    return _read(GENERIC_CASE_TEMPLATE_FILENAME)


def generic_json_template() -> str:
    """ Blank - the canonical context dump is the default structured record """
    return _read(GENERIC_JSON_TEMPLATE_FILENAME)
