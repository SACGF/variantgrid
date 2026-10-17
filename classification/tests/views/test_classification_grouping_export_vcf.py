"""VCFRow (grouping VCF export): the INFO headers and cells each target system writes."""
from types import SimpleNamespace

from django.test import TestCase

from classification.views.exports.classification_export_formatter_vcf import VCFTargetSystem
from classification.views.exports_grouping.classification_grouping_export_formatter_vcf import (
    VCFRow,
)
from library.utils import ExportTweak


class TestClassificationGroupingExportVCF(TestCase):

    def test_generic_writes_link_and_labs(self):
        export_tweak = ExportTweak(categories={"system": (VCFTargetSystem.CORE, VCFTargetSystem.GENERIC)})
        groupings = [SimpleNamespace(lab="Lab B"), SimpleNamespace(lab="Lab A"), SimpleNamespace(lab="Lab A")]
        row = VCFRow(entry=SimpleNamespace(allele_id=42, classification_groupings=groupings), date_str="20260101")

        info_by_id = {}
        for method in VCFRow.get_export_methods(export_tweak=export_tweak):
            header_id = method.vcf_header.header_id
            if header_id in ("link", "labs"):
                info_by_id[header_id] = method.vcf_header.format_info_value(method(row))

        self.assertEqual(info_by_id["labs"], "labs=Lab_A,Lab_B")
        self.assertRegex(info_by_id["link"], r"^link=\S+/a42\?seen=20260101$")
