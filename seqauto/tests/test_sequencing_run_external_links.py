""" SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS: which configured links apply to a run (min_date, per-kit include
    and exclude lists) and how the URL pattern is filled in. """
import datetime

from django.test import SimpleTestCase, override_settings

from seqauto.models import SequencingRun

RAVEN_BATCH = {
    "label": "RAVEN",
    "url_pattern": "http://raven/batchID/%(sequencing_run)s",
    "min_date": "2026-01-01",
    "exclude_enrichment_kits": ["wgs"],
}
RAVEN_WGS = {
    "label": "RAVEN",
    "url_pattern": "http://raven/wgs/run/%(sequencing_run)s",
    "enrichment_kits": ["WGS"],
}


@override_settings(SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS=[RAVEN_BATCH, RAVEN_WGS])
class SequencingRunExternalLinksTest(SimpleTestCase):
    RUN_DATE = datetime.date(2026, 10, 30)

    def test_kit_picks_url(self):
        panel = SequencingRun.get_external_links_for("261030_LH00639_0207_B23N5YFLT3", self.RUN_DATE, "idt_exome")
        self.assertEqual([("RAVEN", "http://raven/batchID/261030_LH00639_0207_B23N5YFLT3")], panel)

        wgs_name = "WGS_26_081_20261030_LH00639_0207_B23N5YFLT3"
        wgs = SequencingRun.get_external_links_for(wgs_name, self.RUN_DATE, "wgs")
        self.assertEqual([("RAVEN", f"http://raven/wgs/run/{wgs_name}")], wgs)

    def test_include_list_needs_a_kit(self):
        links = SequencingRun.get_external_links_for("261030_LH00639_0207_B23N5YFLT3", self.RUN_DATE, None)
        self.assertEqual(["http://raven/batchID/261030_LH00639_0207_B23N5YFLT3"], [url for _, url in links])

    def test_min_date(self):
        self.assertEqual([], SequencingRun.get_external_links_for("251030_LH00639_0207_B23N5YFLT3",
                                                                  datetime.date(2025, 10, 30), "idt_exome"))
        self.assertEqual([], SequencingRun.get_external_links_for("251030_LH00639_0207_B23N5YFLT3",
                                                                  None, "idt_exome"))
