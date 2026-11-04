""" SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS: a per-kit include list and an exclude list choosing which URL a run gets. """
import datetime

from django.test import SimpleTestCase, override_settings

from seqauto.models import SequencingRun

RAVEN_BATCH = {
    "label": "RAVEN",
    "url_pattern": "http://raven/batchID/%(sequencing_run)s",
    "exclude_enrichment_kits": ["wgs"],
}
RAVEN_WGS = {
    "label": "RAVEN",
    "url_pattern": "http://raven/wgs/run/%(sequencing_run)s",
    "enrichment_kits": ["WGS"],
}


@override_settings(SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS=[RAVEN_BATCH, RAVEN_WGS])
class SequencingRunExternalLinksTest(SimpleTestCase):
    def test_kit_picks_url(self):
        run_date = datetime.date(2026, 10, 30)
        panel_name = "261030_LH00639_0207_B23N5YFLT3"
        wgs_name = "WGS_26_081_20261030_LH00639_0207_B23N5YFLT3"

        def urls_for(name, kit):
            return [url for _, url in SequencingRun.get_external_links_for(name, run_date, kit)]

        self.assertEqual([f"http://raven/batchID/{panel_name}"], urls_for(panel_name, "idt_exome"))
        self.assertEqual([f"http://raven/batchID/{panel_name}"], urls_for(panel_name, None))
        self.assertEqual([f"http://raven/wgs/run/{wgs_name}"], urls_for(wgs_name, "wgs"))
