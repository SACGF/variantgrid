import re

from django.template import Context, Template
from django.test import TestCase


class TestPathogenicityThresholdsTemplateTag(TestCase):
    """ Sanity-check that the {% pathogenicity_thresholds %} tag emits the
    Pejaver-calibrated REVEL/CADD bands so variant_details.html spreads them
    into CODE_THRESHOLDS without further edits. """

    def test_revel_and_cadd_phred_calibrated_bands_present(self):
        rendered = Template(
            "{% load pathogenicity_tags %}{% pathogenicity_thresholds %}"
        ).render(Context({}))
        # Tolerate JSON whitespace and trailing-zero stripping.
        self.assertRegex(rendered, r'"revel_score":\s*\[\s*0\.29(?:0)?\s*,\s*0\.644\s*\]')
        self.assertRegex(rendered, r'"cadd_phred":\s*\[\s*22\.7\s*,\s*25\.3\s*\]')
        # Sanity: VARITY_ER (no calibration) must not appear.
        self.assertNotRegex(rendered, r'"varity_er_score"')
