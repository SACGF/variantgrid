from django.test import SimpleTestCase

from classification.views.exports.classification_export_formatter_redcap import redcap_choices


class RedcapChoicesTest(SimpleTestCase):

    def test_pipe_in_label_does_not_split_choice(self):
        options = [
            {"key": "a", "index": 1, "label": "Pathogenic | Likely pathogenic"},
            {"key": "not_indexed"},
            {"key": "likely_benign", "index": 2},
        ]
        self.assertEqual(redcap_choices(options), "1, Pathogenic / Likely pathogenic | 2, Likely benign")
