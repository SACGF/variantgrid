""" Certify a configured MME node connection.

POSTs a fixed synthetic profile - `"test": true`, no real patient data - to a node's /match
endpoint and prints the response, verifying a connection end to end before any real
classification is submitted.
"""
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from mme.client import MMEClient
from mme.disclaimers import effective_disclaimer

TEST_PROFILE = {
    "id": "vg-test-1",
    "contact": {"name": "VariantGrid MME test", "href": "mailto:test@example.org",
                "institution": "VariantGrid"},
    "species": "NCBITaxon:9606",
    "test": True,
    "genomicFeatures": [{"gene": {"id": "NGLY1"}}],
    "features": [{"id": "HP:0001250", "label": "Seizure", "observed": "yes"}],
}


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("node_id")

    def handle(self, *args, **options):
        node_id = options["node_id"]
        if node_id not in (settings.MME_NODES or {}):
            raise CommandError(f"Unknown MME node '{node_id}'. "
                               f"Configured: {sorted(settings.MME_NODES or {})}")
        data = MMEClient(node_id).match(TEST_PROFILE)
        self.stdout.write(f"{len(data.get('results', []))} result(s) from {node_id}")
        self.stdout.write(str(effective_disclaimer(
            node_id, data.get("disclaimer") or "", data.get("terms") or "")))
