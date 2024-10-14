import json
import os

import requests
from django.conf import settings
from django.core.management.base import BaseCommand
from django.urls import reverse

from library.django_utils import get_url_from_view_path


class Command(BaseCommand):
    """
        Version 1 - Hardcode JSON to recreate what is made via 'scan_run_jobs', just for Haem
        Version 2 - Also send up VCF and properly link it
        Version 3 - Hardcode JSON to recreate exome as well
        Version 4 - optimise, send up minimal JSON
        Version 5 - Start replacing hardcoded files with examining files
        Version 6 - Extract into library / example script
    """
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, 'seqauto', 'test_data')
    API_DATA = os.path.join(TEST_DATA_DIR, 'api_client')
    URLS = {
        "experiment": "api_experiment",
        "enrichment_kit": "api_enrichment_kit",
        "sequencing_run": "api_sequencing_run",
        "sample_sheet": "api_sample_sheet",
        "sample_sheet_combined_vcf_file": "api_sample_sheet_combined_vcf_file",
        "qc_gene_list": "api_qc_gene_list",
        "upload_file": "api_file_upload",
    }

    SEQUENCING_RUNS = {
        "Haem_21_001_210216_M02027_0112_000000000_JFT79": {
            "experiment": os.path.join(API_DATA, "experiment", "haem_21_001.json"),
            "enrichment_kit": os.path.join(API_DATA, "enrichment_kit", "idt_haem.json"),
            "sequencing_run": os.path.join(API_DATA, "sequencing_run", "haem_21_minimal.json"),
            "sample_sheet": os.path.join(API_DATA, "sample_sheet", "haem_21_minimal.json"),
            "sample_sheet_combined_vcf_file": os.path.join(API_DATA, "sample_sheet_combined_vcf_file", "haem_21_sample_1.json"),
            "qc_gene_list": os.path.join(API_DATA, "qc_gene_list", "haem_21_sample_1.json"),
            "upload_file": os.path.join(TEST_DATA_DIR, "clinical_hg38", "idt_haem", "Haem_21_001_210216_M02027_0112_000000000_JFT79", "2_variants", "Haem_21_001_210216_M02027_0112_000000000_JFT79.vardict.hg38.vcf.gz"),
        }
    }

    def add_arguments(self, parser):
        parser.add_argument('--token', help="API token (create via drf_user_api_token)")
        parser.add_argument('--server', help="default = obtain from current site config")
        parser.add_argument('--step', help="Step to run (default all)")
        parser.add_argument('--sequencing_run', help="To just run 1 (default all)")

    def handle(self, *args, **options):
        step = options["step"]
        server = options['server']
        sequencing_run = options['sequencing_run']
        if sequencing_run is None:
            sequencing_run = list(self.SEQUENCING_RUNS.keys())

        token = options['token']
        headers = {"Authorization": f"Token {token}"}

        for sr in sequencing_run:
            for name, filename in self.SEQUENCING_RUNS[sr].items():
                if step:
                    if name != step:
                        continue
                print(f"{name=}")

                view_path = self.URLS[name]
                kwargs = {}
                if name == "upload_file":
                    kwargs["files"] = {"file": open(filename, 'rb')}
                    kwargs["params"] = {"path": filename}
                else:
                    with open(filename, "r") as f:
                        json_data = json.load(f)
                        kwargs["json"] = json_data
                    view_path += "-list"

                if server:
                    url = server + reverse(view_path)
                else:
                    url = get_url_from_view_path(view_path)
                print(f"{url=}")
                response = requests.post(url, headers=headers, **kwargs)
                if response.ok:
                    json_response = response.json()
                else:
                    json_response = None
                print(f"{response.status_code=} - {json_response=}")
                print("-" * 50)