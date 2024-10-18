import json
import os
from datetime import date, datetime
from dataclasses import dataclass
from typing import Optional, List

import requests
from django.conf import settings
from django.core.management.base import BaseCommand
from django.urls import reverse

from library.django_utils import get_url_from_view_path

@dataclass
class EnrichmentKit:
    name: str
    version: int


@dataclass
class SequencingRun:
    name: str
    date: date
    sequencer: str
    experiment: str
    enrichment_kit: EnrichmentKit


@dataclass
class SequencingSample:
    sample_id: str
    sample_number: int
    lane: int
    barcode: str
    enrichment_kit: EnrichmentKit
    sample_project: Optional[str] = None
    is_control: bool = False
    failed: bool = False
    data = List[dict]

@dataclass
class SampleSheet:
    path: str
    sequencing_run: SequencingRun
    file_last_modified: datetime
    hash: str
    sequencing_samples: List[SequencingSample]

@dataclass
class Aligner:
    name: str
    version: str

@dataclass
class VariantCaller:
    name: str
    version: str
    run_params: Optional[str] = None

@dataclass
class SampleSheetCombinedVCFFile:
    path: str
    sample_sheet: SampleSheet
    variant_caller: VariantCaller

@dataclass
class BamFile:
    path: str
    aligner: Aligner


@dataclass
class VCFFile:
    path: str
    variant_caller: VariantCaller


@dataclass
class SequencingFiles:
    sample_name: str
    fastq_r1: str
    fastq_r2: str
    bam_file: BamFile
    vcf_file: VCFFile

@dataclass
class QCData:
    sample_name: str
    bam_file: BamFile
    vcf_file: VCFFile
    gene_list: List[str]

@dataclass
class QCExecStats:
    pass

@dataclass
class QCGeneCoverage:
    pass

# TODO: Look into https://github.com/lidatong/dataclasses-json


class VGAPI:
    def __init__(self, server, api_token):
        self.server = server
        self.headers = {"Authorization": f"Token {api_token}"}

    def _get_url(self, url):
        return os.path.join(self.server, url)

    def create_experiment(self, experiment: str):
        url = self._get_url("seqauto/api/v1/experiment/")

    def create_enrichment_kit(self, enrichment_kit: EnrichmentKit):
        url = self._get_url("seqauto/api/v1/enrichment_kit/")

    def create_sequencing_run(self, sequencing_run: SequencingRun):
        url = self._get_url("seqauto/api/v1/sequencing_run/")

    def create_sample_sheet(self, sample_sheet: SampleSheet):
        url = self._get_url("seqauto/api/v1/sample_sheet/")

    def create_sample_sheet_combined_vcf_file(self, sample_sheet_combined_vcf_file):
        url = self._get_url("seqauto/api/v1/sample_sheet_combined_vcf_file/")

    def create_sequencing_data(self, sample_sheet: SampleSheet, sequencing_files: List[SequencingFiles]):
        url = self._get_url("/seqauto/api/v1/sequencing_files/bulk_create")

    def create_qc_gene_list(self):
        # TODO
        pass

    def create_multiple_qc_gene_lists(self, sample_sheet: SampleSheet, qc_gene_lists: List[QCData]):
        url = self._get_url("seqauto/api/v1/qc_gene_list/bulk_create")

    def create_qc_exec_stats(self):
        # TODO
        pass

    def create_multiple_qc_exec_stats(self, sample_sheet: SampleSheet, qc_exec_stats: List[QCExecStats]):
        pass

    def create_multiple_qc_gene_coverage(self, sample_sheet: SampleSheet, qc_fene_lists: List[QCGeneCoverage]):
        pass

    def upload_vcf_file(self, vcf_filename):
        url = self._get_url("upload/api/v1/file_upload")
        with open(vcf_filename, "rb") as f:
            kwargs = {
                "files": {"file": f},
                "params": {"path": vcf_filename}
            }
            response = requests.post(url, headers=self.headers, **kwargs)


class Command(BaseCommand):
    """
        Version 5 - Start replacing hardcoded files with examining files
        Version 6 - Extract into library / example script
    """
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, 'seqauto', 'test_data')
    API_DATA = os.path.join(TEST_DATA_DIR, 'api_client')
    URLS = {
        "experiment": "api_experiment-list",
        "enrichment_kit": "api_enrichment_kit-list",
        "sequencing_run": "api_sequencing_run-list",
        "sample_sheet": "api_sample_sheet-list",
        "sample_sheet_combined_vcf_file": "api_sample_sheet_combined_vcf_file-list",
        "qc_gene_lists": "api_qc_gene_list_bulk_create",
        "sequencing_data": "api_sequencing_files_bulk_create",
        "upload_file": "api_file_upload",
    }

    SEQUENCING_RUNS = {
        "Haem_21_001_210216_M02027_0112_000000000_JFT79": {
            "experiment": os.path.join(API_DATA, "experiment", "haem_21_001.json"),
            "enrichment_kit": os.path.join(API_DATA, "enrichment_kit", "idt_haem.json"),
            "sequencing_run": os.path.join(API_DATA, "sequencing_run", "haem_21_minimal.json"),
            "sample_sheet": os.path.join(API_DATA, "sample_sheet", "haem_21_minimal.json"),
            "sample_sheet_combined_vcf_file": os.path.join(API_DATA, "sample_sheet_combined_vcf_file", "haem_21_bulk_exec_summaries.json"),
            "sequencing_data": os.path.join(API_DATA, "sequencing_data", "haem_21_bulk_sequencing_files.json"),
            "qc_gene_lists": os.path.join(API_DATA, "qc_gene_list", "haem_21_bulk_samples.json"),
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

                if server:
                    url = server + reverse(view_path)
                else:
                    url = get_url_from_view_path(view_path)
                print(f"{url=}")
                response = requests.post(url, headers=headers, **kwargs)
                try:
                    json_response = response.json()
                except Exception as e:
                    json_response = f"Couldn't convert JSON: {e}"
                print(f"{response.status_code=} - {json_response=}")
                print("-" * 50)