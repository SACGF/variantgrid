import json
import logging
import os
from datetime import date, datetime
from dataclasses import dataclass, field
from typing import Optional, List

import requests
from dataclasses_json import dataclass_json, config
from django.conf import settings
from django.core.management.base import BaseCommand

@dataclass_json
@dataclass
class EnrichmentKit:
    name: str
    version: int

@dataclass_json
@dataclass
class SequencingRun:
    name: str
    date: date
    sequencer: str
    experiment: str
    enrichment_kit: EnrichmentKit


@dataclass_json
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
    data: List[dict] = field(default_factory=lambda: [], metadata=config(field_name="sequencingsampledata_set"))


@dataclass_json
@dataclass
class SampleSheet:
    path: str
    sequencing_run: SequencingRun
    file_last_modified: datetime
    hash: str
    sequencing_samples: List[SequencingSample] = field(metadata=config(field_name="sequencingsample_set"))

    def to_simplified_dict(self):
        return {
            "sequencing_run": self.sequencing_run.name,
            "hash": self.hash,
        }


@dataclass_json
@dataclass
class Aligner:
    name: str
    version: str

@dataclass_json
@dataclass
class VariantCaller:
    name: str
    version: str
    run_params: Optional[str] = None


@dataclass_json
@dataclass
class SampleSheetCombinedVCFFile:
    path: str
    sample_sheet: SampleSheet
    variant_caller: VariantCaller


@dataclass_json
@dataclass
class BamFile:
    path: str
    aligner: Aligner


@dataclass_json
@dataclass
class VCFFile:
    path: str
    variant_caller: VariantCaller


@dataclass_json
@dataclass
class SequencingFile:
    sample_name: str
    fastq_r1: str
    fastq_r2: str
    bam_file: BamFile
    vcf_file: VCFFile


@dataclass_json
@dataclass
class QC:
    sequencing_sample: SequencingSample
    bam_file: BamFile
    vcf_file: VCFFile

@dataclass
class SampleGeneList:
    """ This is never sent to API (it needs QC)
        But it can be convenient to create it in this object, and create QCs later """
    path: str
    sample_name: str
    gene_list: List[str]


@dataclass_json
@dataclass
class QCGeneList:
    path: str
    qc: QC
    gene_list: List[str]


@dataclass_json
@dataclass
class QCExecStats:
    pass

@dataclass_json
@dataclass
class QCGeneCoverage:
    pass



class VariantGridAPI:
    def __init__(self, server, api_token):
        self.server = server
        self.headers = {"Authorization": f"Token {api_token}"}

    def _get_url(self, url):
        return os.path.join(self.server, url)

    def _post(self, path, json_data):
        url = self._get_url(path)
        logging.info("url='%s', JSON data:", url)
        logging.info(json.dumps(json_data))
        response = requests.post(url, headers=self.headers, json=json_data)
        if not response.ok:
            response.raise_for_status()
        try:
            json_response = response.json()
        except Exception as e:
            json_response = f"Couldn't convert JSON: {e}"
        return json_response

    def create_experiment(self, experiment: str):
        json_data = {
            "name": experiment
        }
        return self._post("seqauto/api/v1/experiment/", json_data)

    def create_enrichment_kit(self, enrichment_kit: EnrichmentKit):
        return self._post("seqauto/api/v1/enrichment_kit/",
                          enrichment_kit.to_dict())

    def create_sequencing_run(self, sequencing_run: SequencingRun):
        return self._post("seqauto/api/v1/sequencing_run/",
                          sequencing_run.to_dict())

    def create_sample_sheet(self, sample_sheet: SampleSheet):
        json_data = sample_sheet.to_dict()
        # We don't want all sequencing_run just the name
        sequencing_run = json_data.pop("sequencing_run")
        json_data["sequencing_run"] = sequencing_run["name"]
        return self._post("seqauto/api/v1/sample_sheet/",
                          json_data)

    def create_sample_sheet_combined_vcf_file(self, sample_sheet_combined_vcf_file):
        json_data = sample_sheet_combined_vcf_file.to_dict()
        json_data["sample_sheet"] = sample_sheet_combined_vcf_file.sample_sheet.to_simplified_dict()
        return self._post("seqauto/api/v1/sample_sheet_combined_vcf_file/",
                          json_data)

    def create_sequencing_data(self, sample_sheet: SampleSheet, sequencing_files: List[SequencingFile]):
        records = []
        for sf in sequencing_files:
            data = sf.to_dict()
            # put into hierarchial JSON DRF expects
            fastq_r1 = data.pop("fastq_r1")
            fastq_r2 = data.pop("fastq_r2")
            data["unaligned_reads"] = {
                "fastq_r1": {"path": fastq_r1},
                "fastq_r2": {"path": fastq_r2}
            }
            records.append(data)

        json_data = {
            "sample_sheet": sample_sheet.to_simplified_dict(),
            "records": records
        }
        return self._post("seqauto/api/v1/sequencing_files/bulk_create",
                          json_data)

    def create_qc_gene_list(self):
        # TODO
        pass

    def create_multiple_qc_gene_lists(self, qc_gene_lists: List[QCGeneList]):
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
    def add_arguments(self, parser):
        parser.add_argument('--token', help="API token (create via drf_user_api_token)")
        parser.add_argument('--server', help="default = obtain from current site config")
        parser.add_argument('--step', help="Step to run (default all)")

    def handle(self, *args, **options):
        api_token = options["token"]
        step = options["step"]
        server = options['server']

        vg_api = VariantGridAPI(server, api_token)

        experiment = "HAEM_21_001"
        enrichment_kit = EnrichmentKit(name='idt_haem', version=1)
        sequencing_run = SequencingRun(name="Haem_21_001_210216_M02027_0112_000000000_JFT79",
                                       date="2021-02-16",
                                       sequencer="SN1101",
                                       experiment=experiment,
                                       enrichment_kit=enrichment_kit)

        sequencing_samples = [
            SequencingSample(sample_id="fake_sample_1",
                             sample_project=None,
                             sample_number=1,
                             lane= 3,
                             barcode= "GCCAAT",
                             enrichment_kit=enrichment_kit,
                             is_control=False,
                             failed=False,
                             data=[
                                 {
                                    "column": "SAPOrderNumber",
                                    "value": "SAP1000001"
                                 }
                            ]),
            SequencingSample(sample_id="fake_sample_2",
                             sample_project=None,
                             sample_number=2,
                             lane=3,
                             barcode="CAGATC",
                             enrichment_kit=enrichment_kit,
                             is_control=False,
                             failed=False,
                             data=[
                                 {
                                     "column": "SAPOrderNumber",
                                     "value": "SAP1000002"
                                 }
                             ])
        ]

        sample_sheet = SampleSheet(path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/SampleSheet.csv",
                                   sequencing_run=sequencing_run,
                                   file_last_modified=1725941707.0033002,
                                   hash="f0ac87bcae3f0e56b3f65b70fd6389ce",
                                   sequencing_samples=sequencing_samples)

        variant_caller_var_dict = VariantCaller(name="VarDict", version="1.8.2")
        sample_sheet_combined_vcf_file = SampleSheetCombinedVCFFile(path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/2_variants/Haem_21_001_210216_M02027_0112_000000000_JFT79.vardict.hg38.vcf.gz",
                                                                    sample_sheet=sample_sheet,
                                                                    variant_caller=variant_caller_var_dict)

        aligner = Aligner(name='BWA', version="0.7.18")
        variant_caller_gatk = VariantCaller(name="GATK", version="4.1.9.0")
        bam_file_1 = BamFile(path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/1_BAM/fake_sample_1.hg38.bam",
                             aligner=aligner)
        vcf_file_1 = VCFFile(path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/2_variants/gatk_per_sample/fake_sample_1.gatk.hg38.vcf.gz",
                             variant_caller=variant_caller_gatk)

        bam_file_2 = BamFile(path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/1_BAM/fake_sample_2.hg38.bam",
                             aligner=aligner)
        vcf_file_2 = VCFFile(path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/2_variants/gatk_per_sample/fake_sample_2.gatk.hg38.vcf.gz",
                             variant_caller=variant_caller_gatk)

        sequencing_files = [
            SequencingFile(sample_name="fake_sample_1",
                           fastq_r1="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/0_fastq/fake_sample_1_R1.fastq.gz",
                           fastq_r2="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/0_fastq/fake_sample_1_R2.fastq.gz",
                           bam_file=bam_file_1,
                           vcf_file=vcf_file_1),
            SequencingFile(sample_name="fake_sample_2",
                           fastq_r1="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/0_fastq/fake_sample_2_R1.fastq.gz",
                           fastq_r2="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_21_001_210216_M02027_0112_000000000_JFT79/0_fastq/fake_sample_2_R1.fastq.gz",
                           bam_file=bam_file_2,
                           vcf_file=vcf_file_2)
        ]

        gene_list = [
            "TUBA1A",
            "TUBA8",
            "FLNA",
            "TUBB2B",
            "TUBB3",
            "COL4A1",
            "KIAA1279"
        ]
        sample_gene_lists = [
            SampleGeneList(sample_name="fake_sample_1", gene_list=gene_list),
            SampleGeneList(sample_name="fake_sample_2", gene_list=gene_list)
        ]
        qc_gene_lists = self._convert_sample_to_qc_gene_lists(sample_sheet, sequencing_files, sample_gene_lists)


        vcf_filename = os.path.join(self.TEST_DATA_DIR, "clinical_hg38", "idt_haem",
                                    "Haem_21_001_210216_M02027_0112_000000000_JFT79", "2_variants",
                                    "Haem_21_001_210216_M02027_0112_000000000_JFT79.vardict.hg38.vcf.gz")

        API_STEPS = {
            "experiment": lambda: vg_api.create_experiment(experiment),
            "enrichment_kit": lambda: vg_api.create_enrichment_kit(enrichment_kit),
            "sequencing_run": lambda: vg_api.create_sequencing_run(sequencing_run),
            "sample_sheet": lambda: vg_api.create_sample_sheet(sample_sheet),
            "sample_sheet_combined_vcf_file": lambda: vg_api.create_sample_sheet_combined_vcf_file(sample_sheet_combined_vcf_file),
            "sequencing_data": lambda: vg_api.create_sequencing_data(sample_sheet, sequencing_files),
            "qc_gene_lists": lambda: vg_api.create_multiple_qc_gene_lists(qc_gene_lists),
            "qc_exec_summaries": lambda: None,
            "qc_gene_coverage": lambda: None,
            "upload_qc_gene_coverage_file": lambda: None,
            "upload_file": lambda: logging.info("TODO: qc_gene_lists"),
        }

        for name, func in API_STEPS.items():
            if step:
                if name != step:
                    continue
            print(f"{name=}")
            result = func()
            print(f"{result=}")
            print("-" * 50)

    @staticmethod
    def _convert_sample_to_qc_gene_lists(sample_sheet: SampleSheet, sequencing_files: List[SequencingFile],
                                         sample_gene_lists: List[SampleGeneList]):
        # Make a lookup
        # Map into QC objects
        return []



