import copy
import dataclasses
import json
import logging
import os
from datetime import date, datetime
from dataclasses import dataclass, field
from typing import Optional, List, Dict

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
    path: str
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


@dataclass_json
@dataclass
class SampleSheetLookup:
    """ Only used as arguments to find existing sample sheets on server - not enough details to create one """
    sequencing_run: str
    hash: str

    @staticmethod
    def from_sample_sheet(sample_sheet: SampleSheet) -> 'SampleSheetLookup':
        return SampleSheetLookup(sequencing_run=sample_sheet.sequencing_run.name, hash=sample_sheet.hash)


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
    sample_sheet_lookup: SampleSheetLookup = field(metadata=config(field_name="sample_sheet"))
    variant_caller: VariantCaller


@dataclass_json
@dataclass
class BamFile:
    path: str
    aligner: Optional[Aligner] = field(default=None, metadata=config(exclude=lambda x: x is None))


@dataclass_json
@dataclass
class VCFFile:
    path: str
    variant_caller: Optional[VariantCaller] = field(default=None, metadata=config(exclude=lambda x: x is None))


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
class SequencingSampleLookup:
    """ Only used as arguments to find existing sequencing sample on server - not enough details to create one """
    sample_sheet_lookup: SampleSheetLookup = field(metadata=config(field_name="sample_sheet"))
    sample_name: str


@dataclass_json
@dataclass
class QC:
    sequencing_sample_lookup: SequencingSampleLookup = field(metadata=config(field_name="sequencing_sample"))
    bam_file: BamFile
    vcf_file: VCFFile


@dataclass_json
@dataclass
class QCGeneList:
    path: str
    qc: QC
    gene_list: List[str]


@dataclass_json
@dataclass
class QCExecStats:
    path: str
    qc: QC
    created: datetime
    modified: datetime
    hash: str
    is_valid: bool
    deduplicated_reads: int
    indels_dbsnp_percent: float
    mean_coverage_across_genes: float
    mean_coverage_across_kit: float
    median_insert: float
    number_indels: int
    number_snps: int
    percent_10x_goi: float
    percent_20x_goi: float
    percent_20x_kit: float
    percent_error_rate: float
    percent_map_to_diff_chr: float
    percent_read_enrichment: float
    percent_reads: float
    percent_softclip: float
    percent_duplication: float
    reads: int
    sample_id_lod: float
    sex_match: str
    snp_dbsnp_percent: float
    ts_to_tv_ratio: float
    uniformity_of_coverage: float
    percent_100x_goi: Optional[float] = None
    percent_100x_kit: Optional[float] = None
    percent_250x_goi: Optional[float] = None
    percent_250x_kit: Optional[float] = None
    percent_500x_goi: Optional[float] = None
    percent_500x_kit: Optional[float] = None


@dataclass_json
@dataclass
class QCGeneCoverage:
    """ We send this up to associate coverage file path with sequencing sample
        Then, later we upload the coverage file (plus path) """
    qc: QC
    path: str


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

    def create_sample_sheet_combined_vcf_file(self, sample_sheet_combined_vcf_file: SampleSheetCombinedVCFFile):
        json_data = sample_sheet_combined_vcf_file.to_dict()
        return self._post("seqauto/api/v1/sample_sheet_combined_vcf_file/",
                          json_data)

    def create_sequencing_data(self, sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]):
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
            "sample_sheet": sample_sheet_lookup.to_dict(),
            "records": records
        }
        return self._post("seqauto/api/v1/sequencing_files/bulk_create",
                          json_data)

    def create_qc_gene_list(self, qc_gene_list: QCGeneList):
        json_data = qc_gene_list.to_dict()
        return self._post("seqauto/api/v1/qc_gene_list/",
                          json_data)


    def create_multiple_qc_gene_lists(self, qc_gene_lists: List[QCGeneList]):
        json_data = {
            "records": [
                qcgl.to_dict() for qcgl in qc_gene_lists
            ]
        }
        return self._post("seqauto/api/v1/qc_gene_list/bulk_create",
                          json_data)

    def create_qc_exec_stats(self, qc_exec_stats: QCExecStats):
        json_data = qc_exec_stats.to_dict()
        return self._post("seqauto/api/v1/qc_exec_summary/",
                          json_data)

    def create_multiple_qc_exec_stats(self, qc_exec_stats: List[QCExecStats]):
        json_data = {
            "records": [
                qces.to_dict() for qces in qc_exec_stats
            ]
        }
        return self._post("seqauto/api/v1/qc_exec_summary/bulk_create",
                          json_data)

    def create_multiple_qc_gene_coverage(self, qc_gene_coverage_list: List[QCGeneCoverage]):
        json_data = {
            "records": [
                qcgc.to_dict() for qcgc in qc_gene_coverage_list
            ]
        }
        return self._post("seqauto/api/v1/qc_gene_coverage/bulk_create",
                          json_data)

    def upload_file(self, filename: str):
        url = self._get_url("upload/api/v1/file_upload")
        with open(filename, "rb") as f:
            kwargs = {
                "files": {"file": f},
                "params": {"path": filename}
            }
            return requests.post(url, headers=self.headers, **kwargs)


class Command(BaseCommand):
    """
        Version 5 - Start replacing hardcoded files with examining files
        Version 6 - Extract into library / example script
    """
    TEST_DATA_DIR = os.path.join(settings.BASE_DIR, 'seqauto', 'test_data')
    HAEM_DIR = os.path.join(TEST_DATA_DIR, "clinical_hg38", "idt_haem",
                            "Haem_20_999_201231_M02027_0112_000000000_JFT79")

    def add_arguments(self, parser):
        parser.add_argument('--token', help="API token (create via drf_user_api_token)")
        parser.add_argument('--server', help="default = obtain from current site config")
        parser.add_argument('--step', help="Step to run (default all)")

    def handle(self, *args, **options):
        api_token = options["token"]
        step = options["step"]
        server = options['server']

        vg_api = VariantGridAPI(server, api_token)

        experiment = "HAEM_20_999"
        enrichment_kit = EnrichmentKit(name='idt_haem', version=1)
        sequencing_run = SequencingRun(path=self.HAEM_DIR,
                                       name="Haem_20_999_201231_M02027_0112_000000000_JFT79",
                                       date="2020-12-31",
                                       sequencer="SN1101",
                                       experiment=experiment,
                                       enrichment_kit=enrichment_kit)

        sequencing_samples = [
            SequencingSample(sample_id="fake_sample_1",
                             sample_project=None,
                             sample_number=1,
                             lane=3,
                             barcode="GCCAAT",
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

        sample_sheet = SampleSheet(
            path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/SampleSheet.csv",
            sequencing_run=sequencing_run,
            file_last_modified=1725941707.0033002,
            hash="f0ac87bcae3f0e56b3f65b70fd6389ce",
            sequencing_samples=sequencing_samples)

        sample_sheet_lookup = SampleSheetLookup.from_sample_sheet(sample_sheet)

        variant_caller_var_dict = VariantCaller(name="VarDict", version="1.8.2")
        sample_sheet_combined_vcf_file = SampleSheetCombinedVCFFile(
            path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/2_variants/Haem_20_999_201231_M02027_0112_000000000_JFT79.vardict.hg38.vcf.gz",
            sample_sheet_lookup=sample_sheet_lookup,
            variant_caller=variant_caller_var_dict)

        aligner = Aligner(name='BWA', version="0.7.18")
        variant_caller_gatk = VariantCaller(name="GATK", version="4.1.9.0")
        bam_file_1 = BamFile(
            path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/1_BAM/fake_sample_1.hg38.bam",
            aligner=aligner)
        vcf_file_1 = VCFFile(
            path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/2_variants/gatk_per_sample/fake_sample_1.gatk.hg38.vcf.gz",
            variant_caller=variant_caller_gatk)

        bam_file_2 = BamFile(
            path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/1_BAM/fake_sample_2.hg38.bam",
            aligner=aligner)
        vcf_file_2 = VCFFile(
            path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/2_variants/gatk_per_sample/fake_sample_2.gatk.hg38.vcf.gz",
            variant_caller=variant_caller_gatk)

        sequencing_files = [
            SequencingFile(sample_name="fake_sample_1",
                           fastq_r1="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/0_fastq/fake_sample_1_R1.fastq.gz",
                           fastq_r2="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/0_fastq/fake_sample_1_R2.fastq.gz",
                           bam_file=bam_file_1,
                           vcf_file=vcf_file_1),
            SequencingFile(sample_name="fake_sample_2",
                           fastq_r1="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/0_fastq/fake_sample_2_R1.fastq.gz",
                           fastq_r2="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/0_fastq/fake_sample_2_R1.fastq.gz",
                           bam_file=bam_file_2,
                           vcf_file=vcf_file_2)
        ]

        qc_by_sample_name = self._get_qc_by_sample_name(sample_sheet_lookup, sequencing_files)
        gene_list = [
            "TUBA1A",
            "TUBA8",
            "FLNA",
            "TUBB2B",
            "TUBB3",
            "COL4A1",
            "KIAA1279"
        ]
        qc_gene_lists = [
            QCGeneList(
                path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/0_goi/Haem_20_999_201231_M02027_0112_000000000_JFT79_fake_sample_1.txt",
                qc=qc_by_sample_name["fake_sample_1"],
                gene_list=gene_list),
            QCGeneList(
                path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/0_goi/Haem_20_999_201231_M02027_0112_000000000_JFT79_fake_sample_2.txt",
                qc=qc_by_sample_name["fake_sample_2"],
                gene_list=gene_list)
        ]

        qc_exec_stats = [
            QCExecStats(
                qc=qc_by_sample_name["fake_sample_1"],
                path='/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/4_QC/exec_stats/fake_sample_1_qc_summary.txt',
                created='2024-10-18T11:45:26.826823+10:30', modified='2024-10-18T11:45:26.826844+10:30', hash='',
                is_valid=True, deduplicated_reads=165107, indels_dbsnp_percent=95.75, mean_coverage_across_genes=162.84,
                mean_coverage_across_kit=201.43, median_insert=153.0, number_indels=923, number_snps=363,
                percent_10x_goi=100.0, percent_20x_goi=100.0, percent_20x_kit=98.35, percent_error_rate=0.91,
                percent_map_to_diff_chr=0.75, percent_read_enrichment=54.28, percent_reads=3.552, percent_softclip=0.02,
                percent_duplication=3.42, reads=120760, sample_id_lod=16.6, sex_match='M=yes', snp_dbsnp_percent=96.62,
                ts_to_tv_ratio=2.1, uniformity_of_coverage=84.69),
            QCExecStats(
                qc=qc_by_sample_name["fake_sample_2"],
                path='/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/4_QC/exec_stats/fake_sample_2_qc_summary.txt',
                created='2024-10-18T11:45:26.838244+10:30', modified='2024-10-18T11:45:26.838262+10:30', hash='',
                is_valid=True, deduplicated_reads=275107, indels_dbsnp_percent=88.75, mean_coverage_across_genes=162.84,
                mean_coverage_across_kit=150.43, median_insert=222.0, number_indels=853, number_snps=1213,
                percent_10x_goi=100.0, percent_20x_goi=100.0, percent_20x_kit=87.35, percent_error_rate=0.55,
                percent_map_to_diff_chr=0.75, percent_read_enrichment=55.28, percent_reads=3.32, percent_softclip=0.02,
                percent_duplication=5.42, reads=410760, sample_id_lod=16.6, sex_match='M=yes', snp_dbsnp_percent=88.62,
                ts_to_tv_ratio=2.1, uniformity_of_coverage=83.69)
        ]

        qc_gene_coverage_list = [
            QCGeneCoverage(qc=qc_by_sample_name["fake_sample_1"],
                           path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/4_QC/bam_stats/samples/fake_sample_1.per_gene_coverage.tsv.gz"),
            QCGeneCoverage(qc=qc_by_sample_name["fake_sample_2"],
                           path="/home/dlawrence/localwork/variantgrid/seqauto/test_data/clinical_hg38/idt_haem/Haem_20_999_201231_M02027_0112_000000000_JFT79/4_QC/bam_stats/samples/fake_sample_2.per_gene_coverage.tsv.gz")
        ]

        gene_coverage_file = os.path.join(self.HAEM_DIR, "4_QC", "bam_stats", "samples",
                                          "fake_sample_1.per_gene_coverage.tsv.gz")
        vcf_filename = os.path.join(self.HAEM_DIR, "2_variants",
                                    "Haem_20_999_201231_M02027_0112_000000000_JFT79.vardict.hg38.vcf.gz")

        API_STEPS = {
            "experiment": lambda: vg_api.create_experiment(experiment),
            "enrichment_kit": lambda: vg_api.create_enrichment_kit(enrichment_kit),
            "sequencing_run": lambda: vg_api.create_sequencing_run(sequencing_run),
            "sample_sheet": lambda: vg_api.create_sample_sheet(sample_sheet),
            "sample_sheet_combined_vcf_file": lambda: vg_api.create_sample_sheet_combined_vcf_file(
                sample_sheet_combined_vcf_file),
            "sequencing_data": lambda: vg_api.create_sequencing_data(sample_sheet_lookup, sequencing_files),
            "qc_gene_list": lambda: vg_api.create_qc_gene_list(qc_gene_lists[0]),
            "qc_gene_lists": lambda: vg_api.create_multiple_qc_gene_lists(qc_gene_lists),
            "qc_exec_summary": lambda: vg_api.create_qc_exec_stats(qc_exec_stats[0]),
            "qc_exec_summaries": lambda: vg_api.create_multiple_qc_exec_stats(qc_exec_stats),
            "qc_gene_coverage": lambda: vg_api.create_multiple_qc_gene_coverage(qc_gene_coverage_list),
            "upload_qc_gene_coverage_file": lambda: vg_api.upload_file(gene_coverage_file),
            "upload_vcf_file": lambda: vg_api.upload_file(vcf_filename),
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
    def _get_qc_by_sample_name(sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]) -> Dict[
        str, QC]:
        bam_and_vcf_by_name = {}
        for sf in sequencing_files:
            bam_file = dataclasses.replace(sf.bam_file, aligner=None)
            vcf_file = dataclasses.replace(sf.vcf_file, variant_caller=None)
            bam_and_vcf_by_name[sf.sample_name] = (bam_file, vcf_file)

        qc_by_name = {}
        for sample_name, (bam_file, vcf_file) in bam_and_vcf_by_name.items():
            sequencing_sample_lookup = SequencingSampleLookup(sample_sheet_lookup=sample_sheet_lookup,
                                                              sample_name=sample_name)
            qc_by_name[sample_name] = QC(sequencing_sample_lookup=sequencing_sample_lookup,
                                         bam_file=bam_file,
                                         vcf_file=vcf_file)
        return qc_by_name
