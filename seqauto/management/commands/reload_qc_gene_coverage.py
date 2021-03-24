"""
    GeneCoverage used to only store genes that matched

    Now we store everything, so we can get Gold standard GeneCoverage even for
    genes where the canonical transcript is not matched across refseq/ensembl

    So - find the old stuff and reload it.
"""

from django.core.management.base import BaseCommand

from genes.models import GeneCoverageCollection
from genes.tasks.gene_coverage_tasks import reload_gene_coverage_collection
from library.log_utils import log_traceback
from seqauto.models import SequencingRun, QCGeneCoverage
from snpdb.models import DataState


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', required=False)
        parser.add_argument('--gold', action='store_true', required=False)

    def handle(self, *args, **options):
        do_all = options["all"]
        gold = options["gold"]

        sequencing_run_kwargs = {}
        if gold:
            sequencing_run_kwargs["gold_standard"] = True
        for sr in SequencingRun.objects.filter(**sequencing_run_kwargs):
            sheet = sr.sequencingruncurrentsamplesheet.sample_sheet

            for ss in sheet.sequencingsample_set.filter(is_control=False):
                try:
                    qc = ss.get_single_qc()
                    try:
                        qc_gene_coverage = qc.qcgenecoverage
                    except:
                        other_qc_path = QCGeneCoverage.get_path_from_qc(qc)
                        gcc = GeneCoverageCollection.objects.create(path=other_qc_path,
                                                                    data_state=DataState.COMPLETE,
                                                                    genome_build=qc.genome_build)

                        qc_gene_coverage = QCGeneCoverage.objects.create(path=other_qc_path,
                                                                         qc=qc,
                                                                         data_state=DataState.COMPLETE,
                                                                         gene_coverage_collection=gcc)

                    gene_coverage_collection = qc_gene_coverage.gene_coverage_collection
                    has_unmatched_genes = gene_coverage_collection.genecoverage_set.filter(gene__isnull=True).exists()

                    if do_all or not has_unmatched_genes:
                        print(f"Need to reload gene coverage for {ss}")
                        task = reload_gene_coverage_collection.si(gene_coverage_collection.pk)  # @UndefinedVariable
                        task.apply_async()
                    else:
                        print(f"{ss} is ok")

                except:
                    log_traceback()
