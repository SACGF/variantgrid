import logging
import re

from genes.canonical_transcripts.canonical_transcript_manager import CanonicalTranscriptManager
from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneCoverageCollection
from library.file_utils import name_from_filename
from library.log_utils import log_traceback
import pandas as pd
from seqauto.models_enums import DataState
from snpdb.models import Sample, GenomeBuild
from upload.models import UploadedGeneCoverage
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportGeneCoverageTask(ImportTask):
    EXPECTED_SUFFIX = "_genes.txt"

    @staticmethod
    def can_process_file(filename):
        # TAU coverage files only at the moment
        EXPECTED_COLUMNS = ['Gene/Exon', 'Min Coverage', 'Mean Coverage',
                            'Standard Deviation Coverage', '% bases >20x', '% bases <10x',
                            '% bases @ 0x', 'Estimated Sensitivity of Detection']

        if filename.endswith(ImportGeneCoverageTask.EXPECTED_SUFFIX):
            try:
                df = pd.read_csv(filename, sep='\t', nrows=1)
                return list(df.columns) == EXPECTED_COLUMNS
            except:
                pass
        return False

    @staticmethod
    def match_sample(user, filename):
        sample = None
        # Need write permission as it's a 1-to-1 field so only 1 can be assigned
        samples_qs = Sample.filter_for_user(user, has_write_permission=True)
        samples_qs = samples_qs.filter(uploadedgenecoverage__isnull=True)
        try:
            name = name_from_filename(filename)
            pattern = "(.*)%s" % name_from_filename(ImportGeneCoverageTask.EXPECTED_SUFFIX)
            logging.info("name: %s, pattern: %s", name, pattern)
            if m := re.match(pattern, name):
                sample_name = m.group(1)
                logging.info("Looking for sample named: %s", sample_name)
                sample = samples_qs.get(name=sample_name)
        except:
            log_traceback()
        return sample

    def process_items(self, uploaded_file):
        filename = uploaded_file.get_filename()

        try:
            uploaded_gene_coverage = UploadedGeneCoverage.objects.get(uploaded_file=uploaded_file)
        except UploadedGeneCoverage.DoesNotExist:
            uploaded_gene_coverage = UploadedGeneCoverage(uploaded_file=uploaded_file)

        uploaded_gene_coverage.sample = self.match_sample(uploaded_file.user, filename)
        uploaded_gene_coverage.save()
        gene_coverage_collection = GeneCoverageCollection.objects.create(path=filename,
                                                                         data_state=DataState.RUNNING,
                                                                         genome_build=GenomeBuild.legacy_build())
        uploaded_gene_coverage.gene_coverage_collection = gene_coverage_collection
        uploaded_gene_coverage.save()

        gene_matcher = GeneSymbolMatcher()
        canonical_transcript_manager = CanonicalTranscriptManager()
        gene_coverage_collection.load_from_file(None, gene_matcher=gene_matcher, canonical_transcript_manager=canonical_transcript_manager)
        gene_coverage_collection.data_state = DataState.COMPLETE
        gene_coverage_collection.save()

        num_records = gene_coverage_collection.genecoverage_set.count()
        return num_records


ImportGeneCoverageTask = app.register_task(ImportGeneCoverageTask())  # @UndefinedVariable
