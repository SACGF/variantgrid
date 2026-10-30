import logging
import os
import re

import pandas as pd

from genes.canonical_transcripts.canonical_transcript_manager import CanonicalTranscriptManager
from genes.enums import AnnotationConsortium
from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneCoverageCollection, TranscriptVersion
from genes.models.models_gene_coverage import GeneCoverage
from library.django_utils.django_queryset_sql_transformer import get_queryset_with_transformer_hook
from library.log_utils import log_traceback
from library.utils.file_utils import name_from_filename
from snpdb.enums import DataState
from snpdb.models import GenomeBuild, Sample
from upload.models import UploadedGeneCoverage
from upload.signals.signals import gene_coverage_upload_signal
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportGeneCoverageTask(ImportTask):
    EXPECTED_SUFFIX = ".per_gene_coverage.tsv.gz"

    @staticmethod
    def can_process_file(filename):
        # the filename can get mangled if there are duplicate names uploaded, so we can't check for exact suffix
        # TAU coverage files only at the moment
        EXPECTED_COLUMNS = ['Size', 'Mean', 'StdDev', 'Max', 'Min', '% >=1X', '% >=2X', '% >=5X']
        if "per_gene_coverage.tsv" in filename:
            try:
                df = pd.read_csv(filename, sep='\t', nrows=1)
                return all([c in df.columns for c in EXPECTED_COLUMNS])
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

    def process_items(self, file_upload):
        # This was originally written for RUNX1
        genome_build = GenomeBuild.legacy_build()
        enrichment_kit = None

        # A file sent for a sequencing record (seqauto's QCGeneCoverage) takes its build and kit from there,
        # so that we can get the right canonical transcripts etc
        source = None
        for _, response in gene_coverage_upload_signal.send(sender=os.path.basename(__file__), file_upload=file_upload):
            if response:
                source = response
                genome_build = source.genome_build
                enrichment_kit = source.enrichment_kit

        filename = file_upload.get_filename()

        try:
            uploaded_gene_coverage = UploadedGeneCoverage.objects.get(file_upload=file_upload)
        except UploadedGeneCoverage.DoesNotExist:
            uploaded_gene_coverage = UploadedGeneCoverage(file_upload=file_upload)

        # Linking this way was done for RUNX1
        auto_match_samples = False
        if not source and auto_match_samples:
            uploaded_gene_coverage.sample = self.match_sample(file_upload.user, filename)
            uploaded_gene_coverage.save()

        # TODO: This is hardcoded - need to be able to pass this in or select in GUI
        annotation_consortium = AnnotationConsortium.REFSEQ
        gene_coverage_collection = GeneCoverageCollection.objects.create(path=filename,
                                                                         data_state=DataState.RUNNING,
                                                                         genome_build=genome_build)
        uploaded_gene_coverage.gene_coverage_collection = gene_coverage_collection
        uploaded_gene_coverage.save()

        if source:
            source.link_collection(gene_coverage_collection)

        gene_matcher = GeneSymbolMatcher()
        canonical_transcript_manager = CanonicalTranscriptManager()
        transcript_versions_by_id = TranscriptVersion.transcript_versions_by_id(genome_build, annotation_consortium)
        gene_coverage_collection.load_from_file(enrichment_kit, gene_matcher=gene_matcher,
                                                canonical_transcript_manager=canonical_transcript_manager,
                                                transcript_versions_by_id=transcript_versions_by_id)
        gene_coverage_collection.data_state = DataState.COMPLETE
        gene_coverage_collection.save()

        # Count from the collection's own partition: a query on the inheritance parent locks every child table
        qs = get_queryset_with_transformer_hook(klass=GeneCoverage)
        qs.add_sql_transformer(gene_coverage_collection.sql_partition_transformer)
        num_records = qs.filter(gene_coverage_collection=gene_coverage_collection).count()
        return num_records


ImportGeneCoverageTask = app.register_task(ImportGeneCoverageTask())  # @UndefinedVariable
