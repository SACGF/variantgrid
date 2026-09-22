"""
Import of a DRAGEN TSO 500 MetricsOutput.tsv - a run's per-library QC, with a column for every
sample of every pair on the run.

Entry point is ImportDragenTSO500MetricsOutputTask, a single-shot ImportTask (the file has no
variants and no coordinates, so there is no VCF pipeline) returning the number of LibraryQC rows
written. The file is read by upload.tso500.dragen_metrics_output_parser and written by
upload.tso500.dragen_metrics_output_records.

The run is the upload's 'sequencing_run' metadata and is part of the row key, so a file that names
none cannot be keyed and fails the import saying so.

A file that lands before its specimens have been accessioned still writes its rows, with the
extraction claims parked - reconcile_pending_extractions is fired afterwards so a claim that only
became resolvable with this run's CombinedVariantOutput attaches straight away.
"""
from patients.tasks.extraction_matching_tasks import reconcile_pending_extractions
from upload.models import UploadedDragenTSO500MetricsOutput
from upload.tasks.import_task import ImportTask
from upload.tso500.dragen_metrics_output_parser import read_library_qc, read_sections
from upload.tso500.dragen_metrics_output_records import (
    metrics_measured_date,
    metrics_method,
    write_library_qc,
)
from upload.upload_metadata import SEQUENCING_RUN, get_metadata_sequencing_run_name
from variantgrid.celery import app


class ImportDragenTSO500MetricsOutputTask(ImportTask):

    def process_items(self, file_upload):
        sequencing_run_name = get_metadata_sequencing_run_name(file_upload)
        if not sequencing_run_name:
            raise ValueError(f"MetricsOutput.tsv is written per sequencing run and names none itself - "
                             f"upload it with '{SEQUENCING_RUN}' metadata naming the run it came off")
        UploadedDragenTSO500MetricsOutput.objects.get_or_create(file_upload=file_upload)

        sections = read_sections(file_upload.get_filename())
        rows = write_library_qc(read_library_qc(sections), file_upload.user,
                                sequencing_run_name=sequencing_run_name,
                                method=metrics_method(sections),
                                date=metrics_measured_date(sections),
                                file_upload=file_upload)
        reconcile_pending_extractions.delay()
        return len(rows)


ImportDragenTSO500MetricsOutputTask = app.register_task(ImportDragenTSO500MetricsOutputTask())  # @UndefinedVariable
