from analysis.analysis_import_export import analysis_import
from upload.models import UploadedAnalysis
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportAnalysisTask(ImportTask):

    def process_items(self, file_upload):
        analysis = analysis_import(file_upload.user, genome_build=None, filename=file_upload.get_filename())
        UploadedAnalysis.objects.create(file_upload=file_upload, analysis=analysis)
        return 1


ImportAnalysisTask = app.register_task(ImportAnalysisTask())  # @UndefinedVariable
