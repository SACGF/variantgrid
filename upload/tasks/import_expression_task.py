import logging

from expression.cuffdiff.import_cuffdiff import import_cuffdiff
from upload.models import UploadedExpressionFile
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportExpressionTask(ImportTask):
    def process_items(self, uploaded_file):
        logging.debug("ImportExpressionTask: process items")

        if uploaded_file.file_type == UploadedExpressionFile.CUFFDIFF:
            cuff_diff_file = import_cuffdiff(uploaded_file.get_file(), uploaded_file.name, uploaded_file.user)
            uploaded_expression = UploadedExpressionFile.objects.create(uploaded_file=uploaded_file,
                                                                        format=UploadedExpressionFile.CUFFDIFF,
                                                                        annotation_level=cuff_diff_file.annotation_level,
                                                                        cuff_diff_file=cuff_diff_file)
        else:
            msg = f"Unknown Expression format '{uploaded_expression.format}'"
            raise ValueError(msg)
        return cuff_diff_file.imported_records

ImportExpressionTask = app.register_task(ImportExpressionTask())  # @UndefinedVariable
