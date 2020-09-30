from celery.app.task import Task
import logging
import subprocess
import time

from library.log_utils import log_traceback
from upload.models import UploadPipeline


class ImportTask(Task):
    """ Subclass this to be able to perform imports """
    abstract = True

    def process_items(self, uploaded_file):
        raise NotImplementedError("Need to override ImportTask.process_items()")

    def run(self, upload_pipeline_id):
        upload_pipeline = UploadPipeline.objects.get(pk=upload_pipeline_id)

        upload_pipeline.start()

        error_message = None
        try:
            start = time.time()

            items_processed = self.process_items(upload_pipeline.uploaded_file)
            if items_processed is None:
                msg = "%s.process_items() returned None!" % str(self.__class__)
                raise ValueError(msg)

            logging.info("Import task processed %d items", items_processed)

            end = time.time()
            processing_seconds_wall_time = end - start
            processing_seconds_cpu_time = processing_seconds_wall_time
            upload_pipeline.success(items_processed,
                                    processing_seconds_wall_time=processing_seconds_wall_time,
                                    processing_seconds_cpu_time=processing_seconds_cpu_time)
        except subprocess.CalledProcessError as e:
            error_message = e.output
        except Exception as e:
            error_message = str(e)

        if error_message:
            logging.error(error_message)
            log_traceback()
            upload_pipeline.error(error_message)
