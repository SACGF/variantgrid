import logging
import subprocess
import time

from celery.app.task import Task

from library.log_utils import get_traceback
from library.utils import format_called_process_error
from upload.models.models import UploadPipeline, UploadUserError


class ImportRequiresUserInputException(Exception):
    """ Raise from ImportTask.process_items when the file can't finish importing until the user sets
        something it doesn't declare (eg a BED file's genome build) - the message says what """


class ImportTask(Task):
    """ Subclass this to be able to perform imports """
    abstract = True

    def process_items(self, file_upload):
        raise NotImplementedError("Need to override ImportTask.process_items()")

    def run(self, upload_pipeline_id):
        upload_pipeline = UploadPipeline.objects.get(pk=upload_pipeline_id)

        upload_pipeline.start()

        error_message = None
        try:
            start = time.time()

            items_processed = self.process_items(upload_pipeline.file_upload)
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
        except ImportRequiresUserInputException as e:
            upload_pipeline.terminate_early_for_user_input(str(e))
        except UploadUserError as e:
            upload_pipeline.user_error(e)
        except subprocess.CalledProcessError as e:
            error_message = format_called_process_error(e)
        except:
            error_message = get_traceback()

        if error_message:
            logging.error(error_message)
            upload_pipeline.error(error_message)
