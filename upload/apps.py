from django.apps import AppConfig

class UploadConfig(AppConfig):
    name = 'upload'

    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from annotation.signals.manual_signals import annotation_run_complete_signal
        from upload.signals.signal_handlers import annotation_run_complete_signal_handler
        # pylint: enable=import-outside-toplevel,unused-import

        annotation_run_complete_signal.connect(annotation_run_complete_signal_handler)
