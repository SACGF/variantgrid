from django.apps import AppConfig

class UploadConfig(AppConfig):
    name = 'upload'

    def ready(self):
        from annotation.signals import annotation_run_complete_signal
        from upload.signals.signal_handlers import annotation_run_complete_signal_handler

        annotation_run_complete_signal.connect(annotation_run_complete_signal_handler)
