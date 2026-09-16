import logging
import os

from django.core.management.base import BaseCommand, CommandError

from snpdb.models import ProcessingStatus
from upload.models import UploadedFileTypes, UploadPipeline
from upload.uploaded_file_type import retry_upload_pipeline


class Command(BaseCommand):
    """ Re-runs upload pipelines. The retry button on the pipeline page is the usual way, but it is
        hidden on some deployments (Shariant) and for internally generated files (classification
        import, liftover) there's no page a user would go looking at anyway """
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--uploaded_file_type', help="One of the UploadedFileTypes codes")
        parser.add_argument('--upload_pipeline_id', type=int, action='append',
                            help="Pipeline to reload, may be given multiple times")
        parser.add_argument('--status', help="Only reload pipelines in this ProcessingStatus, eg 'E' for Error")
        parser.add_argument('--dry-run', action='store_true', help="List what would be reloaded")

    def handle(self, *args, **options):
        uploaded_file_type = options['uploaded_file_type']
        upload_pipeline_ids = options['upload_pipeline_id']
        status = options['status']
        dry_run = options['dry_run']

        if not (uploaded_file_type or upload_pipeline_ids):
            raise CommandError("Provide --uploaded_file_type and/or --upload_pipeline_id")

        kwargs = {}
        if uploaded_file_type:
            uft_dict = dict(UploadedFileTypes.choices)
            uft_description = uft_dict.get(uploaded_file_type)
            if uft_description is None:
                script_name = os.path.basename(__file__)
                valid_ufts = ','.join(sorted(uft_dict))
                msg = f"Usage: {script_name} --uploaded_file_type=X (where X is one of {valid_ufts})"
                raise CommandError(msg)
            logging.info("Reloading UFPPs of type '%s'", uft_description)
            kwargs["file_upload__file_type"] = uploaded_file_type

        if upload_pipeline_ids:
            kwargs["pk__in"] = upload_pipeline_ids

        if status:
            if status not in dict(ProcessingStatus.choices):
                valid_statuses = ','.join(sorted(dict(ProcessingStatus.choices)))
                raise CommandError(f"--status must be one of {valid_statuses}")
            kwargs["status"] = status

        for upload_pipeline in UploadPipeline.objects.filter(**kwargs).order_by("pk"):
            filename = upload_pipeline.file_upload.get_filename()
            if filename and not os.path.exists(filename):
                # Nothing to feed the pipeline - for a classification import, admin variant re-match
                # is the way back (@see classification.classification_import.reattempt_variant_matching)
                logging.warning("Skipping pipeline %d - file '%s' no longer on disk",
                                upload_pipeline.pk, filename)
                continue

            if dry_run:
                logging.info("Would reload pipeline %d (%s, %s)", upload_pipeline.pk,
                             upload_pipeline.get_file_type_display(), upload_pipeline.get_status_display())
                continue

            retry_upload_pipeline(upload_pipeline)
