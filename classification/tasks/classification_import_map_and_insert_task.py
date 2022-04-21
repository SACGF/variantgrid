import subprocess
from pathlib import Path
from typing import List, Dict
import ijson
import json
from celery import Task
from django.conf import settings
from classification.enums import SubmissionSource
from classification.models import UploadedClassificationsUnmapped, UploadedClassificationsUnmappedStatus
from classification.models.classification_import_run import ClassificationImportRun
from library.log_utils import report_message, NotificationBuilder, report_exc_info
from library.utils import batch_iterator, pretty_label
import pathlib
from variantgrid.celery import app


class ClassificationImportMapInsertTask(Task):
    """
    Sends a file to the omni_importer and just asks it to map the file to variantgrid's native JSON classification format.
    Grabs the return file, then imports it.
    (Note that omni_importer is capable of uploading to variantgrid, but more efficient to call omni_importer on cmd
    and get response, then trigger omni_impoter to import and then have to worry about authentication etc)
    """

    @property
    def data_dir(self) -> Path:
        if data_dir_str := settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_DATA_DIR:
            data_dir = pathlib.Path(data_dir_str)
            data_dir.mkdir(parents=True, exist_ok=True)
            if not data_dir.exists():
                raise ValueError(f"{data_dir_str} does not exist")
            return data_dir

    @property
    def omni_importer_dir(self) -> Path:
        if omni_importer_dir_str := settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_APP_DIR:
            omni_importer_dir = pathlib.Path(omni_importer_dir_str)
            if not omni_importer_dir.exists():
                raise ValueError(f"{omni_importer_dir_str} does not exist")
            if not (omni_importer_dir / "main.py").exists():
                raise ValueError(f"main.py doesn't exist inside OMNI IMPORTER DIR {settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_APP_DIR}")
            return omni_importer_dir
        else:
            raise ValueError("Require settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_APP_DIR to be set to do automated imports")

    @staticmethod
    def update_status(upload_file: UploadedClassificationsUnmapped, status: UploadedClassificationsUnmappedStatus):
        upload_file.status = status
        upload_file.save()

    @staticmethod
    def _rm_dir(path: Path):
        if path.is_dir():
            for sub_file in path.iterdir():
                ClassificationImportMapInsertTask._rm_dir(sub_file)
            path.rmdir()
        else:
            path.unlink()

    @staticmethod
    def cleanup_dir(path: Path, delete_dir: bool = False):
        for sub_path in path.iterdir():
            if sub_path.is_dir() and sub_path.name not in {"source", "output"}:
                raise ValueError(f"Working directory has unexpected sub-directory {sub_path.absolute()}")

        if delete_dir:
            ClassificationImportMapInsertTask._rm_dir(path)
        else:
            for sub_path in path.iterdir():
                ClassificationImportMapInsertTask._rm_dir(sub_path)

    def run(self, upload_classifications_unmapped_id: int, import_records: bool = True):
        """
        :param upload_classifications_unmapped_id: ID of the unmapped file
        :param import_records: If False, we run and record the validation, but don't import any of the classifications
        """
        upload_file = UploadedClassificationsUnmapped.objects.get(pk=upload_classifications_unmapped_id)

        from classification.models.classification_inserter import BulkClassificationInserter
        try:
            user = upload_file.user
            ClassificationImportMapInsertTask.update_status(upload_file, UploadedClassificationsUnmappedStatus.Mapping)

            file_data = upload_file.file_data
            sub_folder = f"upload_{upload_file.pk}"
            working_sub_folder = self.data_dir / sub_folder
            working_sub_folder.mkdir(exist_ok=True)
            # if working folder already exists, delete files inside it
            ClassificationImportMapInsertTask.cleanup_dir(working_sub_folder)

            source_dir = working_sub_folder / "source"
            source_dir.mkdir()

            upload_file.file_data.download_to_dir(source_dir, extract_zip=True)
            # next step is to trigger the omni importer to map the file
            # read the mapped file back
            # and import it
            # --file data/path_west/sharmvl_grch38molecular_variants20220326_060243.json --publish logged_in_users --org path_west --lab unit_1 --env prod

            publish = settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_PUBLISH_LEVEL
            include_source = settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_INCLUDE_SOURCE

            args: List[str] = [
                settings.PYTHON_COMMAND, "main.py",
                "--dir", working_sub_folder.absolute(),
                "--publish", publish,
                "--org", upload_file.lab.organization.group_name,
                "--lab", upload_file.lab.group_name.split("/")[1],
                "--env", f"file"
            ]
            if include_source:
                args.append("--include_source")

            # could make it returned the mapped file to stdout
            # but it's handy having the mapped file present
            process = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                cwd=self.omni_importer_dir
            )

            stdout, _ = process.communicate()
            _ = stdout.decode()
            if error_code := process.returncode:
                raise ValueError(f"cwd {self.omni_importer_dir : {' '.join(args)}} returned {error_code}")

            output_dir = working_sub_folder / "output"
            classifications_file = output_dir / "classifications.json"
            validation_summary_file = output_dir / "validation_summary.json"
            validation_list_file = output_dir / "validation_list.json"
            fatal_error = None

            with open(validation_summary_file, 'r') as validation_handle:
                validation_json = json.load(validation_handle)
                fatal_error = validation_json.get('fatal_error')
                upload_file.validation_summary = validation_json

                nb = NotificationBuilder(message="Import Mapped")
                nb.add_header(":currency_exchange: Import Mapped")
                nb.add_field("File ID", f"*{upload_file.pk}* {upload_file.filename}")
                nb.add_field("Lab", str(upload_file.lab))
                nb.add_divider()
                for key, value in validation_json.items():
                    formatted_key = pretty_label(str(key))
                    if key == "fatal_error":
                        formatted_key = ":bangbang: Fatal Error"
                    if key == "message_counts":
                        key = "messages"
                        if not value:
                            value = "- no messages -"
                        else:
                            value_lst = list()
                            for message_type, message_count in value.items():
                                value_lst.append(f"{message_type}: {message_count}")
                            value = "\n".join(value_lst)
                    nb.add_field(formatted_key, value)
                nb.send()

            with open(validation_list_file, 'r') as validation_handle:
                validation_list = json.load(validation_handle)
                upload_file.validation_list = validation_list

            upload_file.save()
            if fatal_error:
                report_message(f"Could not map file for UploadedClassificationsUnmapped ({upload_file.pk})", extra_data={"target": fatal_error})
                ClassificationImportMapInsertTask.update_status(upload_file, UploadedClassificationsUnmappedStatus.Error)
                return

            if import_records:
                import_id = upload_file.file_data.filename + "_" + upload_file.lab.group_name
                with open(classifications_file, 'r') as file_handle:
                    ClassificationImportMapInsertTask.update_status(upload_file, UploadedClassificationsUnmappedStatus.Importing)

                    def row_generator() -> Dict:
                        nonlocal file_handle
                        for record in ijson.items(file_handle, 'records.item'):
                            yield record

                    for batch in batch_iterator(row_generator(), batch_size=50):
                        # record the import
                        import_run = ClassificationImportRun.record_classification_import(
                            identifier=import_id,
                        )
                        import_run.from_file = upload_file

                        bci = BulkClassificationInserter(user=user)
                        for row in batch:
                            response = bci.insert(data=row, submission_source=SubmissionSource.API)
                            import_run.increment_status(response.status)
                        import_run.save()
                        bci.finish()

                    ClassificationImportRun.record_classification_import(
                        identifier=import_id,
                        is_complete=True
                    )

                    # tidy up if everything is working
                    # if things failed (so we never reached this code), we'll still have this directory to investigate
                    ClassificationImportMapInsertTask.cleanup_dir(working_sub_folder, delete_dir=True)

                    ClassificationImportMapInsertTask.update_status(upload_file, UploadedClassificationsUnmappedStatus.Processed)
            else:  # if not import_records
                ClassificationImportMapInsertTask.cleanup_dir(working_sub_folder, delete_dir=True)
                ClassificationImportMapInsertTask.update_status(upload_file, UploadedClassificationsUnmappedStatus.Validated)
        except:
            ClassificationImportMapInsertTask.update_status(upload_file, UploadedClassificationsUnmappedStatus.Error)
            raise


ClassificationImportMapInsertTask = app.register_task(ClassificationImportMapInsertTask())
