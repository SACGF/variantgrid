import subprocess
from pathlib import Path
from typing import List, Dict
import ijson
from celery import Task
from django.conf import settings
from classification.enums import SubmissionSource
from classification.models import UploadedFileLab, UploadedFileLabStatus
from classification.models.classification_import_run import ClassificationImportRun
from library.utils import batch_iterator
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
    def update_status(upload_file: UploadedFileLab, status: UploadedFileLabStatus):
        upload_file.status = status
        upload_file.save()

    @staticmethod
    def cleanup_dir(path: Path, delete_dir: bool = False):
        if existing_files := list(path.iterdir()):
            if len(existing_files) > 4:
                raise ValueError(f"Working subdirectory {path} has {len(existing_files)} files/folders in it, too dangerous to delete")
            if sub_dirs := [ef for ef in existing_files if ef.is_dir()]:
                raise ValueError(
                    f"Working subdirectory {path} has subdirectories {sub_dirs}, too dangerous to delete")

            for sub_file in existing_files:
                sub_file.unlink()
        if delete_dir:
            path.rmdir()

    def run(self, upload_file_id: int):
        upload_file = UploadedFileLab.objects.get(pk=upload_file_id)

        from classification.models.classification_inserter import BulkClassificationInserter
        try:
            user = upload_file.user
            ClassificationImportMapInsertTask.update_status(upload_file, UploadedFileLabStatus.Mapping)

            file_data = upload_file.file_data
            filename = file_data.filename

            sub_folder = f"upload_{upload_file.pk}"
            working_sub_folder = self.data_dir / sub_folder
            working_sub_folder.mkdir(exist_ok=True)
            # if working folder already exists, delete files inside it
            ClassificationImportMapInsertTask.cleanup_dir(working_sub_folder)

            file_handle = working_sub_folder / filename
            mapped_file = working_sub_folder / f"mapped_classifications.json"

            upload_file.file_data.download_to(file_handle)
            # next step is to trigger the omni importer to map the file
            # read the mapped file back
            # and import it
            # --file data/path_west/sharmvl_grch38molecular_variants20220326_060243.json --publish logged_in_users --org path_west --lab unit_1 --env prod
            args: List[str] = [
                settings.PYTHON_COMMAND, "main.py",
                "--file", file_handle,
                "--publish", "logged_in_users",
                "--org", upload_file.lab.organization.group_name,
                "--lab", upload_file.lab.group_name.split("/")[1],
                "--env", f"file:{mapped_file}"
            ]
            # could make it returned the mapped file to stdout
            # but it's handy having the mapped file present
            process = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                cwd=self.omni_importer_dir
            )

            stdout, _ = process.communicate()
            stdout_str = stdout.decode()
            if error_code := process.returncode:
                raise ValueError(f"cwd {self.omni_importer_dir : {' '.join(args)}} returned {error_code}")

            if mapped_file.exists():
                import_id = upload_file.file_data.filename + "_" + upload_file.lab.group_name
                with open(mapped_file, 'r') as file_handle:
                    ClassificationImportMapInsertTask.update_status(upload_file, UploadedFileLabStatus.Importing)

                    def row_generator() -> Dict:
                        nonlocal file_handle
                        for record in ijson.items(file_handle, 'records.item'):
                            yield record

                    for batch in batch_iterator(row_generator(), batch_size=50):
                        # record the import
                        ClassificationImportRun.record_classification_import(
                            identifier=import_id,
                            add_row_count=len(batch))

                        bci = BulkClassificationInserter(user=user)
                        for row in batch:
                            bci.insert(data=row, submission_source=SubmissionSource.API)
                        bci.finish()

                    ClassificationImportRun.record_classification_import(
                        identifier=import_id,
                        add_row_count=0,
                        is_complete=True
                    )

                    # tidy up if everything is working
                    # if things failed (so we never reached this code), we'll still have this directory to investigate
                    ClassificationImportMapInsertTask.cleanup_dir(working_sub_folder, delete_dir=True)

                    ClassificationImportMapInsertTask.update_status(upload_file, UploadedFileLabStatus.Processed)
            else:
                ClassificationImportMapInsertTask.update_status(upload_file, UploadedFileLabStatus.Error)
        except:
            ClassificationImportMapInsertTask.update_status(upload_file, UploadedFileLabStatus.Error)
            raise


ClassificationImportMapInsertTask = app.register_task(ClassificationImportMapInsertTask())