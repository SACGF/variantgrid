import logging
import os
import subprocess

import celery
from django.conf import settings
from django.utils import timezone

from library.utils.file_utils import name_from_filename, mk_path
from library.log_utils import get_traceback
from seqauto.models import SeqAutoRun
from seqauto.models.models_enums import SequencingFileType
from seqauto.pbs.create_jobs import create_jobs_and_launch_script
from seqauto.sequencing_files.create_resource_models import create_resource_models

SEQAUTO_SCRIPTS = [
    (SequencingFileType.SAMPLE_SHEET, settings.SEQAUTO_FLOWCELL_SCRIPT),
    (SequencingFileType.ILLUMINA_FLOWCELL_QC, settings.SEQAUTO_ILLUMINATE_QC),
    (SequencingFileType.FASTQ, settings.SEQAUTO_FASTQ_SCRIPT),
    (SequencingFileType.FASTQC, settings.SEQAUTO_FASTQC_SCRIPT),
    (SequencingFileType.BAM, settings.SEQAUTO_BAM_SCRIPT),
    (SequencingFileType.FLAGSTATS, settings.SEQAUTO_FLAGSTATS_SCRIPT),
    (SequencingFileType.VCF, settings.SEQAUTO_VCF_SCRIPT),
    (SequencingFileType.QC, settings.SEQAUTO_QC_SCRIPT),
]


def get_seqauto_scripts(only_process_file_types):
    seqauto_scripts = []
    for file_type, script_name in SEQAUTO_SCRIPTS:
        if only_process_file_types and file_type not in only_process_file_types:
            continue
        if script_name:
            seqauto_scripts.append((file_type, script_name))

    return seqauto_scripts


def scan_resources(seqauto_run, seqauto_scripts):
    for _, script_name in seqauto_scripts:
        script = os.path.join(settings.SEQAUTO_SCRIPTS_DIR, script_name)
        output_filename = os.path.join(seqauto_run.scan_resources_dir, "%s.txt" % name_from_filename(script_name))
        with open(output_filename, "w") as out_f:
            subprocess.check_call([script, seqauto_run.scan_resources_dir], stdout=out_f)

    #shutil.rmtree(scan_resources_dir)


@celery.shared_task(track_started=True)
def scan_run_jobs(only_process_file_types=None, only_launch_file_types=None, run_launch_script=None, fake_data=None):
    if run_launch_script is None:
        run_launch_script = settings.SEQAUTO_SCAN_RUN_SCRIPTS
    if only_launch_file_types is None:
        only_launch_file_types = [SequencingFileType.ILLUMINA_FLOWCELL_QC]

    params = (only_process_file_types, only_launch_file_types, run_launch_script)
    logging.info("only_process_file_types=%s, only_launch_file_types=%s, run_launch_script=%s", *params)

    task_id = scan_run_jobs.request.id
    seqauto_run = SeqAutoRun.objects.create(task_id=task_id, fake_data=fake_data)
    exception = None
    try:
        # Scan resources
        scan_resources_dir = seqauto_run.get_scan_resources_dir()
        mk_path(scan_resources_dir)

        seqauto_run.scan_start = timezone.now()
        seqauto_run.scan_resources_dir = scan_resources_dir
        seqauto_run.save()

        process_seqauto_scripts = get_seqauto_scripts(only_process_file_types)
        scan_resources(seqauto_run, process_seqauto_scripts)

        # Create Models
        seqauto_run.create_models_start = timezone.now()
        seqauto_run.save()

        create_resource_models(seqauto_run, process_seqauto_scripts)

        # Jobs
        all_sequencing_file_types = set(dict(SEQAUTO_SCRIPTS)) | {SequencingFileType.COMBINED_VCF}
        launch_file_types = only_launch_file_types or all_sequencing_file_types
        seqauto_run.scripts_and_jobs_start = timezone.now()

        job_launch_script_filename = create_jobs_and_launch_script(seqauto_run, launch_file_types)
        if job_launch_script_filename:
            logging.info("Wrote launch script %s", job_launch_script_filename)

            seqauto_run.job_launch_script_filename = job_launch_script_filename
            if run_launch_script:
                try:
                    logging.info("Running launch script %s", job_launch_script_filename)
                    subprocess.check_output([job_launch_script_filename], stderr=subprocess.STDOUT)
                except Exception as e:
                    error_exception = "Error running launch script: %s" % get_traceback()
                    output = getattr(e, "output", None)  # CalledProcessError
                    if output:
                        error_exception += "\n" + output.decode()

                    seqauto_run.error_exception = error_exception
                    exception = e

        seqauto_run.finish_date = timezone.now()
    except Exception as e:
        seqauto_run.error_exception = get_traceback()
        exception = e

    if seqauto_run.error_exception:
        logging.error(seqauto_run.error_exception)
    elif not settings.DEBUG:
        seqauto_run.remove_scan_resources_dir()

    seqauto_run.save()

    if exception:
        raise exception  # for other error handlers (eg Rollbar) to report
