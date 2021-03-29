from django.conf import settings
from django.template.loader import render_to_string
import collections
import logging
import os

from seqauto.models import SequencingFileType

BASH_TEMPLATE = os.path.join(settings.SEQAUTO_DIR, "templates", "bash_script.template")

COMMAND_PATTERNS = {SequencingFileType.SAMPLE_SHEET: settings.SEQAUTO_SAMPLE_SHEET_COMMAND,
                    SequencingFileType.ILLUMINA_FLOWCELL_QC: settings.SEQAUTO_ILLUMINA_FLOWCELL_QC_COMMAND,
                    SequencingFileType.FASTQC: settings.SEQAUTO_FASTQC_COMMAND,
                    SequencingFileType.BAM: settings.SEQAUTO_BAM_COMMAND,
                    SequencingFileType.FLAGSTATS: settings.SEQAUTO_FLAGSTATS_COMMAND,
                    SequencingFileType.VCF: settings.SEQAUTO_VCF_COMMAND,
                    SequencingFileType.COMBINED_VCF: settings.SEQAUTO_COMBINED_VCF_COMMAND,
                    SequencingFileType.QC: settings.SEQAUTO_QC_COMMAND,
                    SequencingFileType.DATA_MIGRATION: settings.SEQAUTO_MIGRATE_COMMAND}

DEFAULT_CORES = 1
CORES = {SequencingFileType.SAMPLE_SHEET: 16,
         SequencingFileType.BAM: 4,
         SequencingFileType.VCF: 4,
         SequencingFileType.QC: 4}
DEFAULT_MEM = 16
MEM = {SequencingFileType.ILLUMINA_FLOWCELL_QC: 4,
       SequencingFileType.FASTQC: 4,
       SequencingFileType.FLAGSTATS: 4}


def create_bash_script(name, filename, out_file, command, cores, mem, job_script_pk):
    pythonpaths = [settings.BASE_DIR]
    pbs_pythonpath = settings.SEQAUTO_SCRIPT_PARAMS.get("pythonpath")
    if pbs_pythonpath:
        pythonpaths.append(pbs_pythonpath)
    pythonpath = ':'.join(pythonpaths)

    paths = [os.path.join(settings.BASE_DIR, "scripts"),
             os.path.join(settings.SEQAUTO_DIR, "scripts")]
    path = ':'.join(paths)

    job_complete = settings.SEQAUTO_JOB_COMPLETE
    if job_complete and settings.SEQAUTO_VIRTUALENV_RUNNER:
        job_complete = f"{settings.SEQAUTO_VIRTUALENV_RUNNER} {job_complete}"

    context = {"name": name,
               "out": out_file,
               "mem": mem,
               "cores": cores,
               "path": path,
               "pythonpath": pythonpath,
               "scratch_base_dir": settings.SEQAUTO_SCRATCH_BASE_DIR,
               "command": command,
               "job_script_complete": job_complete,
               "job_script_pk": job_script_pk}
    script_contents = render_to_string(BASH_TEMPLATE, context)

    with open(filename, 'w') as f:
        f.write(script_contents)


def get_job_data(seqauto_run, file_type, qs):
    logging.info("Create PBS Scripts for %s", SequencingFileType(file_type).label)
    job_data = {}
    pattern = COMMAND_PATTERNS.get(file_type)
    cores = CORES.get(file_type, DEFAULT_CORES)
    mem = MEM.get(file_type, DEFAULT_MEM)

    scripts_dir = seqauto_run.get_job_scripts_dir()

    if pattern is not None:
        for record in qs:
            params = settings.SEQAUTO_SCRIPT_PARAMS.copy()
            params["base_dir"] = settings.BASE_DIR
            params.update(record.get_params())

            if isinstance(pattern, str):
                command = pattern % params
            elif isinstance(pattern, collections.Callable):
                command = pattern(params)
            else:
                msg = f"Unknown pattern for {file_type}: {type(pattern)}"
                raise ValueError(msg)

            name = f"{file_type}_{record.pk}"
            if settings.SEQAUTO_USE_PBS:
                extension = ".csh"
            else:
                extension = ".sh"

            script_filename = os.path.join(scripts_dir, name + extension)
            out_file = os.path.join(settings.SEQAUTO_JOB_SCRIPTS_OUT_DIR, name)
            jd = {"seqauto_run": seqauto_run,
                  "path": script_filename,
                  "file_type": file_type,
                  "out_file": out_file,
                  "record": record,
                  "name": name,
                  "command": command,
                  "cores": cores,
                  "mem": mem}
            job_data[record.pk] = jd
    else:
        logging.info("Skipped scripts for %s", SequencingFileType(file_type).label)

    logging.info("Generated scripts for %d records", len(job_data))
    if job_data:
        data = next(iter(job_data.values()))
        print(f"A job was written to {data['path']}")

    return job_data
