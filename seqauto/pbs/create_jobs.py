from collections import defaultdict
from django.conf import settings
import logging
import os
import stat

from library.file_utils import add_permissions_to_file, mk_path_for_file
from seqauto.job_scripts import get_job_data, create_bash_script
from seqauto.models import VCFFile, SampleSheet, BamFile, \
    SequencingFileType, QC, SampleSheetCombinedVCFFile, IlluminaFlowcellQC, \
    FastQC, Flagstats, JobScript
from snpdb.models import DataState
from seqauto.pbs.pbs_scripts import get_dependency_flags, create_pbs_script


def create_jobs_and_launch_script(seqauto_run, launch_file_types):
    data_states = [DataState.NON_EXISTENT]

    sample_sheets = SampleSheet.objects.filter(sequencingsample__fastq__data_state__in=data_states).distinct()
    illumina_flowcell_qcs = IlluminaFlowcellQC.objects.filter(data_state__in=data_states)
    fastqcs = FastQC.objects.filter(data_state__in=data_states)
    bams = BamFile.objects.filter(data_state__in=data_states)
    flagstats = Flagstats.objects.filter(data_state__in=data_states)
    vcfs = VCFFile.objects.filter(data_state__in=data_states)
    combined_vcfs = SampleSheetCombinedVCFFile.objects.filter(data_state__in=data_states)
    qc = QC.objects.filter(data_state__in=data_states)

    data = {SequencingFileType.SAMPLE_SHEET: sample_sheets,
            SequencingFileType.ILLUMINA_FLOWCELL_QC: illumina_flowcell_qcs,
            SequencingFileType.FASTQC: fastqcs,
            SequencingFileType.BAM: bams,
            SequencingFileType.FLAGSTATS: flagstats,
            SequencingFileType.VCF: vcfs,
            SequencingFileType.COMBINED_VCF: combined_vcfs,
            SequencingFileType.QC: qc,
            SequencingFileType.DATA_MIGRATION: set([b.unaligned_reads.sequencing_sample.sample_sheet for b in bams])}

    job_data_by_file_type = defaultdict(dict)
    for file_type, qs in data.items():
        if file_type in launch_file_types:
            s = get_job_data(seqauto_run, file_type, qs)
            job_data_by_file_type[file_type].update(s)

    launch_script_filename = get_launch_script_with_dependencies(seqauto_run, job_data_by_file_type, launch_file_types,
                                                                 allow_dependencies=True)
    return launch_script_filename


class ScriptWriter:

    def __init__(self, seqauto_run):
        self.seqauto_run = seqauto_run
        self.jobs = 0
        self.lines = []

    def add_line(self, l):
        self.lines.append(l)

    def create_job_script_from_data(self, job_data):
        file_type = job_data["file_type"]
        logging.info(f"Creating job {file_type} - {job_data['path']}")
        copy_fields = ["seqauto_run", "path", "file_type", "out_file"]
        kwargs = {f: job_data[f] for f in copy_fields}
        kwargs["record"] = job_data["record"]
        job_script = JobScript.objects.create(**kwargs)

        name = job_data["name"]
        script_filename = job_data["path"]
        out_file = job_data["out_file"]
        command = job_data["command"]
        cores = job_data["cores"]
        mem = job_data["mem"]

        mk_path_for_file(script_filename)
        mk_path_for_file(out_file)

        if settings.SEQAUTO_USE_PBS:
            create_pbs_script(name, script_filename, out_file, command, cores, mem, job_script.pk)
        else:
            create_bash_script(name, script_filename, out_file, command, cores, mem, job_script.pk)

        return job_script

    def add_job(self, job_data, dependency_scripts=None, allow_dependencies=True):
        if dependency_scripts is None:
            dependency_scripts = []

        self.jobs += 1
        job_script = self.create_job_script_from_data(job_data)

        if allow_dependencies or not dependency_scripts:
            record = job_script.get_record()
            job_submitted_command = f"{settings.SEQAUTO_JOB_SUBMITTED} --job_script_id={job_script.pk}"

            lines = [f"# {record.path}"]

            if settings.SEQAUTO_USE_PBS:
                variable = job_script.get_variable_name()
                dependency_flags = get_dependency_flags(*tuple(dependency_scripts))

                submit_job_cmd = f"{variable}=$(qsub {dependency_flags}{job_script.path})"
                lines.append(submit_job_cmd)

                job_submitted_command += " --job_script_id=${%s}" % variable
                lines.append(job_submitted_command)  # After job is submitted (with the ID)
            else:
                lines.append(job_submitted_command)  # Before job is submitted

                run_job_cmd = f"bash {job_script.path}"
                lines.append(run_job_cmd)

            lines.append('')  # Blank spacer

            self.lines.extend(lines)
        else:
            logging.debug("Skipped job '%s' as it has dependencies", job_script)

        return job_script

    def get_launch_script(self):
        launch_script = None
        if self.jobs:
            job_scripts_dir = self.seqauto_run.get_job_scripts_dir()
            launch_script = os.path.join(job_scripts_dir, "launch_script.sh")
            logging.info(f"Writing {launch_script} for {self.jobs} jobs")
            with open(launch_script, "w") as f:
                for line in self.lines:
                    f.write(line + "\n")

            add_permissions_to_file(launch_script, stat.S_IXUSR)  # Make executable

        return launch_script


def get_launch_script_with_dependencies(seqauto_run, job_data_by_file_type, launch_file_types, allow_dependencies):
    sample_sheets_dict = job_data_by_file_type[SequencingFileType.SAMPLE_SHEET]
    illumina_qcs_dict = job_data_by_file_type[SequencingFileType.ILLUMINA_FLOWCELL_QC]
    fastqcs_dict = job_data_by_file_type[SequencingFileType.FASTQC]
    bams_dict = job_data_by_file_type[SequencingFileType.BAM]
    flagstats_dict = job_data_by_file_type[SequencingFileType.FLAGSTATS]
    vcfs_dict = job_data_by_file_type[SequencingFileType.VCF]
    combined_vcfs_dict = job_data_by_file_type[SequencingFileType.COMBINED_VCF]
    data_migrations_dict = job_data_by_file_type[SequencingFileType.DATA_MIGRATION]
    qc_dict = job_data_by_file_type[SequencingFileType.QC]

    run_sample_sheets = SequencingFileType.SAMPLE_SHEET in launch_file_types
    run_illumina_qc = SequencingFileType.ILLUMINA_FLOWCELL_QC in launch_file_types
    run_fastqc = SequencingFileType.FASTQC in launch_file_types
    run_bams = SequencingFileType.BAM in launch_file_types
    run_flagstats = SequencingFileType.FLAGSTATS in launch_file_types
    run_vcf = SequencingFileType.VCF in launch_file_types
    run_combo_vcf = SequencingFileType.COMBINED_VCF in launch_file_types
    run_data_migration = SequencingFileType.DATA_MIGRATION in launch_file_types
    run_qc = SequencingFileType.QC in launch_file_types

    sw = ScriptWriter(seqauto_run)

    sw.add_line("#!/bin/bash\n\n")
    sw.add_line("set -e\n\n")
    if settings.SEQAUTO_LAUNCH_SCRIPT_PATHS_LIST:
        new_path = ':'.join(settings.SEQAUTO_LAUNCH_SCRIPT_PATHS_LIST)
        sw.add_line("export PATH=${PATH}:%s\n\n" % new_path)

    if settings.SEQAUTO_JOB_COMPLETE:
        cmd = " ".join(settings.MANAGE_COMMAND)
        sw.add_line(f"{cmd} > /dev/null\n\n")  # This will crash out whole script if further vg scripts won't run

    if sample_sheets_dict and run_sample_sheets:
        sw.add_line('echo "Basecalling"\n')
        for job_data in sample_sheets_dict.values():
            sw.add_job(job_data)
        sw.add_line('\n')

    if illumina_qcs_dict and run_illumina_qc:
        sw.add_line('echo "IlluminaFlowcellQCs"\n')
        for job_data in illumina_qcs_dict.values():
            sw.add_job(job_data)
        sw.add_line('\n')

    if fastqcs_dict and run_fastqc:
        sw.add_line('echo "FastQC"\n')
        for job_data in fastqcs_dict.values():
            fastqc = job_data["record"]
            sample_sheet_pk = fastqc.fastq.sequencing_sample.sample_sheet.pk
            dependency_script = sample_sheets_dict.get(sample_sheet_pk)
            fastqc_allow_deps = allow_dependencies and run_sample_sheets
            sw.add_job(job_data, [dependency_script], fastqc_allow_deps)
        sw.add_line('\n')

    bams_scripts_for_sample_sheet = defaultdict(list)
    if bams_dict and run_bams:
        sw.add_line('echo "BAMs"\n')
        for job_data in bams_dict.values():
            bam_file = job_data["record"]
            sample_sheet_pk = bam_file.sample_sheet.pk
            dependency_script = sample_sheets_dict.get(sample_sheet_pk)
            bams_allow_deps = allow_dependencies and run_sample_sheets
            job_script = sw.add_job(job_data, [dependency_script], bams_allow_deps)
            bams_scripts_for_sample_sheet[sample_sheet_pk].append(job_script)
        sw.add_line('\n')

    if flagstats_dict and run_flagstats:
        sw.add_line('echo "Flagstats"\n')
        for job_data in flagstats_dict.values():
            flagstats = job_data["record"]
            bam_pk = flagstats.bam_file.pk
            dependency_script = bams_dict.get(bam_pk)
            flagstats_allow_deps = allow_dependencies and run_bams
            sw.add_job(job_data, [dependency_script], flagstats_allow_deps)
        sw.add_line('\n')

    # Migrations require multiple files to be produced, so store everything under sample sheet PK
    migration_dependencies = {}
    if vcfs_dict and run_vcf:
        sw.add_line('echo "VCFs"\n')
        for job_data in vcfs_dict.values():
            vcf_file = job_data["record"]
            bam_pk = vcf_file.bam_file.pk
            dependency_script = bams_dict.get(bam_pk)
            vcf_allow_deps = allow_dependencies and run_bams
            sw.add_job(job_data, [dependency_script], vcf_allow_deps)
        sw.add_line('\n')

    combined_vcf_script_for_sample_sheet = {}
    if combined_vcfs_dict and run_combo_vcf:
        sw.add_line('echo "Combined VCFs"\n')
        for job_data in combined_vcfs_dict.values():
            combined_vcf_file = job_data["record"]
            sample_sheet = combined_vcf_file.sample_sheet
            dependency_scripts = bams_scripts_for_sample_sheet.get(sample_sheet.pk, [])
            combined_vcf_allow_deps = allow_dependencies and run_bams and run_vcf
            job_script = sw.add_job(job_data, dependency_scripts, combined_vcf_allow_deps)
            combined_vcf_script_for_sample_sheet[sample_sheet.pk] = job_script
        sw.add_line('\n')

    if data_migrations_dict and run_data_migration:
        sw.add_line('echo "Migrations"\n')
        for job_data in data_migrations_dict.values():
            sample_sheet = job_data["record"]
            dependency_script = combined_vcf_script_for_sample_sheet.get(sample_sheet.pk)
            data_migrations_allow_deps = allow_dependencies and run_bams and run_vcf and run_combo_vcf
            job_script = sw.add_job(job_data, [dependency_script], data_migrations_allow_deps)
            migration_dependencies[sample_sheet.pk] = job_script
        sw.add_line('\n')

    if qc_dict and run_qc:
        sw.add_line('echo "QC"\n')
        for job_data in qc_dict.values():
            qc = job_data["record"]
            sample_sheet = qc.bam_file.sample_sheet
            dependency_script = migration_dependencies.get(sample_sheet.pk)
            qc_allow_deps = allow_dependencies and run_bams and run_vcf and run_combo_vcf and run_data_migration
            sw.add_job(job_data, [dependency_script], allow_dependencies=qc_allow_deps)
        sw.add_line('\n')

    return sw.get_launch_script()
