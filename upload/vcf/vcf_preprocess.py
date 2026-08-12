import glob
import logging
import os
from dataclasses import dataclass
from subprocess import PIPE, CalledProcessError, Popen
from typing import Optional

import pandas as pd
from django.conf import settings
from django.utils import timezone

from library.django_utils.django_file_utils import get_import_processing_filename
from library.genomics.vcf_enums import (
    INFO_LIFTOVER_SWAPPED_REF_ALT,
    UNDECLARED_FILTERS_INFO,
    UNDECLARED_FILTERS_SEPARATOR,
)
from library.genomics.vcf_utils import write_cleaned_vcf_header
from library.utils.file_utils import name_from_filename
from upload.models import (
    ModifiedImportedVariant,
    ModifiedImportedVariants,
    SimpleVCFImportInfo,
    ToolVersion,
    UploadStep,
    UploadStepMultiFileOutput,
    UploadStepTaskType,
    VCFPipelineStage,
    VCFSkippedContig,
    VCFSkippedContigs,
)
from upload.tasks.vcf.unknown_variants_task import (
    AnnotateImportedVCFTask,
    SeparateUnknownVariantsTask,
)

MAX_STDERR_OUTPUT = 5000  # How much stderr output per process to store in DB

VCF_CLEAN_AND_FILTER_SUB_STEP = "vcf_clean_and_filter"
SORT_VCF_SUB_STEP = "sort_vcf"
VCF_CLEAN_ALTS_SUB_STEP = "vcf_clean_alts"
REMOVE_HEADER_SUB_STEP = "remove_header"
SPLIT_VCF_SUB_STEP = "split_vcf"

UNSORTED_VCF_MESSAGE = "VCF was not sorted - records have been sorted for import"


def get_bcftools_tool_version(bcftools_command):
    cmd = [bcftools_command, "--version"]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    if p.returncode:
        raise CalledProcessError(p.returncode, cmd, output=stdout, stderr=stderr)

    output_list = stdout.split("\n")
    bcftools_version = output_list[0]
    if not bcftools_version.startswith("bcftools"):
        raise ValueError(f"Expected to find bcftools on 1st line of '{' '.join(cmd)}' output: {bcftools_version}")
    htslib_version = output_list[1] if len(output_list) > 1 else ""
    if "htslib" not in htslib_version:
        raise ValueError(f"Expected to find htslib version on 2nd line of '{' '.join(cmd)}' output: {htslib_version}")

    version = f"{bcftools_version}, {htslib_version}"
    tool_version, _ = ToolVersion.objects.get_or_create(name='bcftools', version=version)
    return tool_version


def create_sub_step(upload_step, sub_step_name, sub_step_commands):
    start_date = timezone.now()
    tool_version = get_bcftools_tool_version(settings.BCFTOOLS_COMMAND)
    command_line = ' '.join(sub_step_commands)
    return UploadStep.objects.create(name=sub_step_name,
                                     script=command_line,
                                     upload_pipeline=upload_step.upload_pipeline,
                                     sort_order=upload_step.sort_order,
                                     parent_upload_step=upload_step,
                                     tool_version=tool_version,
                                     start_date=start_date,
                                     task_type=UploadStepTaskType.TOOL)


def _write_cleaned_header(genome_build, upload_pipeline, vcf_filename) -> str:
    # We clean/replace header. Used for initial read of vcf_clean_and_filter then added on top of every splice vcf
    clean_vcf_dir = upload_pipeline.get_pipeline_processing_subdir("clean_vcf_dir")
    cleaned_vcf_header_filename = os.path.join(clean_vcf_dir, name_from_filename(vcf_filename, remove_gz=True) + ".header.vcf")
    tag = ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG
    HEADERS = [
        f'##INFO=<ID={tag},Number=1,Type=String,Description="Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX">',
        f'##INFO=<ID={INFO_LIFTOVER_SWAPPED_REF_ALT},Number=0,Type=Flag,Description="Whether we have switched REF/ALT due to SWAP">',
        f'##INFO=<ID={UNDECLARED_FILTERS_INFO},Number=1,Type=String,Description="FILTER values undeclared in the source header, moved here by vcf_clean_and_filter. Separator: {UNDECLARED_FILTERS_SEPARATOR}">',
    ]
    write_cleaned_vcf_header(genome_build, vcf_filename, cleaned_vcf_header_filename,
                             new_info_lines=HEADERS, standard_contigs_only=True)
    return cleaned_vcf_header_filename


@dataclass
class PreprocessFiles:
    """ Paths the pipe stages read/write. Grouped so a sort retry can clear everything the failed
        attempt left behind before running again """
    vcf_filename: str
    cleaned_vcf_header: str
    split_vcf_dir: str
    unsorted_marker: str
    skipped_contigs_stats: str
    skipped_records_stats: str
    skipped_filters_stats: str
    skipped_alts_stats: str
    converted_alts_stats: str
    swap_skipped: Optional[str] = None

    @property
    def written_by_pipe(self) -> list[str]:
        filenames = [self.unsorted_marker, self.skipped_contigs_stats, self.skipped_records_stats,
                     self.skipped_filters_stats, self.skipped_alts_stats, self.converted_alts_stats]
        if self.swap_skipped:
            filenames.append(self.swap_skipped)
        return filenames


def _get_preprocess_files(upload_pipeline, vcf_filename: str, cleaned_vcf_header: str,
                          disable_swap: bool) -> PreprocessFiles:
    vcf_name = name_from_filename(vcf_filename, remove_gz=True)

    def _processing_filename(suffix):
        return get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.{suffix}")

    swap_skipped_file = _processing_filename("swap_skipped.tsv") if disable_swap else None
    return PreprocessFiles(vcf_filename=vcf_filename,
                           cleaned_vcf_header=cleaned_vcf_header,
                           split_vcf_dir=upload_pipeline.get_pipeline_processing_subdir("split_vcf"),
                           unsorted_marker=_processing_filename("unsorted.txt"),
                           skipped_contigs_stats=_processing_filename("skipped_contigs.tsv"),
                           skipped_records_stats=_processing_filename("skipped_records.tsv"),
                           skipped_filters_stats=_processing_filename("skipped_filters.tsv"),
                           skipped_alts_stats=_processing_filename("skipped_alts.tsv"),
                           converted_alts_stats=_processing_filename("converted_alts.tsv"),
                           swap_skipped=swap_skipped_file)


def _build_pipe_commands(upload_step, files: PreprocessFiles, disable_swap=False, sort_vcf=False) -> tuple[dict, dict]:
    """ Returns (pipe_commands, sub_steps) - pipe_commands dict order is the pipeline order.
        sort_vcf: splice 'bcftools sort' in after the clean/filter stage, which then skips its sorted check """
    genome_build = upload_step.genome_build
    upload_pipeline = upload_step.upload_pipeline
    vcf_name = name_from_filename(files.vcf_filename, remove_gz=True)

    pipe_commands = {}  # Use dict order to decide pipeline order
    if files.vcf_filename.endswith(".gz"):
        pipe_commands["zcat"] = [settings.BASH_ZCAT, files.vcf_filename]
    else:
        pipe_commands["cat"] = ["cat", files.vcf_filename]

    manage_command = settings.MANAGE_COMMAND
    read_variants_cmd = manage_command + ["vcf_clean_and_filter",
                                          "--genome-build",
                                          genome_build.name,
                                          "--replace-header",
                                          files.cleaned_vcf_header,
                                          "--skipped-contigs-stats-file",
                                          files.skipped_contigs_stats,
                                          "--skipped-records-stats-file",
                                          files.skipped_records_stats,
                                          "--skipped-filters-stats-file",
                                          files.skipped_filters_stats]
    if sort_vcf:
        read_variants_cmd.append("--allow-unsorted")
    else:
        read_variants_cmd += ["--unsorted-marker-file", files.unsorted_marker]
    pipe_commands[VCF_CLEAN_AND_FILTER_SUB_STEP] = read_variants_cmd

    if sort_vcf:
        # Sorts into the order of the cleaned header, whose contigs the stage above has already renamed to.
        # Trailing separator so bcftools treats this as a directory rather than a temp file prefix
        sort_temp_dir = os.path.join(upload_pipeline.get_pipeline_processing_subdir("sort_tmp"), "")
        pipe_commands[SORT_VCF_SUB_STEP] = [settings.BCFTOOLS_COMMAND, "sort", "--temp-dir", sort_temp_dir]

    pipe_commands[UploadStep.NORMALIZE_SUB_STEP] = [
        settings.BCFTOOLS_COMMAND, "norm",
        "--multiallelics=-",
        # We don't remove duplicates due to:
        # * https://github.com/samtools/bcftools/issues/2225 - rmdup removes --old-rec-tag so lose normalize info
        # "--rm-dup=exact",
        "--check-ref=s",  # Set ref (ie replace N with actual ref base)
        f"--old-rec-tag={ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG}",
        f"--fasta-ref={genome_build.reference_fasta}", "-",
    ]
    pipe_commands[REMOVE_HEADER_SUB_STEP] = [settings.BCFTOOLS_COMMAND, "view", "--no-header", "-"]

    # We have performed multi-allelic splitting above, so that we can now filter out bad alts (but leave rest of record)
    clean_alts_cmd = manage_command + [
        "vcf_clean_alts",
        "--skipped-records-stats-file", files.skipped_alts_stats,
        "--converted-records-stats-file", files.converted_alts_stats,
    ]
    if disable_swap:
        clean_alts_cmd += ["--disable-swap", "--swap-skipped-file", files.swap_skipped]
    pipe_commands[VCF_CLEAN_ALTS_SUB_STEP] = clean_alts_cmd

    # Split up the VCF
    split_file_rows = upload_step.split_file_rows or settings.VCF_IMPORT_FILE_SPLIT_ROWS
    pipe_commands[SPLIT_VCF_SUB_STEP] = ["split", "-", vcf_name, "--additional-suffix=.vcf.gz", "--numeric-suffixes",
                                         "--lines", str(split_file_rows),
                                         "--filter='bash -c \"set -eo pipefail; { cat $VG_HEADER_FILE; cat; } | bgzip -c > $VG_SPLIT_VCF_DIR/$FILE\"'"]

    sub_steps = {}
    for sub_step_name in [VCF_CLEAN_AND_FILTER_SUB_STEP, SORT_VCF_SUB_STEP, UploadStep.NORMALIZE_SUB_STEP]:
        if sub_step_commands := pipe_commands.get(sub_step_name):
            sub_steps[sub_step_name] = create_sub_step(upload_step, sub_step_name, sub_step_commands)

    return pipe_commands, sub_steps


def _run_pipe(pipe_commands: dict, sub_steps: dict, split_env: dict, upload_pipeline):
    piped_command = ' | '.join([' '.join(command_with_args) for command_with_args in pipe_commands.values()])

    # if POPEN_SHELL=True we're running this all as one command as native shell text
    # this is not recommended, less portable and if commands have errors
    if settings.VCF_IMPORT_PREPROCESS_POPEN_SHELL:
        logging.info("single_commands: %s" % " | ".join([' '.join(x) for x in pipe_commands.values()]))
        # Use pipefail so that a non-zero exit from any stage in the pipe is reported,
        # not just the last command (see #3813 - bcftools errors were silently swallowed)
        shell_command = f"set -o pipefail; {piped_command}"
        p = Popen(
            shell_command,
            shell=True,
            executable="/bin/bash",  # pipefail requires bash, not sh
            stdout=PIPE,
            stderr=PIPE,
            env=split_env,
        )
        p_stdout, p_stderr = p.communicate()
        logging.info("single command pipe/shell completed - return code: %d", p.returncode)

        p_stdout = p_stdout.decode()
        if p_stdout:
            logging.info(p_stdout)

        p_stderr = p_stderr.decode()
        if p_stderr:
            logging.error(p_stderr)

        if p.returncode:
            raise CalledProcessError(p.returncode, piped_command, stderr=p_stderr)

    else:
        pipes = {}
        stderr_filenames = {}
        p = None
        for sub_step_name, cmd in pipe_commands.items():
            stderr_filename = get_import_processing_filename(upload_pipeline.pk, f"stderr_out_{sub_step_name}.log")
            stderr_f = open(stderr_filename, "wb")
            kwargs = {"stdout": PIPE,
                      "stderr": stderr_f}
            if p:
                kwargs["stdin"] = p.stdout
            if sub_step_name == SPLIT_VCF_SUB_STEP:
                kwargs["env"] = split_env
            p = Popen(cmd, **kwargs)
            pipes[sub_step_name] = p
            stderr_filenames[sub_step_name] = stderr_filename

        logging.info("pipe_commands: %s", " | ".join([' '.join(x) for x in pipe_commands.values()]))

        for sub_step_name, p_cmd in reversed(pipe_commands.items()):
            p = pipes[sub_step_name]
            p.communicate()
            stderr_output = None
            stderr_filename = stderr_filenames.get(sub_step_name)
            if stderr_filename:
                with open(stderr_filename) as sf:
                    stderr_output = sf.read()
                if stderr_output and len(stderr_output) > MAX_STDERR_OUTPUT:
                    half = MAX_STDERR_OUTPUT // 2
                    stderr_output = stderr_output[:half] + "\n ... \n" + stderr_output[-half:]
                sub_step = sub_steps.get(sub_step_name)
                if sub_step:
                    sub_step.output_text = stderr_output
                    sub_step.save()
                else:
                    logging.warning("%s stderr:", p_cmd)
                    logging.warning(stderr_output)

            if p.returncode:
                raise CalledProcessError(p.returncode, p_cmd, stderr=stderr_output)


def _reset_for_retry(files: PreprocessFiles, sub_steps: dict):
    """ Put the pipeline back to how it was before the failed attempt, so the next one starts clean """
    for sub_step in sub_steps.values():
        sub_step.delete()

    for split_filename in glob.glob(os.path.join(files.split_vcf_dir, "*")):
        os.remove(split_filename)

    for filename in files.written_by_pipe:
        if os.path.exists(filename):
            os.remove(filename)


def preprocess_vcf(upload_step, annotate_gnomad_af=False, disable_swap=False):
    genome_build = upload_step.genome_build
    _ = genome_build.reference_fasta  # Fails if not available.

    vcf_filename = upload_step.input_filename
    if not os.path.exists(vcf_filename):
        raise FileNotFoundError(f"Can't access vcf: '{vcf_filename}'")

    upload_pipeline = upload_step.upload_pipeline
    cleaned_vcf_header_filename = _write_cleaned_header(genome_build, upload_pipeline, vcf_filename)
    files = _get_preprocess_files(upload_pipeline, vcf_filename, cleaned_vcf_header_filename, disable_swap)

    # Pass paths via env vars rather than interpolating them into the shell --filter string,
    # so that metacharacters in paths cannot affect shell parsing.
    # split runs the filter via sh -c; $VG_HEADER_FILE/$VG_SPLIT_VCF_DIR are expanded there.
    split_env = {**os.environ,
                 'VG_HEADER_FILE': files.cleaned_vcf_header,
                 'VG_SPLIT_VCF_DIR': files.split_vcf_dir}

    if os.path.exists(files.unsorted_marker):
        os.remove(files.unsorted_marker)  # Left behind by an earlier attempt at this step

    pipe_commands, sub_steps = _build_pipe_commands(upload_step, files, disable_swap=disable_swap)
    sorted_vcf = False
    try:
        _run_pipe(pipe_commands, sub_steps, split_env, upload_pipeline)
    except CalledProcessError:
        if not os.path.exists(files.unsorted_marker):
            raise
        with open(files.unsorted_marker) as f:
            logging.info("Re-running through '%s sort': %s", settings.BCFTOOLS_COMMAND, f.read().strip())
        _reset_for_retry(files, sub_steps)
        pipe_commands, sub_steps = _build_pipe_commands(upload_step, files, disable_swap=disable_swap, sort_vcf=True)
        _run_pipe(pipe_commands, sub_steps, split_env, upload_pipeline)
        sorted_vcf = True

    clean_sub_step = sub_steps[VCF_CLEAN_AND_FILTER_SUB_STEP]
    if sorted_vcf:
        SimpleVCFImportInfo.objects.create(upload_step=clean_sub_step, message_string=UNSORTED_VCF_MESSAGE)

    if os.path.exists(files.skipped_contigs_stats):
        df = pd.read_csv(files.skipped_contigs_stats, header=None, sep='\t', index_col=0)
        if not df.empty:
            import_info = VCFSkippedContigs.objects.create(upload_step=clean_sub_step)

            for contig, count in df.iloc[:, 0].items():
                VCFSkippedContig.objects.create(import_info=import_info,
                                                contig=contig,
                                                num_skipped=count)

    _store_vcf_stats(files.skipped_records_stats, clean_sub_step, "records")
    _store_vcf_stats(files.skipped_filters_stats, clean_sub_step, "FILTER")
    _store_vcf_stats(files.skipped_alts_stats, clean_sub_step, "ALTs")
    _store_vcf_stats(files.converted_alts_stats, clean_sub_step, "ALTs", operation='Converted')

    _schedule_split_file_steps(upload_step, files.split_vcf_dir, annotate_gnomad_af=annotate_gnomad_af)

    return files.swap_skipped


def _schedule_split_file_steps(upload_step, split_vcf_dir: str, annotate_gnomad_af=False):
    """ One SeparateUnknownVariantsTask per split VCF, plus the UploadStepMultiFileOutput rows that
        ScheduleMultiFileOutputTasksTask fans the data-insertion tasks out over. Shared by every
        preprocess variant - a pipeline that skips the bcftools stages still needs all of this. """
    upload_pipeline = upload_step.upload_pipeline
    # Create this here so downstream tasks (running in parallel) can all link against the same one
    ModifiedImportedVariants.get_for_pipeline(upload_pipeline)
    vcf_import_annotate_dir = upload_pipeline.get_pipeline_processing_subdir("vcf_import_annotate")
    sort_order = upload_pipeline.get_max_step_sort_order()
    for split_vcf_filename in glob.glob(f"{split_vcf_dir}/*.vcf.gz"):
        sort_order += 1
        separate_upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                         name="Separate Unknown Variants Task",
                                                         sort_order=sort_order,
                                                         task_type=UploadStepTaskType.CELERY,
                                                         pipeline_stage=VCFPipelineStage.PRE_DATA_INSERTION,
                                                         input_filename=split_vcf_filename)
        separate_upload_step.launch_task(SeparateUnknownVariantsTask)
        output_filename = split_vcf_filename

        # If we annotate, that file will be processed in UploadStepMultiFileOutput
        if annotate_gnomad_af:
            sort_order += 1
            name = name_from_filename(split_vcf_filename, remove_gz=True)
            output_filename = os.path.join(vcf_import_annotate_dir, f"{name}.annotated.vcf.gz")
            gnomad_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                    name="Annotate gnomAD AF",
                                                    sort_order=sort_order,
                                                    task_type=UploadStepTaskType.CELERY,
                                                    pipeline_stage=VCFPipelineStage.PRE_DATA_INSERTION,
                                                    input_filename=split_vcf_filename,
                                                    output_filename=output_filename)
            gnomad_step.launch_task(AnnotateImportedVCFTask)

        # We don't know how big the last split file is, so leave items_to_process as null so no check
        UploadStepMultiFileOutput.objects.create(upload_step=upload_step,
                                                 output_filename=output_filename)


def _store_vcf_stats(filename, upload_step, description, operation='Skipped'):
    if os.path.exists(filename):
        df = pd.read_csv(filename, header=None, sep='\t', index_col=0)
        if not df.empty:
            for name, count in df.iloc[:, 0].items():
                message_string = f"{operation} {count} '{name}' {description}"
                SimpleVCFImportInfo.objects.create(upload_step=upload_step, message_string=message_string)


def preprocess_gene_level_vcf(upload_step):
    """ Preprocess for gene-level variants - split into sub-files and nothing else.

        The bcftools stages every other import runs all need a reference sequence: 'bcftools norm
        --check-ref=s --fasta-ref' reads the base at the position, and vcf_clean_and_filter maps
        chroms through standard_contigs_only (which excludes the gene-level contig by role) before
        replacing the header with one built the same way. A gene-level locus has no reference base
        and its position is a gene id, so none of that applies - @see snpdb.gene_level_variants.

        Nothing is lost by skipping them: the VCF is written by us (GeneFusionCreateVCFTask) rather
        than supplied by a lab, so it already has our contig names, no undeclared FILTERs, one alt
        per record and is written in sorted order. Everything downstream of the split - unknown
        variant insert, the parallel data-insertion tasks, max-variant tracking - is unchanged. """

    vcf_filename = upload_step.input_filename
    if not os.path.exists(vcf_filename):
        raise FileNotFoundError(f"Can't access vcf: '{vcf_filename}'")

    upload_pipeline = upload_step.upload_pipeline
    split_vcf_dir = upload_pipeline.get_pipeline_processing_subdir("split_vcf")
    header_filename = _write_gene_level_header(upload_pipeline, vcf_filename)
    vcf_name = name_from_filename(vcf_filename, remove_gz=True)

    split_env = {**os.environ,
                 'VG_HEADER_FILE': header_filename,
                 'VG_SPLIT_VCF_DIR': split_vcf_dir}
    split_file_rows = upload_step.split_file_rows or settings.VCF_IMPORT_FILE_SPLIT_ROWS
    pipe_commands = {
        "cat": ["cat", vcf_filename],
        REMOVE_HEADER_SUB_STEP: [settings.BCFTOOLS_COMMAND, "view", "--no-header", "-"],
        SPLIT_VCF_SUB_STEP: ["split", "-", vcf_name, "--additional-suffix=.vcf.gz", "--numeric-suffixes",
                             "--lines", str(split_file_rows),
                             "--filter='bash -c \"set -eo pipefail; { cat $VG_HEADER_FILE; cat; } | bgzip -c > $VG_SPLIT_VCF_DIR/$FILE\"'"],
    }
    _run_pipe(pipe_commands, {}, split_env, upload_pipeline)
    _schedule_split_file_steps(upload_step, split_vcf_dir)


def _write_gene_level_header(upload_pipeline, vcf_filename) -> str:
    """ The header prepended to each split file. Taken from the VCF as written rather than rebuilt
        from the build's contigs, which would drop the gene-level contig on its SequenceRole. """
    clean_vcf_dir = upload_pipeline.get_pipeline_processing_subdir("clean_vcf_dir")
    header_filename = os.path.join(clean_vcf_dir,
                                   name_from_filename(vcf_filename, remove_gz=True) + ".header.vcf")
    with open(vcf_filename) as in_f, open(header_filename, "w") as out_f:
        for line in in_f:
            if not line.startswith("#"):
                break
            out_f.write(line)
    return header_filename
