import glob
import logging
import os
from subprocess import Popen, PIPE, CalledProcessError

import pandas as pd
from django.conf import settings
from django.utils import timezone

from library.django_utils.django_file_utils import get_import_processing_filename
from library.genomics.vcf_enums import INFO_LIFTOVER_SWAPPED_REF_ALT
from library.genomics.vcf_utils import write_cleaned_vcf_header
from library.utils.file_utils import name_from_filename
from upload.models import ModifiedImportedVariants, ToolVersion, UploadStep, \
    UploadStepTaskType, VCFSkippedContigs, VCFSkippedContig, UploadStepMultiFileOutput, VCFPipelineStage, \
    SimpleVCFImportInfo, ModifiedImportedVariant
from upload.tasks.vcf.unknown_variants_task import SeparateUnknownVariantsTask, AnnotateImportedVCFTask


def get_bcftools_tool_version(bcftools_command):
    p = Popen([bcftools_command, "--version"], stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    if p.returncode:
        raise CalledProcessError(f"Error running '{bcftools_command} --version': {stderr=}, returned: {p.returncode}")

    output_list = stdout.split("\n")
    bcftools_version = output_list[0]
    if not bcftools_version.startswith("bcftools"):
        raise CalledProcessError(f"Expected to find bcftools on 1st line of output of '{bcftools_command} --version' output: {bcftools_version}")
    htslib_version = output_list[1]
    if not "htslib" in htslib_version:
        raise CalledProcessError(f"Expected to find htslib version on 2nd line of '{bcftools_command} --version' output: {htslib_version}")

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
        f'##INFO=<ID={INFO_LIFTOVER_SWAPPED_REF_ALT},Number=0,Type=Flag,Description="Whether we have switched REF/ALT due to SWAP=1">'
    ]
    write_cleaned_vcf_header(genome_build, vcf_filename, cleaned_vcf_header_filename,
                             new_info_lines=HEADERS, standard_contigs_only=True)
    return cleaned_vcf_header_filename


def preprocess_vcf(upload_step, annotate_gnomad_af=False):
    MAX_STDERR_OUTPUT = 5000  # How much stderr output per process to store in DB

    VCF_CLEAN_AND_FILTER_SUB_STEP = "vcf_clean_and_filter"
    VCF_CLEAN_ALTS_SUB_STEP = "vcf_clean_alts"
    REMOVE_HEADER_SUB_STEP = "remove_header"
    SPLIT_VCF_SUB_STEP = "split_vcf"

    genome_build = upload_step.genome_build
    _ = genome_build.reference_fasta  # Fails if not available.

    vcf_filename = upload_step.input_filename
    if not os.path.exists(vcf_filename):
        raise FileNotFoundError(f"Can't access vcf: '{vcf_filename}'")

    pipe_commands = {}  # Use dict order to decide pipeline order
    sub_steps = {}
    if vcf_filename.endswith(".gz"):
        pipe_commands["zcat"] = [settings.BASH_ZCAT, vcf_filename]
    else:
        pipe_commands["cat"] = ["cat", vcf_filename]

    upload_pipeline = upload_step.upload_pipeline
    vcf_name = name_from_filename(vcf_filename, remove_gz=True)

    skipped_contigs_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.skipped_contigs.tsv")
    skipped_records_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.skipped_records.tsv")
    skipped_filters_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.skipped_filters.tsv")
    skipped_alts_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.skipped_alts.tsv")
    converted_alts_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.converted_alts.tsv")

    cleaned_vcf_header_filename = _write_cleaned_header(genome_build, upload_pipeline, vcf_filename)

    manage_command = settings.MANAGE_COMMAND
    read_variants_cmd = manage_command + ["vcf_clean_and_filter",
                                          "--genome-build",
                                          genome_build.name,
                                          "--replace-header",
                                          cleaned_vcf_header_filename,
                                          "--skipped-contigs-stats-file",
                                          skipped_contigs_stats_file,
                                          "--skipped-records-stats-file",
                                          skipped_records_stats_file,
                                          "--skipped-filters-stats-file",
                                          skipped_filters_stats_file]
    pipe_commands[VCF_CLEAN_AND_FILTER_SUB_STEP] = read_variants_cmd
    sub_steps[VCF_CLEAN_AND_FILTER_SUB_STEP] = create_sub_step(upload_step, VCF_CLEAN_AND_FILTER_SUB_STEP, read_variants_cmd)

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
    pipe_commands[VCF_CLEAN_ALTS_SUB_STEP] = manage_command + [
         "vcf_clean_alts",
         "--skipped-records-stats-file",
         skipped_alts_stats_file,
        "--converted-records-stats-file",
        converted_alts_stats_file,
    ]

    norm_substep_names = [UploadStep.NORMALIZE_SUB_STEP]

    # Split up the VCF
    split_file_rows = upload_step.split_file_rows or settings.VCF_IMPORT_FILE_SPLIT_ROWS
    split_vcf_dir = upload_pipeline.get_pipeline_processing_subdir("split_vcf")
    pipe_commands[SPLIT_VCF_SUB_STEP] = ["split", "-", vcf_name, "--additional-suffix=.vcf.gz", "--numeric-suffixes",
                                         "--lines", str(split_file_rows),
                                         f"--filter='sh -c \"{{ cat {cleaned_vcf_header_filename}; cat; }} | bgzip -c > {split_vcf_dir}/$FILE\"'"]

    for sub_step_name in norm_substep_names:
        sub_step_commands = pipe_commands[sub_step_name]
        sub_steps[sub_step_name] = create_sub_step(upload_step, sub_step_name, sub_step_commands)

    piped_command = ' | '.join([' '.join(command_with_args) for command_with_args in pipe_commands.values()])

    # if POPEN_SHELL=True we're running this all as one command as native shell text
    # this is not recommended, less portable and if commands have errors
    if settings.VCF_IMPORT_PREPROCESS_POPEN_SHELL:
        logging.info("single_commands: %s" % " | ".join([' '.join(x) for x in pipe_commands.values()]))
        p = Popen(
            piped_command,
            shell=True,
            stdout=PIPE,
            stderr=PIPE,
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
            raise CalledProcessError(p.returncode, piped_command, output=p_stderr)

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
            p = Popen(cmd, **kwargs)
            pipes[sub_step_name] = p
            stderr_filenames[sub_step_name] = stderr_filename

        print("pipe_commands: %s" % " | ".join([' '.join(x) for x in pipe_commands.values()]))

        for sub_step_name, p_cmd in reversed(pipe_commands.items()):
            p = pipes[sub_step_name]
            p.communicate()
            stderr_output = None
            stderr_filename = stderr_filenames.get(sub_step_name)
            if stderr_filename:
                with open(stderr_filename, "r") as sf:
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
                raise CalledProcessError(p.returncode, p_cmd, output=stderr_output)

    clean_sub_step = sub_steps[VCF_CLEAN_AND_FILTER_SUB_STEP]
    if os.path.exists(skipped_contigs_stats_file):
        df = pd.read_csv(skipped_contigs_stats_file, header=None, sep='\t', index_col=0)
        if not df.empty:
            import_info = VCFSkippedContigs.objects.create(upload_step=clean_sub_step)

            for contig, count in df.iloc[:, 0].items():
                VCFSkippedContig.objects.create(import_info=import_info,
                                                contig=contig,
                                                num_skipped=count)

    _store_vcf_stats(skipped_records_stats_file, clean_sub_step, "records")
    _store_vcf_stats(skipped_filters_stats_file, clean_sub_step, "FILTER")
    _store_vcf_stats(skipped_alts_stats_file, clean_sub_step, "ALTs")
    _store_vcf_stats(converted_alts_stats_file, clean_sub_step, "ALTs", operation='Converted')

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
