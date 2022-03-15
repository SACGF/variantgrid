import glob
import logging
import os
from collections import OrderedDict
from subprocess import Popen, PIPE, CalledProcessError

import cyvcf2
import pandas as pd
from django.conf import settings
from django.utils import timezone

from library.django_utils.django_file_utils import get_import_processing_filename
from library.file_utils import name_from_filename
from upload.models import ModifiedImportedVariants, ToolVersion, UploadStep, \
    UploadStepTaskType, VCFSkippedContigs, VCFSkippedContig, UploadStepMultiFileOutput, VCFPipelineStage, \
    SimpleVCFImportInfo
from upload.tasks.vcf.unknown_variants_task import SeparateUnknownVariantsTask, AnnotateImportedVCFTask


def get_vt_tool_version(vt_command):
    p = Popen([vt_command, "--version"], stderr=PIPE)
    _, stderr = p.communicate()
    vt_version = stderr.decode().split('\n')[0]
    tool_version, _ = ToolVersion.objects.get_or_create(name='vt', version=vt_version)
    return tool_version


def create_sub_step(upload_step, sub_step_name, sub_step_commands):
    start_date = timezone.now()
    tool_version = get_vt_tool_version(settings.VCF_IMPORT_VT_COMMAND)
    command_line = ' '.join(sub_step_commands)
    return UploadStep.objects.create(name=sub_step_name,
                                     script=command_line,
                                     upload_pipeline=upload_step.upload_pipeline,
                                     sort_order=upload_step.sort_order,
                                     parent_upload_step=upload_step,
                                     tool_version=tool_version,
                                     start_date=start_date,
                                     task_type=UploadStepTaskType.TOOL)


def _write_split_headers(vcf_filename, split_headers_filename):
    """ Remove any INFO lines (they can break bcftools annotate) and write in VT normalize ones """
    VT_HEADERS = [
        '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">',
        '##INFO=<ID=OLD_VARIANT,Number=.,Type=String,Description="Original chr:pos:ref:alt encoding">',
    ]
    with open(split_headers_filename, "w") as f:
        written_vt_headers = False
        for line in cyvcf2.Reader(vcf_filename).raw_header.split("\n"):
            info_line = line.startswith("##INFO")
            if info_line or line.startswith("#CHROM"):
                if not written_vt_headers:
                    for vt_line in VT_HEADERS:
                        f.write(vt_line + "\n")
                    written_vt_headers = True

                if info_line:
                    continue
            if line:
                f.write(line + "\n")


def preprocess_vcf(upload_step, annotate_gnomad_af=False):
    MAX_STDERR_OUTPUT = 5000  # How much stderr output per process to store in DB

    VCF_CLEAN_AND_FILTER_SUB_STEP = "vcf_clean_and_filter"
    DECOMPOSE_SUB_STEP = "decompose"
    NORMALIZE_SUB_STEP = "normalize"
    UNIQ_SUB_STEP = "uniq"
    REMOVE_HEADER_SUB_STEP = "remove_header"
    SPLIT_VCF_SUB_STEP = "split_vcf"

    genome_build = upload_step.genome_build
    _ = genome_build.reference_fasta  # Fails if not available.

    vcf_filename = upload_step.input_filename
    if not os.path.exists(vcf_filename):
        raise FileNotFoundError(f"Can't access vcf: '{vcf_filename}'")

    pipe_commands = OrderedDict()  # Insert in pipeline order
    sub_steps = {}
    if vcf_filename.endswith(".gz"):
        pipe_commands["zcat"] = [settings.BASH_ZCAT, vcf_filename]
    else:
        pipe_commands["cat"] = ["cat", vcf_filename]

    upload_pipeline = upload_step.upload_pipeline
    vcf_name = name_from_filename(vcf_filename, remove_gz=True)
    skipped_contigs_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.skipped_contigs.tsv")
    skipped_records_stats_file = get_import_processing_filename(upload_pipeline.pk, f"{vcf_name}.skipped_records.tsv")

    manage_command = settings.MANAGE_COMMAND
    read_variants_cmd = manage_command + ["vcf_clean_and_filter",
                                          "--genome-build",
                                          genome_build.name,
                                          "--skipped-contigs-stats-file",
                                          skipped_contigs_stats_file,
                                          "--skipped-records-stats-file",
                                          skipped_records_stats_file]
    pipe_commands[VCF_CLEAN_AND_FILTER_SUB_STEP] = read_variants_cmd
    sub_steps[VCF_CLEAN_AND_FILTER_SUB_STEP] = create_sub_step(upload_step, VCF_CLEAN_AND_FILTER_SUB_STEP, read_variants_cmd)

    # VT isn't the bottleneck here, it's my programs - so no speed advantage to using "+" for Uncompressed BCF streams
    pipe_commands[DECOMPOSE_SUB_STEP] = [settings.VCF_IMPORT_VT_COMMAND, "decompose", "-s", "-"]
    pipe_commands[NORMALIZE_SUB_STEP] = [settings.VCF_IMPORT_VT_COMMAND, "normalize", "-n", "-r", genome_build.reference_fasta, "-"]
    pipe_commands[UNIQ_SUB_STEP] = [settings.VCF_IMPORT_VT_COMMAND, "uniq", "-"]

    # Split up the VCF
    split_vcf_dir = upload_pipeline.get_pipeline_processing_subdir("split_vcf")
    # We'll cat split_headers_filename on top of every split VCF
    split_headers_filename = os.path.join(split_vcf_dir, name_from_filename(vcf_filename, remove_gz=True))
    _write_split_headers(vcf_filename, split_headers_filename)

    pipe_commands[REMOVE_HEADER_SUB_STEP] = [settings.VCF_IMPORT_VT_COMMAND, "view", "-"]
    pipe_commands[SPLIT_VCF_SUB_STEP] = ["split", "-", vcf_name, "--additional-suffix=.vcf.gz", "--numeric-suffixes",
                                         "--lines", str(settings.VCF_IMPORT_FILE_SPLIT_ROWS),
                                         f"--filter='sh -c \"{{ cat {split_headers_filename}; cat; }} | bgzip -c > {split_vcf_dir}/$FILE\"'"]

    for sub_step_name in [DECOMPOSE_SUB_STEP, NORMALIZE_SUB_STEP, UNIQ_SUB_STEP]:
        sub_step_commands = pipe_commands[sub_step_name]
        sub_steps[sub_step_name] = create_sub_step(upload_step, sub_step_name, sub_step_commands)

    piped_command = ' | '.join([' '.join(command_with_args) for command_with_args in pipe_commands.values()])

    # if POPEN_SHELL=True we're running this all as one command as native shell text
    # this is not recommended, less portable and if commands have errors
    if settings.POPEN_SHELL:
        print("single_commands: %s" % " | ".join([' '.join(x) for x in pipe_commands.values()]))
        p = Popen(
            piped_command,
            shell=True
        )
        _, p_stderr = p.communicate()
        if p_stderr:
            logging.error("%s stderr:", piped_command)
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

            for contig, count in df.iloc[:, 0].iteritems():
                VCFSkippedContig.objects.create(import_info=import_info,
                                                contig=contig,
                                                num_skipped=count)
    if os.path.exists(skipped_records_stats_file):
        df = pd.read_csv(skipped_records_stats_file, header=None, sep='\t', index_col=0)
        if not df.empty:
            for name, count in df.iloc[:, 0].iteritems():
                message_string = f"Skipped {count} '{name}' records."
                SimpleVCFImportInfo.objects.create(upload_step=clean_sub_step, message_string=message_string)

    # Create this here so downstream tasks can add modified imported variant messages
    import_info, _ = ModifiedImportedVariants.objects.get_or_create(upload_step=upload_step)
    vcf_import_annotate_dir = upload_pipeline.get_pipeline_processing_subdir("vcf_import_annotate")
    sort_order = upload_pipeline.get_max_step_sort_order()
    for split_vcf_filename in glob.glob(f"{split_vcf_dir}/*.vcf.gz"):
        sort_order += 1
        separate_upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                         name="Separate Unknown Variants Task",
                                                         sort_order=sort_order,
                                                         task_type=UploadStepTaskType.CELERY,
                                                         pipeline_stage=VCFPipelineStage.INSERT_UNKNOWN_VARIANTS,
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
                                                    pipeline_stage=VCFPipelineStage.INSERT_UNKNOWN_VARIANTS,
                                                    input_filename=split_vcf_filename,
                                                    output_filename=output_filename)
            gnomad_step.launch_task(AnnotateImportedVCFTask)

        # We don't know how big the last split file is, so leave items_to_process as null so no check
        UploadStepMultiFileOutput.objects.create(upload_step=upload_step,
                                                 output_filename=output_filename)
