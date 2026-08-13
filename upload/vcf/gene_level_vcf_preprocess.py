"""
Preprocess for gene-level variants - split the VCF into sub-files, and nothing else.

@see snpdb.gene_level_variants for why these are Variants at all.

Every stage upload.vcf.vcf_preprocess runs needs a reference sequence: 'bcftools norm --check-ref=s
--fasta-ref' reads the base at the position, and vcf_clean_and_filter maps chroms through
standard_contigs_only (which excludes the gene-level contig by role) before replacing the header with
one built the same way. A gene-level locus has no reference base and its position is a gene id, so
none of that applies.

Nothing is lost by skipping them: the VCF is written by us rather than supplied by a lab, so it
already has our contig names, no undeclared FILTERs, one alt per record, and is written sorted.
Everything downstream of the split - unknown variant insert, the parallel data-insertion tasks,
max-variant tracking - is the shared code path.
"""
import os

from django.conf import settings

from library.utils.file_utils import name_from_filename
from upload.vcf.vcf_preprocess import (
    REMOVE_HEADER_SUB_STEP,
    SPLIT_VCF_SUB_STEP,
    create_sub_step,
    run_pipe,
    schedule_split_file_steps,
)


def preprocess_gene_level_vcf(upload_step):
    vcf_filename = upload_step.input_filename
    if not os.path.exists(vcf_filename):
        raise FileNotFoundError(f"Can't access vcf: '{vcf_filename}'")

    upload_pipeline = upload_step.upload_pipeline
    split_vcf_dir = upload_pipeline.get_pipeline_processing_subdir("split_vcf")
    header_filename = _write_gene_level_header(upload_pipeline, vcf_filename)
    vcf_name = name_from_filename(vcf_filename, remove_gz=True)

    # Paths go via env vars rather than into the shell --filter string, so metacharacters in a path
    # cannot affect shell parsing - split runs the filter via sh -c, which expands them there.
    split_env = {**os.environ,
                 'VG_HEADER_FILE': header_filename,
                 'VG_SPLIT_VCF_DIR': split_vcf_dir}
    split_file_rows = upload_step.split_file_rows or settings.VCF_IMPORT_FILE_SPLIT_ROWS
    pipe_commands = {
        "cat": ["cat", vcf_filename],
        # sed rather than 'bcftools view --no-header' so no stage in this pipe needs a reference
        REMOVE_HEADER_SUB_STEP: ["sed", "/^#/d"],
        SPLIT_VCF_SUB_STEP: ["split", "-", vcf_name, "--additional-suffix=.vcf.gz", "--numeric-suffixes",
                             "--lines", str(split_file_rows),
                             "--filter='bash -c \"set -eo pipefail; { cat $VG_HEADER_FILE; cat; } | bgzip -c > $VG_SPLIT_VCF_DIR/$FILE\"'"],
    }
    # A sub step so a failed split reports on the pipeline the way the bcftools stages do. No tool
    # version - coreutils split is the only tool in this pipe
    sub_steps = {SPLIT_VCF_SUB_STEP: create_sub_step(upload_step, SPLIT_VCF_SUB_STEP,
                                                     pipe_commands[SPLIT_VCF_SUB_STEP],
                                                     tool_version=None)}
    run_pipe(pipe_commands, sub_steps, split_env, upload_pipeline)
    schedule_split_file_steps(upload_step, split_vcf_dir)


def _write_gene_level_header(upload_pipeline, vcf_filename) -> str:
    """ The header prepended to each split file. Taken from the VCF as written rather than rebuilt
        from the build's contigs, which would drop the gene-level contig on its SequenceRole. """
    header_dir = upload_pipeline.get_pipeline_processing_subdir("gene_level_header")
    header_filename = os.path.join(header_dir,
                                   name_from_filename(vcf_filename, remove_gz=True) + ".header.vcf")
    with open(vcf_filename) as in_f, open(header_filename, "w") as out_f:
        for line in in_f:
            if not line.startswith("#"):
                break
            out_f.write(line)
    return header_filename
