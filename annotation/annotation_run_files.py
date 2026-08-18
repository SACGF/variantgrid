"""
An AnnotationRun's on-disk artifacts: where each attempt writes, and how its variants get dumped.

Split out of annotation.tasks.annotate_variants so the pipeline runners (annotation.pipelines) and the
cleanup receiver (annotation.signals.annotation_run_cleanup) can share it without importing the celery
tasks - the tasks import the runners, so the dependency has to run the other way.
"""
import gzip
import os

from django.conf import settings

from library.utils.file_utils import name_from_filename
from snpdb.variants_to_vcf import VARIANT_GRID_INFO_DICT, write_contig_sorted_values_to_vcf_file


def get_annotated_filename(annotation_run, vcf_dump_filename) -> str:
    """ Path VEP writes its annotated VCF to for a given dump. Derived from the dump stem, which #1658
        makes per-task, so each attempt's annotated output is a file of its own. """
    name = name_from_filename(vcf_dump_filename)
    vcf_annotated_basename = f"{name}.vep_annotated_{annotation_run.genome_build.name}.vcf.gz"
    return os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, vcf_annotated_basename)


def get_annotsv_dir(annotation_run) -> str:
    """ AnnotSV output dir for a run. Deliberately keyed on the run, not the task - shared between
        attempts, so any party holding the run can name it. See the #720 note: the TSV inside is named
        from the per-task dump stem, so concurrent attempts write side by side and neither can truncate
        the other, while a per-task *directory* would be nameable only by the attempt that created it. """
    return os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, f"annotsv_{annotation_run.pk}")


def write_qs_to_vcf(vcf_filename, genome_build, qs, info_dict=VARIANT_GRID_INFO_DICT, use_accession=False,
                    samples=None) -> int:
    # We had an issue with writing accessions in VEP, so use chrom names and the default VEP fasta instead
    # @see https://github.com/Ensembl/ensembl-vep/issues/1635
    # Contigs are shared between builds (eg GRCh37/hg19) so the ordering join needs restricting to this
    # build, otherwise a variant is written once per build its contig belongs to
    qs = qs.filter(locus__contig__genomebuildcontig__genome_build=genome_build)
    qs = qs.order_by("locus__contig__genomebuildcontig__order", "locus__position")
    if use_accession:
        chrom_key = "locus__contig__refseq_accession"
    else:
        chrom_key = "locus__contig__name"

    sorted_values = qs.values("id", chrom_key, "locus__position",
                              "locus__ref__seq", "alt__seq", "end", "svlen")

    # External dumps are written .vcf.gz (see AnnotationRun.get_dump_filename) - gzip on the way out
    if vcf_filename.endswith(".gz"):
        f = gzip.open(vcf_filename, "wt", compresslevel=6)
    else:
        f = open(vcf_filename, "w")
    with f:
        return write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict=info_dict,
                                                      use_accession=use_accession, samples=samples)
