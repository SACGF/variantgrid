"""
An AnnotationRun's on-disk artifacts: where each attempt writes, and how its variants get dumped.

Split out of annotation.tasks.annotate_variants so the pipeline runners (annotation.pipelines) and the
cleanup receiver (annotation.signals.annotation_run_cleanup) can share it without importing the celery
tasks - the tasks import the runners, so the dependency has to run the other way.
"""
import io
import os
from typing import Optional

from bgzip import BGZipWriter
from django.conf import settings

from annotation.annotation_pipeline_routing import symbolic_annotated_as_small
from library.utils.file_utils import name_from_filename
from snpdb.models.models_variant import VariantCoordinate
from snpdb.variants_to_vcf import VARIANT_GRID_INFO_DICT, write_contig_sorted_values_to_vcf_file

# Prefix of a run's scratch dir under settings.IMPORT_PROCESSING_DIR (where the bulk inserter writes the
# CSVs it SQL COPYs from). Here rather than on the inserter so AnnotationRun can name the dir it owns
# without importing the annotation pipeline.
ANNOTATION_RUN_IMPORT_PROCESSING_PREFIX = "annotation_run"


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


def _small_symbolic_as_explicit(genome_build, sorted_values, structural_variant_min_size: int):
    """ A symbolic del/dup/inv the STANDARD pipeline annotates goes to VEP as its sequence, so VEP and its
        plugins match it as the small variant it is (@see annotation.annotation_pipeline_routing) """
    for data in sorted_values:
        alt = data["alt__seq"]
        if symbolic_annotated_as_small(alt, data["svlen"], structural_variant_min_size):
            vc = VariantCoordinate(chrom=data["locus__contig__name"], position=data["locus__position"],
                                   ref=data["locus__ref__seq"], alt=alt, svlen=data["svlen"])
            explicit = vc.as_external_explicit(genome_build)
            data = {**data, "locus__position": explicit.position, "locus__ref__seq": explicit.ref,
                    "alt__seq": explicit.alt, "end": None, "svlen": None}
        yield data


def write_qs_to_vcf(vcf_filename, genome_build, qs, info_dict=VARIANT_GRID_INFO_DICT, use_accession=False,
                    samples=None, structural_variant_min_size: Optional[int] = None) -> int:
    """ structural_variant_min_size: write a symbolic del/dup/inv shorter than this as its sequence """
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

    # Streamed (server side cursor) rather than materialised - a whole-version dump (#1675) is millions
    # of rows, and each is written and forgotten
    sorted_values = qs.values("id", chrom_key, "locus__position",
                              "locus__ref__seq", "alt__seq", "end", "svlen").iterator(chunk_size=10_000)
    if structural_variant_min_size:
        if use_accession:
            raise ValueError("Writing symbolic as explicit reads the reference by contig name, not accession")
        sorted_values = _small_symbolic_as_explicit(genome_build, sorted_values, structural_variant_min_size)

    # External dumps are written .vcf.gz (see AnnotationRun.get_dump_filename). bgzip rather than plain
    # gzip - reads anywhere gzip did, and tabix can index it, which is what lets a dump be the target of
    # a `bcftools annotate` (@see annotation.backfill_columns)
    if vcf_filename.endswith(".gz"):
        with open(vcf_filename, "wb") as raw:
            with BGZipWriter(raw) as bgzip_f:
                # bgzip is binary; wrap as text so VCFWriter only ever deals with str
                f = io.TextIOWrapper(bgzip_f, encoding="utf-8", write_through=True)
                try:
                    return write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f,
                                                                  info_dict=info_dict,
                                                                  use_accession=use_accession, samples=samples)
                finally:
                    # flush + detach so the BGZipWriter (closed by its 'with') is closed exactly once
                    f.flush()
                    f.detach()

    with open(vcf_filename, "w") as f:
        return write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict=info_dict,
                                                      use_accession=use_accession, samples=samples)
