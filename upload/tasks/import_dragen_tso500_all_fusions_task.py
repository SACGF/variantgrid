"""
Loader for Illumina DRAGEN TSO 500's AllFusions.csv - one vendor's format, not a standard. Anything
here that reads a named column belongs to that format; the fusion identity it resolves to does not
(@see genes.gene_fusions).

The rows become a VCF of gene-level variants which goes through the normal VCF import pipeline, so
the VCF/Sample/Cohort come from the header the way every other import's do, and the CohortGenotype
rows are written by the same SQL COPY path. Only the bcftools stages are skipped, since they all
need a reference base a gene-level locus does not have. Nothing here names a genome build - the
file's '# Source =' line becomes '##source' and VCFSourceSettings says what that caller is run
against (@see upload.vcf.vcf_import.resolve_genome_build). The create-VCF step resolves that build
itself, since the breakpoints it resolves gene names against are positions in one.
@see snpdb.gene_level_variants for why these are Variants at all, and
upload.vcf.gene_level_vcf_preprocess for exactly what is skipped and why.

Two steps:
  DragenTSO500AllFusionsCreateVCFTask - parse the CSV, resolve gene names, write the VCF
  DragenTSO500AllFusionsInsertTask    - once the variants exist, create the GeneFusions against them

Each caller row becomes an observation carried in INFO, so what the caller wrote survives import.
Several rows can name one gene pair (one caller reports ENTPD3::RPL14 three times with three 5'
breakpoints), and those become one Variant with several observations.

A fusion caller asserts the fusion is present, not a diploid genotype, so the VCF has no GT and the
sample has no zygosity to filter on. What it does have is read support, written as the sample's
ALT_READS and REF_READS; the ^FusionProcessor VCFSourceSettings row binds those as alt and ref depth
so the sample node's minimum-reads threshold and allele frequency work on fusions.
"""
import logging
from collections import defaultdict
from typing import Optional

import simplejson

from genes.gene_fusions import GeneFusionResolver, create_gene_fusions_for_variants
from library.genomics.vcf_writer import (
    VCFInfoHeader,
    VCFWriter,
    build_header_lines,
    percent_decode_info_value,
    percent_encode_info_value,
)
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_LENGTH, GENE_LEVEL_CONTIG_NAME
from snpdb.models import GenomeBuild
from upload.models import (
    ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
    ModifiedImportedVariants,
    UploadStep,
)
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.tso500.dragen_all_fusions_parser import (
    FUSION_INFO,
    FUSION_OBSERVATIONS_INFO,
    format_fusion_observations,
    read_all_fusions,
    reference_reads,
    supporting_reads,
)
from upload.vcf.vcf_import import resolve_genome_build_from_source
from variantgrid.celery import app

# The sample's FORMAT fields - read support rather than a genotype
ALT_READS_FORMAT = "ALT_READS"
REF_READS_FORMAT = "REF_READS"
VCF_MISSING_VALUE = "."


def _source_from_comments(comments) -> str:
    """ eg '# Source = FusionProcessor 1.0.0.614' - the same thing '##source' gives a VCF """
    for comment in comments:
        key, _, value = comment.lstrip("# ").partition("=")
        if key.strip() == "Source":
            return value.strip()
    return ""


def _observations_by_variant_coordinate(rows, genome_build: Optional[GenomeBuild]) -> dict:
    """ {variant coordinate: [the rows that named that gene pair, as the caller wrote them]}

        The breakpoints decide which gene each side is where we know the build to look them up in;
        with no build resolvable the names are all there is, which is where this started """
    resolver = GeneFusionResolver()
    observations = defaultdict(list)
    for row in rows:
        gene_a = resolver.resolve_side(row.gene_a, breakpoint=row.gene_a_breakpoint,
                                       genome_build=genome_build)
        gene_b = resolver.resolve_side(row.gene_b, breakpoint=row.gene_b_breakpoint,
                                       genome_build=genome_build)
        resolved_fusion = resolver.resolve_fusion(gene_a, gene_b, row.directionality_known)
        observations[resolved_fusion].append(row.data)
    return observations


def _read_support(observations: list[dict]) -> tuple[Optional[int], Optional[int]]:
    """ (supporting reads, reference reads) for one gene pair. Each observation is a distinct
        breakpoint with its own supporting reads, so those add up; the reference reads across a
        junction are re-reported by every call that shares it, so the largest stands for the pair """
    alt = [r for o in observations if (r := supporting_reads(o)) is not None]
    ref = [r for o in observations if (r := reference_reads(o)) is not None]
    return (sum(alt) if alt else None), (max(ref) if ref else None)


def _sample_call(observations: list[dict]) -> str:
    alt, ref = _read_support(observations)
    return ":".join(VCF_MISSING_VALUE if v is None else str(v) for v in (alt, ref))


def _write_gene_level_vcf(filename: str, observations: dict, sample_name: str, source: str):
    """ Written already-clean and sorted, which is what lets preprocess skip straight to the split.
        END = POS gives svlen 0 through vcf_get_ref_alt_svlen_and_modification, which needs one of
        SVLEN/END for any symbolic alt (and reads SVLEN=0 as absent). """

    meta_lines = [f"##source={source}"] if source else []
    header_lines = build_header_lines(
        meta_lines=meta_lines,
        info=[
            VCFInfoHeader(id="END", type="Integer", description="Stop position of the interval"),
            VCFInfoHeader(id=FUSION_INFO, type="String",
                          description="Gene fusion in VICC gene-level nomenclature"),
            VCFInfoHeader(id=FUSION_OBSERVATIONS_INFO, type="String",
                          description="JSON list of the caller rows this fusion was called from"),
        ],
        formats=[
            f'##FORMAT=<ID={ALT_READS_FORMAT},Number=1,Type=Integer,'
            f'Description="Reads supporting the fusion, summed over its breakpoints">',
            f'##FORMAT=<ID={REF_READS_FORMAT},Number=1,Type=Integer,'
            f'Description="Reads across the junctions that do not support the fusion">',
        ],
        contig_lines=[f"##contig=<ID={GENE_LEVEL_CONTIG_NAME},length={GENE_LEVEL_CONTIG_LENGTH}>"],
        samples=[sample_name],
    )

    with open(filename, "w") as f:
        writer = VCFWriter(f, header_lines, encode_info=percent_encode_info_value)
        for resolved_fusion in sorted(observations, key=lambda r: (r.anchor.pk, r.alt)):
            variant_coordinate = resolved_fusion.variant_coordinate
            info = {
                "END": variant_coordinate.position,
                FUSION_INFO: resolved_fusion.canonical_str,
                FUSION_OBSERVATIONS_INFO: simplejson.dumps(observations[resolved_fusion], ignore_nan=True),
            }
            writer.write_record(variant_coordinate.chrom, variant_coordinate.position,
                                variant_coordinate.ref, variant_coordinate.alt,
                                info=info, fmt=f"{ALT_READS_FORMAT}:{REF_READS_FORMAT}",
                                sample_calls=[_sample_call(observations[resolved_fusion])])


class DragenTSO500AllFusionsCreateVCFTask(ImportVCFStepTask):
    """ Write the fusion variants as a VCF so they go through the normal insert pipeline """

    def process_items(self, upload_step):
        comments, rows = read_all_fusions(upload_step.input_filename)
        source = _source_from_comments(comments)
        file_upload = upload_step.upload_pipeline.file_upload
        # The VCF this step writes is what the build would normally be resolved from, so the source
        # line has to answer it here - @see upload.vcf.vcf_import.resolve_genome_build
        genome_build = resolve_genome_build_from_source(source, file_upload)
        observations = _observations_by_variant_coordinate(rows, genome_build)
        _write_gene_level_vcf(upload_step.output_filename, observations,
                              sample_name=file_upload.name, source=source)
        return len(rows)


class DragenTSO500AllFusionsInsertTask(ImportVCFStepTask):
    """ Runs after data insertion, so the Variants and their CohortGenotypes exist - create the
        GeneFusions against them """

    def process_items(self, upload_step: UploadStep):
        vcf = upload_step.upload_pipeline.uploadedvcf.vcf
        variant_qs = vcf.get_variant_qs()
        created = create_gene_fusions_for_variants(variant_qs)
        _record_merged_rows(upload_step, variant_qs)
        logging.info("Created %d gene fusions for %s", created, vcf)
        return created


def _record_merged_rows(upload_step, variant_qs):
    """ Several rows can share one gene pair and become one Variant, so one CohortGenotype. Every
        row's data is kept in the info blob; this records that the merge happened. """

    cgc = upload_step.upload_pipeline.uploadedvcf.vcf.cohort.cohort_genotype_collection
    info_alias = f"{cgc.cohortgenotype_alias}__info"
    merged = []
    for variant_id, info in variant_qs.values_list("pk", info_alias):
        # INFO values are stored as VCF writes them - htslib doesn't decode, so we do
        encoded = (info or {}).get(FUSION_OBSERVATIONS_INFO) or "[]"
        observations = simplejson.loads(percent_decode_info_value(encoded))
        if len(observations) > 1:
            merged.append((variant_id, len(observations), format_fusion_observations(observations)))

    if not merged:
        return

    import_info = ModifiedImportedVariants.get_for_pipeline(upload_step.upload_pipeline)
    ModifiedImportedVariant.objects.bulk_create([
        ModifiedImportedVariant(import_info=import_info,
                                variant_id=variant_id,
                                operation=ModifiedImportedVariantOperation.MERGED_RECORDS,
                                operation_detail=f"{count} calls merged onto one gene pair: {calls}")
        for variant_id, count, calls in merged
    ])


DragenTSO500AllFusionsCreateVCFTask = app.register_task(DragenTSO500AllFusionsCreateVCFTask())
DragenTSO500AllFusionsInsertTask = app.register_task(DragenTSO500AllFusionsInsertTask())
