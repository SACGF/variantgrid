"""
Loader for Illumina DRAGEN TSO 500's AllFusions.csv - one vendor's format, not a standard. Anything
here that reads a named column belongs to that format; the fusion identity it resolves to does not
(@see genes.gene_fusions).

The rows become a VCF of gene-level variants which goes through the normal VCF import pipeline, so
the VCF/Sample/Cohort come from the header the way every other import's do, and the CohortGenotype
rows are written by the same SQL COPY path. Only the bcftools stages are skipped, since they all
need a reference base a gene-level locus does not have. Nothing here names a genome build - the
file's '# Source =' line becomes '##source' and VCFSourceSettings says what that caller is run
against (@see upload.vcf.vcf_import.resolve_genome_build).
@see snpdb.gene_level_variants for why these are Variants at all, and
upload.vcf.gene_level_vcf_preprocess for exactly what is skipped and why.

Two steps:
  DragenTSO500AllFusionsCreateVCFTask - parse the CSV, resolve gene names, write the VCF
  DragenTSO500AllFusionsInsertTask    - once the variants exist, create the GeneFusions against them

Each caller row becomes an observation carried in INFO, so what the caller wrote survives import.
Several rows can name one gene pair (one caller reports ENTPD3::RPL14 three times with three 5'
breakpoints), and those become one Variant with several observations.
"""
import logging
from collections import defaultdict

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
from upload.tso500.dragen_all_fusions_parser import read_all_fusions
from upload.models import (
    ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
    ModifiedImportedVariants,
    UploadStep,
)
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from variantgrid.celery import app

# INFO fields carrying what the caller reported - these land in CohortGenotype.info via the
# standard bulk importer, which stores every INFO field the header declares
FUSION_INFO = "FUSION"
FUSION_OBSERVATIONS_INFO = "FUSION_OBS"

# A fusion caller asserts the fusion is present, not a diploid genotype, so there is no zygosity to
# report. The importer reads a missing GT as unknown zygosity.
NO_GENOTYPE_CALL = "./."


def _source_from_comments(comments) -> str:
    """ eg '# Source = FusionProcessor 1.0.0.614' - the same thing '##source' gives a VCF """
    for comment in comments:
        key, _, value = comment.lstrip("# ").partition("=")
        if key.strip() == "Source":
            return value.strip()
    return ""


def _observations_by_variant_coordinate(rows) -> dict:
    """ {variant coordinate: [the rows that named that gene pair, as the caller wrote them]} """
    resolver = GeneFusionResolver()
    observations = defaultdict(list)
    for row in rows:
        gene_a = resolver.resolve_side(row.gene_a)
        gene_b = resolver.resolve_side(row.gene_b)
        resolved_fusion = resolver.resolve_fusion(gene_a, gene_b, row.directionality_known)
        observations[resolved_fusion].append(row.data)
    return observations


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
        formats=['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'],
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
                                info=info, fmt="GT", sample_calls=[NO_GENOTYPE_CALL])


class DragenTSO500AllFusionsCreateVCFTask(ImportVCFStepTask):
    """ Write the fusion variants as a VCF so they go through the normal insert pipeline """

    def process_items(self, upload_step):
        comments, rows = read_all_fusions(upload_step.input_filename)
        observations = _observations_by_variant_coordinate(rows)
        file_upload = upload_step.upload_pipeline.file_upload
        _write_gene_level_vcf(upload_step.output_filename, observations,
                              sample_name=file_upload.name, source=_source_from_comments(comments))
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
            calls = "; ".join(f"{o.get('Caller')} {o.get('Gene A Breakpoint')}->{o.get('Gene B Breakpoint')}"
                              for o in observations)
            merged.append((variant_id, len(observations), calls))

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
