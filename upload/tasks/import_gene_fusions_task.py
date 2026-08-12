"""
Loader for AllFusions.csv.

The rows become a VCF of gene-level variants, which then goes through the normal VCF import
pipeline - so the variants are inserted by the same optimised path as everything else, and
max-variant tracking (and therefore annotation gating) comes for free. Only the bcftools stages are
skipped, since they all need a reference base a gene-level locus does not have.
@see snpdb.gene_level_variants for why these are Variants at all, and
upload.vcf.vcf_preprocess.preprocess_gene_level_vcf for exactly what is skipped and why.

Three steps, in pipeline order:
  GeneFusionCreateVCFTask      - parse the CSV, resolve gene names, write the VCF
  GeneFusionCreateVCFModelTask - create the VCF/Sample/Cohort the rows will hang off
  GeneFusionInsertTask         - once the variants exist, create GeneFusions and CohortGenotypes

The CSV is re-read by the last step rather than staged in a model of its own: resolution is
deterministic (the FusionGeneIds already exist by then), and it is a small file.
"""
import logging
from collections import Counter, defaultdict

from django.db import transaction

from genes.gene_fusions import GeneFusionResolver, ResolvedFusion
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models_enums import Zygosity
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_LENGTH, GENE_LEVEL_CONTIG_NAME
from snpdb.models import (
    Cohort,
    CohortGenotype,
    CohortSample,
    GenomeBuild,
    ImportStatus,
    Sample,
    Variant,
)
from snpdb.models.models_enums import CohortGenotypeCollectionType
from snpdb.tasks.cohort_genotype_tasks import create_cohort_genotype_collection
from upload.gene_fusions.all_fusions_parser import read_all_fusions
from upload.models import (
    ModifiedImportedVariant,
    ModifiedImportedVariantOperation,
    ModifiedImportedVariants,
    UploadedVCF,
    UploadStep,
)
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.upload_metadata import get_metadata_genome_build, get_metadata_source
from upload.vcf.vcf_import import assign_sample_extractions, create_vcf_from_uploaded_vcf
from variantgrid.celery import app


def _resolve_genome_build(file_upload) -> GenomeBuild:
    """ The file declares no build - its breakpoints are chrN:pos with nothing to detect from. The
        fusion Variant itself is build-independent (one shared contig), so this is only for the
        VCF/Sample bookkeeping and the VariantAllele. """
    if genome_build := get_metadata_genome_build(file_upload):
        return genome_build
    genome_builds = list(GenomeBuild.builds_with_annotation())
    if len(genome_builds) == 1:
        return genome_builds[0]
    raise ValueError("AllFusions.csv declares no genome build - send one as upload metadata "
                     f"('genome_build'). Available: {GenomeBuild.available_names_or_aliases()}")


def _source_from_comments(comments) -> str:
    """ eg '# Source = FusionProcessor 1.0.0.614' - the same thing '##source' gives a VCF """
    for comment in comments:
        key, _, value = comment.lstrip("# ").partition("=")
        if key.strip() == "Source":
            return value.strip()
    return ""


def _resolve_rows(rows) -> list[tuple[ResolvedFusion, dict]]:
    """ (identity, the row as the caller wrote it) per call. Deterministic, so the create-VCF step
        and the insert step arrive at the same identities """
    resolver = GeneFusionResolver()
    resolved = []
    for row in rows:
        gene_a = resolver.resolve_side(row.gene_a)
        gene_b = resolver.resolve_side(row.gene_b)
        resolved.append((resolver.resolve_fusion(gene_a, gene_b, row.directionality_known), row.data))
    return resolved


def _write_gene_level_vcf(filename: str, variant_coordinates):
    """ Written already-clean and sorted, which is what lets preprocess skip straight to the split.
        END = POS gives svlen 0 through vcf_get_ref_alt_svlen_and_modification, which needs one of
        SVLEN/END for any symbolic alt (and reads SVLEN=0 as absent). """
    header = [
        "##fileformat=VCFv4.2",
        "##source=VariantGrid",
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">',
        f"##contig=<ID={GENE_LEVEL_CONTIG_NAME},length={GENE_LEVEL_CONTIG_LENGTH}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    with open(filename, "w") as f:
        for line in header:
            f.write(line + "\n")
        for vc in sorted(variant_coordinates, key=lambda c: (c.position, c.alt)):
            f.write(f"{vc.chrom}\t{vc.position}\t.\t{vc.ref}\t{vc.alt}\t.\t.\tEND={vc.position}\n")


class GeneFusionCreateVCFTask(ImportVCFStepTask):
    """ Write the fusion variants as a VCF so they go through the normal insert pipeline """

    def process_items(self, upload_step):
        _comments, rows = read_all_fusions(upload_step.input_filename)
        resolved = _resolve_rows(rows)
        # One record per distinct gene pair - several calls can share one, and they become one Variant
        variant_coordinates = {r.variant_coordinate for r, _data in resolved}
        _write_gene_level_vcf(upload_step.output_filename, variant_coordinates)
        return len(rows)


class GeneFusionCreateVCFModelTask(ImportVCFStepTask):
    """ The VCF/Sample/Cohort the rows hang off. Runs where a normal import reads its VCF header -
        this file has no header worth reading, so the build comes from upload metadata instead.

        The file creates its own VCF and Sample rather than joining the RNA arm's SpliceVariants
        sample: a Sample belongs to exactly one VCF, and CohortGenotype hangs off a per-VCF
        collection. Phase 4's extraction linking ties this sample back to the rest of the arm. """

    def process_items(self, upload_step):
        upload_pipeline = upload_step.upload_pipeline
        file_upload = upload_pipeline.file_upload
        comments, _rows = read_all_fusions(file_upload.get_filename())

        uploaded_vcf, _ = UploadedVCF.objects.update_or_create(
            file_upload=file_upload,
            defaults={"upload_pipeline": upload_pipeline, "vcf_importer": None})

        vcf = create_vcf_from_uploaded_vcf(uploaded_vcf, num_genotype_samples=0)
        vcf.genome_build = _resolve_genome_build(file_upload)
        vcf.source = get_metadata_source(file_upload) or _source_from_comments(comments)
        vcf.header = "\n".join(comments)
        vcf.save()
        uploaded_vcf.vcf = vcf
        uploaded_vcf.save()

        sample = Sample.objects.create(vcf=vcf,
                                       vcf_sample_name=file_upload.name,
                                       name=file_upload.name,
                                       import_status=ImportStatus.IMPORTING)
        assign_sample_extractions(vcf, upload_step)

        cohort = Cohort.objects.create(vcf=vcf,
                                       name=f"VCF: {vcf.pk} {vcf}"[:40],
                                       import_status=ImportStatus.IMPORTING,
                                       genome_build=vcf.genome_build,
                                       user=vcf.user)
        assign_permission_to_user_and_groups(vcf.user, cohort)
        CohortSample.objects.create(cohort=cohort, sample=sample,
                                    cohort_genotype_packed_field_index=0, sort_order=0)
        create_cohort_genotype_collection(cohort)
        return 0


class GeneFusionInsertTask(ImportVCFStepTask):
    """ Runs after data insertion, so the Variants exist - create the GeneFusions against them and
        the CohortGenotype rows carrying what each caller actually reported """

    def process_items(self, upload_step: UploadStep):
        file_upload = upload_step.upload_pipeline.file_upload
        _comments, rows = read_all_fusions(file_upload.get_filename())
        resolved = _resolve_rows(rows)

        vcf = file_upload.uploadedvcf.vcf
        sample = vcf.sample_set.get()
        collection = vcf.cohort.cohortgenotypecollection_set.get(
            collection_type=CohortGenotypeCollectionType.UNCOMMON)

        results = Counter()
        observations_by_variant = defaultdict(list)
        fusion_by_variant = {}
        for resolved_fusion, row_data in resolved:
            variant = self._get_variant(resolved_fusion, vcf.genome_build)
            gene_fusion = resolved_fusion.get_gene_fusion(variant)
            fusion_by_variant[variant.pk] = gene_fusion
            observations_by_variant[variant.pk].append(row_data)
            results["rows"] += 1

        merged = {}
        with transaction.atomic():
            for variant_id, observations in observations_by_variant.items():
                gene_fusion = fusion_by_variant[variant_id]
                CohortGenotype.objects.update_or_create(
                    collection=collection,
                    variant_id=variant_id,
                    defaults={
                        # A fusion caller asserts presence, not a diploid genotype, so there's no
                        # zygosity to record - same as the sibling SpliceVariants.vcf, which has no GT
                        "samples_zygosity": Zygosity.UNKNOWN_ZYGOSITY,
                        "info": {
                            "fusion": gene_fusion.canonical_str,
                            "observations": observations,
                        },
                    })
                results["variants"] += 1
                if len(observations) > 1:
                    merged[gene_fusion.variant] = [
                        f"{o.get('Caller')} {o.get('Gene A Breakpoint')}->{o.get('Gene B Breakpoint')}"
                        for o in observations
                    ]

            _record_merged_rows(upload_step, merged)

            sample.import_status = ImportStatus.SUCCESS
            sample.save()
            vcf.import_status = ImportStatus.SUCCESS
            vcf.save()
            cohort = vcf.cohort
            cohort.import_status = ImportStatus.SUCCESS
            cohort.save()

        logging.info("Imported gene fusions: %s", dict(results))
        return results["rows"]

    @staticmethod
    def _get_variant(resolved_fusion: ResolvedFusion, genome_build: GenomeBuild) -> Variant:
        variant_coordinate = resolved_fusion.variant_coordinate
        try:
            return Variant.get_from_variant_coordinate(variant_coordinate, genome_build)
        except Variant.DoesNotExist as dne:
            raise ValueError(f"Gene fusion variant '{variant_coordinate}' was not inserted by the "
                             f"import pipeline") from dne


def _record_merged_rows(upload_step, merged: dict):
    """ Several rows can share one gene pair (one caller reports ENTPD3-RPL14 three times with three
        5' breakpoints) and become one Variant, so one CohortGenotype. Every row's data is kept in
        the info blob; this records that the merge happened, the way SHARED_LOCUS does for VCF rows
        whose depths were summed. """
    if not merged:
        return
    import_info = ModifiedImportedVariants.get_for_pipeline(upload_step.upload_pipeline)
    ModifiedImportedVariant.objects.bulk_create([
        ModifiedImportedVariant(import_info=import_info,
                                variant=variant,
                                operation=ModifiedImportedVariantOperation.SHARED_LOCUS,
                                operation_detail=f"{len(observations)} calls merged onto one gene pair: "
                                                 + "; ".join(observations))
        for variant, observations in merged.items()
    ])


GeneFusionCreateVCFTask = app.register_task(GeneFusionCreateVCFTask())
GeneFusionCreateVCFModelTask = app.register_task(GeneFusionCreateVCFModelTask())
GeneFusionInsertTask = app.register_task(GeneFusionInsertTask())
