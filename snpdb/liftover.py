"""
    Liftover: convert variants to other genome builds
"""
import itertools
import logging
import os
from collections import defaultdict
from collections.abc import Iterable
from functools import lru_cache
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Prefetch, QuerySet

from genes.hgvs import HGVSMatcher
from library.django_utils.django_file_utils import get_import_processing_dir
from library.genomics.vcf_utils import get_contigs_header_lines, write_vcf_from_variant_coordinates
from library.guardian_utils import admin_bot
from snpdb.bcftools_liftover import bcftools_pre_liftover_error_check
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models.models_enums import (
    AlleleConversionTool,
    AlleleOrigin,
    ImportSource,
    ProcessingStatus,
)
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import (
    Allele,
    AlleleLiftover,
    LiftoverRun,
    Variant,
    VariantAllele,
    VariantCoordinate,
)
from upload.models import FileUpload, UploadedFileTypes, UploadedLiftover, UploadPipeline
from upload.upload_processing import process_upload_pipeline

# VariantCoordinate can be None, with an error string at the end
LIFTOVER_TUPLE = tuple[Optional[AlleleConversionTool], Optional['VariantCoordinate'], Optional[str]]

def create_liftover_pipelines(user: User, alleles: Iterable[Allele],
                              import_source: ImportSource,
                              inserted_genome_build: GenomeBuild,
                              destination_genome_builds: list[GenomeBuild] = None):
    """ Creates and runs a liftover pipeline for each destination GenomeBuild (default = all other builds)

        Alleles are handled in batches of settings.LIFTOVER_BATCH_SIZE - a batch's alleles, coordinates and
        AlleleLiftover records are all held in memory while its VCF is written, and anything that goes wrong
        only takes out the batch rather than the whole run """

    for allele_batch in _batch_alleles(alleles):
        _create_liftover_pipelines_for_batch(user, allele_batch, import_source,
                                             inserted_genome_build, destination_genome_builds)


def _batch_alleles(alleles: Iterable[Allele]) -> Iterable[list[Allele]]:
    """ QuerySets are paged by pk rather than iterated, so we don't hold a cursor open for the
        (potentially hours) it takes to create the pipelines. Re-querying also skips alleles that
        earlier batches have since lifted over """
    batch_size = settings.LIFTOVER_BATCH_SIZE
    if isinstance(alleles, QuerySet):
        allele_qs = alleles.order_by("pk")
        last_pk = 0
        while allele_batch := list(allele_qs.filter(pk__gt=last_pk)[:batch_size]):
            last_pk = allele_batch[-1].pk
            yield allele_batch
    else:
        for allele_batch in itertools.batched(alleles, batch_size):
            yield list(allele_batch)


def _create_liftover_pipelines_for_batch(user: User, alleles: list[Allele],
                                         import_source: ImportSource,
                                         inserted_genome_build: GenomeBuild,
                                         destination_genome_builds: list[GenomeBuild] = None):
    build_liftover_existing_allele_and_variants, build_liftover_allele_variant_coordinate_error = _get_build_liftover_dicts(alleles, inserted_genome_build, destination_genome_builds)
    for genome_build, liftover_tuples in build_liftover_existing_allele_and_variants.items():
        for conversion_tool, av_tuples in liftover_tuples.items():
            liftover = LiftoverRun.objects.create(user=user,
                                                  conversion_tool=conversion_tool,
                                                  genome_build=genome_build)

            if conversion_tool == AlleleConversionTool.SAME_CONTIG:
                _run_liftover_using_same_contig(liftover, av_tuples)

    for genome_build, liftover_tuples in build_liftover_allele_variant_coordinate_error.items():
        for conversion_tool, allele_variant_coordinate_error in liftover_tuples.items():
            liftover = LiftoverRun.objects.create(user=user,
                                                  conversion_tool=conversion_tool,
                                                  genome_build=genome_build)
            # Because we need to normalise / insert etc, it's easier just to write a VCF
            # and run through upload pipeline
            working_dir = get_import_processing_dir(liftover.pk, "liftover")
            liftover_vcf_filename = os.path.join(working_dir, f"liftover_variants.{genome_build.name}.vcf")
            if AlleleConversionTool.vcf_tuples_in_destination_build(conversion_tool):
                vcf_genome_build = genome_build
                vcf_filename = liftover_vcf_filename  # Can write directly
            else:
                vcf_genome_build = inserted_genome_build
                vcf_filename = os.path.join(working_dir, f"source_variants.{inserted_genome_build.name}.vcf")
                liftover.source_vcf = vcf_filename
                liftover.source_genome_build = inserted_genome_build
                liftover.save()

            allele_liftover_records = []
            # These 2 are in sync, used to write VCF
            vcf_ids = []
            variant_coordinates = []
            for allele, variant_coordinate, error_message in allele_variant_coordinate_error:
                if variant_coordinate is not None:
                    if contig_error := _non_standard_contig_error(vcf_genome_build, variant_coordinate):
                        variant_coordinate = None
                        error_message = contig_error

                if variant_coordinate is not None:
                    al = AlleleLiftover(allele=allele,
                                        liftover=liftover,
                                        status=ProcessingStatus.CREATED)
                    vcf_ids.append(allele.pk)
                    variant_coordinates.append(variant_coordinate)
                else:
                    al = AlleleLiftover(allele=allele,
                                        liftover=liftover,
                                        status=ProcessingStatus.ERROR,
                                        error={
                                            AlleleLiftover.ERROR_JSON_MESSAGE_KEY: f"Validation Error: {error_message}",
                                        })
                allele_liftover_records.append(al)

            if allele_liftover_records:
                AlleleLiftover.objects.bulk_create(allele_liftover_records, batch_size=2000)

            if vcf_ids:  # Need to write VCF and run
                # BCFTools uses chromosomes not contigs
                used_chroms = set(vc.chrom for vc in variant_coordinates)
                header_lines = get_contigs_header_lines(vcf_genome_build, use_accession=False,
                                                        contig_allow_list=used_chroms)
                write_vcf_from_variant_coordinates(vcf_filename, variant_coordinates=variant_coordinates,
                                                   vcf_ids=vcf_ids, header_lines=header_lines)
                file_upload = FileUpload.objects.create(path=liftover_vcf_filename,
                                                        import_source=import_source,
                                                        name='Liftover',
                                                        user=user,
                                                        file_type=UploadedFileTypes.LIFTOVER)

                UploadedLiftover.objects.create(file_upload=file_upload,
                                                liftover=liftover)
                upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
                process_upload_pipeline(upload_pipeline)
            else:
                logging.info("LiftoverRun %s doesn't need to be run", liftover)


def _non_standard_contig_error(vcf_genome_build: GenomeBuild, variant_coordinate: VariantCoordinate) -> Optional[str]:
    """ Liftover VCFs are written with standard contigs only (and tools like BCFTools look the chrom up in the
        reference fasta) so an alt/unlocalized contig would fail the entire run - see issue #1197 """
    if variant_coordinate.chrom in vcf_genome_build.chrom_standard_contig_mappings:
        return None

    if contig := vcf_genome_build.chrom_contig_mappings.get(variant_coordinate.chrom):
        description = f"'{contig}' has role '{contig.get_role_display()}'"
    else:
        description = f"'{variant_coordinate.chrom}' is not in {vcf_genome_build}"
    return f"Liftover VCFs only contain standard contigs - {description}"


def _liftover_allele_qs(allele_ids: list[int]) -> QuerySet[Allele]:
    """ Everything the per-allele liftover methods below look at, fetched up front for the whole batch """
    variant_allele_qs = VariantAllele.objects.select_related("genome_build", "variant__locus__contig",
                                                            "variant__locus__ref",
                                                            "variant__alt").order_by("genome_build__name")
    return Allele.objects.filter(pk__in=allele_ids).select_related("clingen_allele") \
        .prefetch_related(Prefetch("variantallele_set", queryset=variant_allele_qs))


@lru_cache
def _genome_build_contig_ids(genome_build: GenomeBuild) -> frozenset[int]:
    """ Cached as it's checked once per allele """
    return frozenset(genome_build.contig_ids)


def _variant_allele_for_build(allele, genome_build: GenomeBuild) -> Optional['VariantAllele']:
    """ Uses the prefetch (@see _liftover_allele_qs) rather than a query per allele """
    for variant_allele in allele.variantallele_set.all():
        if variant_allele.genome_build_id == genome_build.pk:
            return variant_allele
    return None


def _get_build_liftover_dicts(alleles: Iterable[Allele], inserted_genome_build: GenomeBuild,
                              destination_genome_builds: list[GenomeBuild] = None) -> tuple[dict, dict]:
    """ ID column set to allele_id """
    if destination_genome_builds is None:
        destination_genome_builds = GenomeBuild.builds_with_annotation()

    other_builds = set()
    other_build_contig_ids = set()
    for genome_build in destination_genome_builds:
        if genome_build != inserted_genome_build:
            other_builds.add(genome_build)
            other_build_contig_ids.update(genome_build.contig_ids)

    if not other_builds:
        return {}, {}  # Nothing to do

    allele_ids = [allele.pk for allele in alleles]
    alleles = _liftover_allele_qs(allele_ids)
    build_failed_tools = {gb: AlleleLiftover.get_failed_conversion_tools(allele_ids, gb) for gb in other_builds}

    build_liftover_existing_allele_and_variants = defaultdict(lambda: defaultdict(list))  # Already lifted over
    build_liftover_allele_variant_coordinate_error = defaultdict(lambda: defaultdict(list))  # Need to run pipelines

    for allele in alleles:
        # Builds this allele already has variants for
        existing_builds = {va.genome_build_id for va in allele.variantallele_set.all()
                           if va.variant.locus.contig_id in other_build_contig_ids}
        for genome_build in other_builds:
            if genome_build.pk in existing_builds:
                # logging.info("%s already lifted over to %s", allele, genome_build)
                continue

            # Now try different liftover methods (same contig, using dest build coords, using other builds and tool)
            # Then exit loop w/continue on first success

            same_contig_tool, variant = _liftover_using_existing_contig(allele, genome_build)
            if same_contig_tool:
                build_liftover_existing_allele_and_variants[genome_build][same_contig_tool].append((allele, variant))
                continue

            hgvs_matcher = HGVSMatcher.instance(genome_build)
            failed_tools = build_failed_tools[genome_build].get(allele.pk, set())

            for tool_coordinate_error in itertools.chain(
                    _liftover_using_dest_variant_coordinate(allele, genome_build,
                                                            hgvs_matcher=hgvs_matcher,
                                                            failed_tools=failed_tools),
                    _liftover_using_source_variant_coordinate(allele, inserted_genome_build, genome_build,
                                                              failed_tools=failed_tools),
            ):
                conversion_tool, variant_coordinate, error_string = tool_coordinate_error
                allele_vc_err = (allele, variant_coordinate, error_string)
                build_liftover_allele_variant_coordinate_error[genome_build][conversion_tool].append(allele_vc_err)
                if variant_coordinate:
                    # Quit on 1st success
                    break

    return build_liftover_existing_allele_and_variants, build_liftover_allele_variant_coordinate_error


def liftover_alleles(allele_qs, user: User = None):
    """ Creates then runs (async) liftover pipelines for a queryset of alleles """
    if user is None:
        user = admin_bot()

    for genome_build in GenomeBuild.builds_with_annotation():
        variants_qs = Variant.objects.filter(variantallele__allele__in=allele_qs)
        populate_clingen_alleles_for_variants(genome_build, variants_qs)
        create_liftover_pipelines(user, allele_qs, ImportSource.WEB, inserted_genome_build=genome_build)


def _run_liftover_using_same_contig(liftover, av_tuples: list[tuple[Allele, Variant]]):
    """ Special case of e.g. Mitochondria that has the same contig across multiple builds
        we just need to create a VariantAllele object - will already have annotation for both builds """

    variant_alleles = []
    allele_liftovers = []
    for allele, variant in av_tuples:
        va = VariantAllele(variant=variant,
                           genome_build=liftover.genome_build,
                           allele=allele,
                           origin=AlleleOrigin.LIFTOVER,
                           allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
        variant_alleles.append(va)

        al = AlleleLiftover(allele=allele,
                            liftover=liftover,
                            status=ProcessingStatus.SUCCESS)
        allele_liftovers.append(al)

    if variant_alleles:
        VariantAllele.objects.bulk_create(variant_alleles, ignore_conflicts=True, batch_size=2000)

    if allele_liftovers:
        AlleleLiftover.objects.bulk_create(allele_liftovers, batch_size=2000)


def _liftover_using_existing_contig(allele, dest_genome_build: GenomeBuild) -> tuple[AlleleConversionTool, 'Variant']:
    """ For Mito, 37 and 38 contigs are the same so we can re-use a variant """
    conversion_tool = None
    variant = None

    # Check if the other build shares existing contig and the variant already exists
    genome_build_contigs = _genome_build_contig_ids(dest_genome_build)
    # We shouldn't be here if a variant for build is already linked to allele - don't return these
    for variant_allele in allele.variantallele_set.all():
        if variant_allele.genome_build_id == dest_genome_build.pk:
            continue
        if variant_allele.variant.locus.contig_id in genome_build_contigs:
            conversion_tool = AlleleConversionTool.SAME_CONTIG
            # Return variant_id so we can create it directly
            variant = variant_allele.variant
    return conversion_tool, variant


def _liftover_using_dest_variant_coordinate(allele, dest_genome_build: GenomeBuild,
                                            hgvs_matcher=None,
                                            failed_tools: set[str] = None) -> Iterable[LIFTOVER_TUPLE]:
    """ This returns tuples FOR a genome build (if something can look them up)

        Used by to write VCF coordinates during liftover. Can be slow (API call)

        If you know a VariantAllele exists for your build, use variant_for_build(genome_build).as_tuple()

        Optionally pass in hgvs_matcher to save re-instantiating it all the time, and failed_tools
        (@see AlleleLiftover.get_failed_conversion_tools) to save a query per allele """

    from annotation.models import VariantAnnotationVersion
    from genes.hgvs import get_hgvs_variant_coordinate
    from snpdb.models.models_dbsnp import DbSNP

    if failed_tools is None:
        failed_tools = AlleleLiftover.get_failed_conversion_tools([allele.pk], dest_genome_build).get(allele.pk, set())

    conversion_tool = None
    g_hgvs = None
    if AlleleConversionTool.CLINGEN_ALLELE_REGISTRY not in failed_tools:
        clingen_failure_message = None
        if allele.clingen_allele:
            try:
                g_hgvs = allele.clingen_allele.get_g_hgvs(dest_genome_build)
                conversion_tool = AlleleConversionTool.CLINGEN_ALLELE_REGISTRY
            except ValueError:  # Various contig errors all subclass from this
                clingen_failure_message = f"{allele.clingen_allele} did not contain g.HGVS for {dest_genome_build}"
        else:
            # allele.clingen_allele is None — ClinGen was skipped, not tried-and-failed
            # (a real failure would be stored via AlleleLiftover, caught by failed_tools above)
            skip_reason = None
            if va := next(iter(allele.variantallele_set.all()), None):
                skip_reason = va.variant.clingen_allele_skip_reason()
            clingen_failure_message = skip_reason or "ClinGen skipped: no ClinGen allele registered for this variant"

        # Store the fact that we couldn't use ClinGen
        if clingen_failure_message:
            yield AlleleConversionTool.CLINGEN_ALLELE_REGISTRY, None, clingen_failure_message

    if g_hgvs is None:
        if settings.LIFTOVER_DBSNP_ENABLED:
            if AlleleConversionTool.DBSNP not in failed_tools:
                conversion_tool = AlleleConversionTool.DBSNP
                error_message = None
                if va := next(iter(allele.variantallele_set.all()), None):
                    if dbsnp := DbSNP.get_for_variant(va.variant, VariantAnnotationVersion.latest(va.genome_build)):
                        g_hgvs = dbsnp.get_g_hgvs(dest_genome_build, alt=va.variant.alt)
                    else:
                        error_message = f"No dbSNP for {va.variant} ({va.genome_build})"
                else:
                    error_message = "Allele contains no VariantAlleles at all! Cannot liftover"
                if error_message:
                    yield conversion_tool, None, error_message

    if g_hgvs:
        if hgvs_matcher:
            variant_coordinate = hgvs_matcher.get_variant_coordinate(g_hgvs)
        else:
            variant_coordinate = get_hgvs_variant_coordinate(g_hgvs, dest_genome_build)

        yield conversion_tool, variant_coordinate, None


def _liftover_using_source_variant_coordinate(allele, source_genome_build: GenomeBuild,
                                             dest_genome_build: GenomeBuild,
                                             failed_tools: set[str] = None) -> Iterable[LIFTOVER_TUPLE]:
    """ This gets tuples from another build to run through a tool

        Optionally pass in failed_tools (@see AlleleLiftover.get_failed_conversion_tools) to save a query """

    # If we have >=3 builds, sometimes this will be called for us when we don't have the variant in this build
    variant_allele = _variant_allele_for_build(allele, source_genome_build)
    if not variant_allele:
        return
    variant_coordinate = variant_allele.variant.coordinate

    if failed_tools is None:
        failed_tools = AlleleLiftover.get_failed_conversion_tools([allele.pk], dest_genome_build).get(allele.pk, set())

    # BCFTools fails with "Unable to fetch sequence" if any variant is outside contig size
    variant_errors = Variant.validate(source_genome_build, variant_coordinate.chrom,
                                      variant_coordinate.position)
    variant_errors_str = "\n".join(variant_errors)

    # Try tools that write other builds, then run conversion
    options = [
        # Enabled, Tool Enum, Requires reference to match genome build
        (settings.LIFTOVER_BCFTOOLS_ENABLED, AlleleConversionTool.BCFTOOLS_LIFTOVER, bcftools_pre_liftover_error_check),
    ]

    for enabled, conversion_tool, check_liftover_errors in options:
        if enabled:
            if conversion_tool in failed_tools:
                continue  # Skip as already failed liftover method to desired build

            if variant_errors_str:
                yield conversion_tool, None, variant_errors_str
                continue

            if check_liftover_errors is not None:
                if error_message := check_liftover_errors(variant_coordinate, source_genome_build):
                    yield conversion_tool, None, error_message
                    continue

            yield conversion_tool, variant_coordinate, None
            break  # Just want 1st one


def allele_can_attempt_liftover(allele, genome_build) -> bool:
    conversion_tool, variant = _liftover_using_existing_contig(allele, genome_build)
    if conversion_tool and variant:
        return True

    for conversion_tool, variant_coordinate, _error_message in _liftover_using_dest_variant_coordinate(allele, genome_build):
        if conversion_tool and variant_coordinate:
            return True

    for va in allele.variantallele_set.all():
        if va.genome_build_id == genome_build.pk:
            continue
        for conversion_tool, variant_coordinate, _error_message in _liftover_using_source_variant_coordinate(allele, va.genome_build, genome_build):
            if conversion_tool and variant_coordinate:
                return True

    return False
