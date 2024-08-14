"""
    Liftover: convert variants to other genome builds
"""
import logging
import operator
import os
from collections import defaultdict
from functools import reduce
from typing import Iterable, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models.query_utils import Q

from genes.hgvs import HGVSMatcher
from library.django_utils.django_file_utils import get_import_processing_dir
from library.genomics.vcf_utils import write_vcf_from_variant_coordinates, get_contigs_header_lines
from library.guardian_utils import admin_bot
from library.log_utils import log_traceback
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models.models_enums import ImportSource, AlleleConversionTool, AlleleOrigin, ProcessingStatus
from snpdb.models.models_genome import GenomeBuild, Contig, GenomeFasta
from snpdb.models.models_variant import LiftoverRun, Allele, Variant, VariantAllele, AlleleLiftover
from upload.models import UploadedFile, UploadedLiftover, UploadPipeline, UploadedFileTypes
from upload.upload_processing import process_upload_pipeline

LIFTOVER_TOOL_AND_COORDINATE = tuple[Optional[AlleleConversionTool], Optional['VariantCoordinate']]


def create_liftover_pipelines(user: User, alleles: Iterable[Allele],
                              import_source: ImportSource,
                              inserted_genome_build: GenomeBuild,
                              destination_genome_builds: list[GenomeBuild] = None):
    """ Creates and runs a liftover pipeline for each destination GenomeBuild (default = all other builds) """

    build_liftover_existing_allele_and_variants, build_liftover_variant_coordinates = _get_build_liftover_dicts(alleles, inserted_genome_build, destination_genome_builds)
    for genome_build, liftover_tuples in build_liftover_existing_allele_and_variants.items():
        for conversion_tool, av_tuples in liftover_tuples.items():
            liftover = LiftoverRun.objects.create(user=user,
                                                  conversion_tool=conversion_tool,
                                                  genome_build=genome_build)

            if conversion_tool == AlleleConversionTool.SAME_CONTIG:
                _run_liftover_using_same_contig(liftover, av_tuples)

    for genome_build, liftover_tuples in build_liftover_variant_coordinates.items():
        for conversion_tool, allele_and_variant_coordinates in liftover_tuples.items():
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
            used_chroms: set[str] = set()
            for allele, variant_coordinate in allele_and_variant_coordinates:
                used_chroms.add(variant_coordinate.chrom)
                al = AlleleLiftover(allele=allele,
                                    liftover=liftover,
                                    status=ProcessingStatus.CREATED)
                allele_liftover_records.append(al)

            if allele_liftover_records:
                AlleleLiftover.objects.bulk_create(allele_liftover_records, batch_size=2000)

            # BCFTools uses chromosomes not contigs
            header_lines = get_contigs_header_lines(vcf_genome_build, use_accession=False, contig_allow_list=used_chroms)
            vcf_ids = (allele.pk for allele, _ in allele_and_variant_coordinates)
            variant_coordinates = (vc for _, vc in allele_and_variant_coordinates)
            write_vcf_from_variant_coordinates(vcf_filename, variant_coordinates=variant_coordinates,
                                               vcf_ids=vcf_ids, header_lines=header_lines)
            uploaded_file = UploadedFile.objects.create(path=liftover_vcf_filename,
                                                        import_source=import_source,
                                                        name='Liftover',
                                                        user=user,
                                                        file_type=UploadedFileTypes.LIFTOVER)

            UploadedLiftover.objects.create(uploaded_file=uploaded_file,
                                            liftover=liftover)
            upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
            process_upload_pipeline(upload_pipeline)


def _get_build_liftover_dicts(alleles: Iterable[Allele], inserted_genome_build: GenomeBuild,
                              destination_genome_builds: list[GenomeBuild] = None) -> tuple[dict, dict]:
    """ ID column set to allele_id """
    if destination_genome_builds is None:
        destination_genome_builds = GenomeBuild.builds_with_annotation()

    other_build_contigs_q_list = []
    other_builds = set()
    hgvs_matchers = {}
    for genome_build in destination_genome_builds:
        if genome_build != inserted_genome_build:
            other_builds.add(genome_build)
            q = Q(variantallele__variant__locus__contig__in=genome_build.contigs)
            other_build_contigs_q_list.append(q)
        hgvs_matchers[genome_build] = HGVSMatcher(genome_build)

    if not other_builds:
        return {}, {}  # Nothing to do

    other_build_contigs_q = reduce(operator.or_, other_build_contigs_q_list)

    # Store builds where alleles have already been lifted over
    allele_builds = defaultdict(set)
    allele_ids = [allele.pk for allele in alleles]
    qs = Allele.objects.filter(other_build_contigs_q, pk__in=allele_ids)
    for allele_id, genome_build_name in qs.values_list("pk", "variantallele__genome_build"):
        allele_builds[allele_id].add(genome_build_name)

    build_liftover_existing_allele_and_variants = defaultdict(lambda: defaultdict(list))  # Already lifted over
    build_liftover_variant_coordinates = defaultdict(lambda: defaultdict(list))  # Need to run pipelines

    for allele in alleles:
        existing_builds = allele_builds[allele.pk]
        for genome_build in other_builds:
            if genome_build.pk in existing_builds:
                logging.info("%s already lifted over to %s", allele, genome_build)
                continue

            hgvs_matcher = hgvs_matchers[genome_build]
            conversion_tool = None
            variant_coordinate = None
            try:
                conversion_tool, variant = _liftover_using_existing_contig(allele, genome_build)
                if conversion_tool:
                    build_liftover_existing_allele_and_variants[genome_build][conversion_tool].append((allele, variant))
                else:
                    # Try and get coordinates for builds we want
                    conversion_tool, variant_coordinate = _liftover_using_dest_variant_coordinate(allele, genome_build,
                                                                                                  hgvs_matcher=hgvs_matcher)
            except (Contig.ContigNotInBuildError, GenomeFasta.ContigNotInFastaError):
                log_traceback()

            if not variant_coordinate:
                conversion_tool, variant_coordinate = _liftover_using_source_variant_coordinate(allele, inserted_genome_build, genome_build)

            if conversion_tool and variant_coordinate:
                build_liftover_variant_coordinates[genome_build][conversion_tool].append((allele, variant_coordinate))

    return build_liftover_existing_allele_and_variants, build_liftover_variant_coordinates


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
    genome_build_contigs = set(c.pk for c in dest_genome_build.chrom_contig_mappings.values())
    # We shouldn't be here if a variant for build is already linked to allele - don't return these
    for variant_allele in allele.variantallele_set.exclude(genome_build=dest_genome_build):
        if variant_allele.variant.locus.contig_id in genome_build_contigs:
            conversion_tool = AlleleConversionTool.SAME_CONTIG
            # Return variant_id so we can create it directly
            variant = variant_allele.variant
    return conversion_tool, variant


def _liftover_using_dest_variant_coordinate(allele, dest_genome_build: GenomeBuild,
                                            hgvs_matcher=None) -> LIFTOVER_TOOL_AND_COORDINATE:
    """ This returns tuples FOR a genome build (if something can look them up)

        Used by to write VCF coordinates during liftover. Can be slow (API call)

        If you know a VariantAllele exists for your build, use variant_for_build(genome_build).as_tuple()

        Optionally pass in hgvs_matcher to save re-instantiating it all the time """

    from annotation.models import VariantAnnotationVersion
    from snpdb.models.models_dbsnp import DbSNP
    from genes.hgvs import get_hgvs_variant_coordinate

    conversion_tool = None
    g_hgvs = None
    clingen_failure_message = None
    if allele.clingen_allele:
        try:
            g_hgvs = allele.clingen_allele.get_g_hgvs(dest_genome_build)
            conversion_tool = AlleleConversionTool.CLINGEN_ALLELE_REGISTRY
        except ValueError:  # Various contig errors all subclass from this
            clingen_failure_message = f"{allele.clingen_allele} did not contain g.HGVS for {dest_genome_build}"
    else:
        clingen_failure_message = f"No ClinGenAllele for variant"

    # Store the fact that we couldn't use ClinGen
    if clingen_failure_message:
        lr = LiftoverRun.get_clingen_auto_fail_liftover_run(dest_genome_build)
        AlleleLiftover.objects.create(liftover=lr, allele=allele, status=ProcessingStatus.ERROR,
                                      error={"message": clingen_failure_message})

    if g_hgvs is None:
        if settings.LIFTOVER_DBSNP_ENABLED:
            va = allele.variantallele_set.all().first()
            if va is None:
                raise ValueError("Allele contains no VariantAlleles at all! Cannot liftover")
            if dbsnp := DbSNP.get_for_variant(va.variant, VariantAnnotationVersion.latest(va.genome_build)):
                g_hgvs = dbsnp.get_g_hgvs(dest_genome_build, alt=va.variant.alt)
                conversion_tool = AlleleConversionTool.DBSNP

    variant_coordinate = None
    if g_hgvs:
        if hgvs_matcher:
            variant_coordinate = hgvs_matcher.get_variant_coordinate(g_hgvs)
        else:
            variant_coordinate = get_hgvs_variant_coordinate(g_hgvs, dest_genome_build)

    return conversion_tool, variant_coordinate


def _liftover_using_source_variant_coordinate(allele, source_genome_build: GenomeBuild,
                                             dest_genome_build: GenomeBuild) -> LIFTOVER_TOOL_AND_COORDINATE:
    """ This gets tuples from another build to run through a tool """

    # Try tools that write other builds, then run conversion
    options = [
        (settings.LIFTOVER_BCFTOOLS_ENABLED, AlleleConversionTool.BCFTOOLS_LIFTOVER),
    ]

    conversion_tool = None
    variant_coordinate = None
    for enabled, potential_conversion_tool in options:
        if enabled:
            if allele.alleleliftover_set.filter(liftover__genome_build=dest_genome_build,
                                                liftover__conversion_tool=potential_conversion_tool,
                                                status=ProcessingStatus.ERROR).exists():
                continue  # Skip as already failed liftover method to desired build
            conversion_tool = potential_conversion_tool
            break  # Just want 1st one

    if conversion_tool:
        # Return VCF tuples in inserted genome build
        try:
            variant = allele.variant_for_build(source_genome_build)
            variant_coordinate = variant.coordinate
            # BCFTools fails with "Unable to fetch sequence" if any variant is outside contig size
            if errors := Variant.validate(source_genome_build, variant_coordinate.chrom, variant_coordinate.position):
                raise ValueError("\n".join(errors))
        except ValueError as e:  # No variant for source build (merged allele?)
            logging.warning("Skipped %s: %s", allele, e)
            return None, None

    return conversion_tool, variant_coordinate


def allele_can_attempt_liftover(allele, genome_build) -> bool:
    try:
        conversion_tool, _ = allele.liftover_using_existing_variant(genome_build)
        if conversion_tool is not None:
            return True

        conversion_tool, _ = allele.liftover_using_dest_variant_coordinate(genome_build)[0]
        if conversion_tool is not None:
            return True
    except:
        pass

    for va in allele.variantallele_set.exclude(genome_build=genome_build):
        conversion_tool, _ = _liftover_using_source_variant_coordinate(allele, va.genome_build, genome_build)
        if conversion_tool is not None:
            return True

    return False
