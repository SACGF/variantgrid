"""
Liftover via Clingen Allele Registry or NCBI remap
"""
import logging
import operator
import os
from collections import defaultdict
from functools import reduce
from typing import Any, Union, Iterable

from django.contrib.auth.models import User
from django.db.models.query_utils import Q

from genes.hgvs import HGVSMatcher
from library.django_utils.django_file_utils import get_import_processing_dir
from library.genomics.vcf_utils import write_vcf_from_tuples, get_vcf_header_contig_lines
from library.guardian_utils import admin_bot
from library.log_utils import log_traceback
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models.models_enums import ImportSource, AlleleConversionTool, AlleleOrigin, ProcessingStatus
from snpdb.models.models_genome import GenomeBuild, Contig, GenomeFasta
from snpdb.models.models_variant import LiftoverRun, Allele, Variant, VariantAllele, AlleleLiftover
from upload.models import UploadedFile, UploadedLiftover, UploadPipeline, UploadedFileTypes
from upload.upload_processing import process_upload_pipeline


def get_used_contigs_header_lines(genome_build, used_contigs: set) -> list[str]:
    contigs = []
    for contig in genome_build.contigs:
        if contig.name in used_contigs:
            contigs.append((contig.name, contig.length, genome_build.name))
    return get_vcf_header_contig_lines(contigs)


def create_liftover_pipelines(user: User, alleles: Iterable[Allele],
                              import_source: ImportSource,
                              inserted_genome_build: GenomeBuild,
                              destination_genome_builds: list[GenomeBuild] = None):
    """ Creates and runs a liftover pipeline for each destination GenomeBuild (default = all other builds) """

    build_liftover_vcf_tuples = _get_build_liftover_tuples(alleles, inserted_genome_build, destination_genome_builds)
    for genome_build, liftover_tuples in build_liftover_vcf_tuples.items():
        for conversion_tool, av_tuples in liftover_tuples.items():
            liftover = LiftoverRun.objects.create(user=user,
                                                  conversion_tool=conversion_tool,
                                                  genome_build=genome_build)

            if conversion_tool == AlleleConversionTool.SAME_CONTIG:
                _liftover_using_same_contig(liftover, av_tuples)
            else:
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
                used_contigs = set()
                for avt in av_tuples:
                    used_contigs.add(avt[0])
                    al = AlleleLiftover(allele_id=avt[2],
                                        liftover=liftover,
                                        status=ProcessingStatus.CREATED)
                    allele_liftover_records.append(al)

                if allele_liftover_records:
                    AlleleLiftover.objects.bulk_create(allele_liftover_records, batch_size=2000)

                header_lines = get_used_contigs_header_lines(vcf_genome_build, used_contigs)
                write_vcf_from_tuples(vcf_filename, av_tuples, tuples_have_id_field=True, header_lines=header_lines)
                uploaded_file = UploadedFile.objects.create(path=liftover_vcf_filename,
                                                            import_source=import_source,
                                                            name='Liftover',
                                                            user=user,
                                                            file_type=UploadedFileTypes.LIFTOVER)

                UploadedLiftover.objects.create(uploaded_file=uploaded_file,
                                                liftover=liftover)
                upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
                process_upload_pipeline(upload_pipeline)


VCF_ROW = tuple[str, int, int, str, str]
LIFTOVER_TUPLE = list[Union[tuple[int, int], VCF_ROW]]


def _get_build_liftover_tuples(alleles: Iterable[Allele], inserted_genome_build: GenomeBuild,
                               destination_genome_builds: list[GenomeBuild] = None) -> dict[Any, dict[Any, LIFTOVER_TUPLE]]:
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
        return {}  # Nothing to do

    other_build_contigs_q = reduce(operator.or_, other_build_contigs_q_list)

    # Store builds where alleles have already been lifted over
    allele_builds = defaultdict(set)
    allele_ids = [allele.pk for allele in alleles]
    qs = Allele.objects.filter(other_build_contigs_q, pk__in=allele_ids)
    for allele_id, genome_build_name in qs.values_list("pk", "variantallele__genome_build"):
        allele_builds[allele_id].add(genome_build_name)

    build_liftover_vcf_tuples = defaultdict(lambda: defaultdict(list))

    for allele in alleles:
        existing_builds = allele_builds[allele.pk]
        for genome_build in other_builds:
            if genome_build.pk in existing_builds:
                logging.info("%s already lifted over to %s", allele, genome_build)
                continue

            hgvs_matcher = hgvs_matchers[genome_build]
            conversion_tool = None
            variant_id_or_coordinate = None
            try:
                # Try and get coordinates for builds we want
                conversion_tool, variant_id_or_coordinate = allele.get_liftover_tuple(genome_build,
                                                                                      hgvs_matcher=hgvs_matcher)
            except (Contig.ContigNotInBuildError, GenomeFasta.ContigNotInFastaError):
                log_traceback()

            if variant_id_or_coordinate:
                # Converted ok - return VCF tuples in desired genome build
                if conversion_tool == AlleleConversionTool.SAME_CONTIG:
                    avt = (allele.pk, variant_id_or_coordinate)
                else:
                    chrom, position, ref, alt, svlen = variant_id_or_coordinate
                    avt = (chrom, position, allele.pk, ref, alt, svlen)
                build_liftover_vcf_tuples[genome_build][conversion_tool].append(avt)
            else:
                if other_build_data := allele.get_liftover_tuple_from_other_build(inserted_genome_build, genome_build):
                    conversion_tool, avt = other_build_data
                    build_liftover_vcf_tuples[genome_build][conversion_tool].append(avt)

    return build_liftover_vcf_tuples


def liftover_alleles(allele_qs, user: User = None):
    """ Creates then runs (async) liftover pipelines for a queryset of alleles """
    if user is None:
        user = admin_bot()

    for genome_build in GenomeBuild.builds_with_annotation():
        variants_qs = Variant.objects.filter(variantallele__allele__in=allele_qs)
        populate_clingen_alleles_for_variants(genome_build, variants_qs)
        create_liftover_pipelines(user, allele_qs, ImportSource.WEB, inserted_genome_build=genome_build)


def _liftover_using_same_contig(liftover, av_tuples: list[tuple[int, int]]):
    """ Special case of e.g. Mitochondria that has the same contig across multiple builds
        we just need to create a VariantAllele object - will already have annotation for both builds """

    variant_alleles = []
    allele_liftovers = []
    for allele_id, variant_id in av_tuples:
        va = VariantAllele(variant_id=variant_id,
                           genome_build=liftover.genome_build,
                           allele_id=allele_id,
                           origin=AlleleOrigin.LIFTOVER,
                           allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
        variant_alleles.append(va)

        al = AlleleLiftover(allele_id=allele_id,
                            liftover=liftover,
                            status=ProcessingStatus.SUCCESS)
        allele_liftovers.append(al)

    if variant_alleles:
        VariantAllele.objects.bulk_create(variant_alleles, ignore_conflicts=True, batch_size=2000)

    if allele_liftovers:
        AlleleLiftover.objects.bulk_create(allele_liftovers, batch_size=2000)
