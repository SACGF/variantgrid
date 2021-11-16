"""
Liftover via Clingen Allele Registry or NCBI remap
"""
import logging
import operator
import os
from collections import defaultdict
from functools import reduce
from typing import Dict, List, Tuple, Any, Union

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models.query_utils import Q

from library.django_utils.django_file_utils import get_import_processing_dir
from library.guardian_utils import admin_bot
from library.log_utils import log_traceback
from library.vcf_utils import write_vcf_from_tuples
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models.models_enums import ImportSource, AlleleConversionTool, AlleleOrigin
from snpdb.models.models_genome import GenomeBuild, Contig, GenomeFasta
from snpdb.models.models_variant import AlleleSource, Liftover, Allele, Variant, VariantAlleleCollectionSource, \
    VariantAllele, VariantAlleleCollectionRecord
from upload.models import UploadedFile, UploadedLiftover, UploadPipeline, UploadedFileTypes
from upload.upload_processing import process_upload_pipeline


def create_liftover_pipelines(user: User, allele_source: AlleleSource,
                              import_source: ImportSource,
                              inserted_genome_build: GenomeBuild,
                              destination_genome_builds: List[GenomeBuild] = None):
    """ Creates and runs a liftover pipeline for each destination GenomeBuild (default = all other builds) """

    build_liftover_vcf_tuples = _get_build_liftover_tuples(allele_source, inserted_genome_build, destination_genome_builds)
    for genome_build, liftover_tuples in build_liftover_vcf_tuples.items():
        for conversion_tool, av_tuples in liftover_tuples.items():
            if conversion_tool == AlleleConversionTool.SAME_CONTIG:
                _liftover_using_same_contig(genome_build, av_tuples)
            else:
                liftover = Liftover.objects.create(user=user,
                                                   allele_source=allele_source,
                                                   conversion_tool=conversion_tool,
                                                   genome_build=genome_build)

                # Because we need to normalise / insert etc, it's easier just to write a VCF
                # and run through upload pipeline
                working_dir = get_import_processing_dir(liftover.pk, "liftover")
                liftover_vcf_filename = os.path.join(working_dir, f"liftover_variants.{genome_build.name}.vcf")
                if AlleleConversionTool.vcf_tuples_in_destination_build(conversion_tool):
                    vcf_filename = liftover_vcf_filename  # Can write directly
                else:
                    vcf_filename = os.path.join(working_dir, f"source_variants.{inserted_genome_build.name}.vcf")
                    liftover.source_vcf = vcf_filename
                    liftover.source_genome_build = inserted_genome_build
                    liftover.save()

                write_vcf_from_tuples(vcf_filename, av_tuples, tuples_have_id_field=True)
                uploaded_file = UploadedFile.objects.create(path=liftover_vcf_filename,
                                                            import_source=import_source,
                                                            name='Liftover',
                                                            user=user,
                                                            file_type=UploadedFileTypes.LIFTOVER,
                                                            visible=False)

                UploadedLiftover.objects.create(uploaded_file=uploaded_file,
                                                liftover=liftover)
                upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
                process_upload_pipeline(upload_pipeline)


VCF_ROW = Tuple[str, int, int, str, str]
LIFTOVER_TUPLE = List[Union[Tuple[int,int], VCF_ROW]]


def _get_build_liftover_tuples(allele_source: AlleleSource, inserted_genome_build: GenomeBuild,
                               destination_genome_builds: List[GenomeBuild] = None) -> Dict[Any, Dict[Any, LIFTOVER_TUPLE]]:
    """ ID column set to allele_id """
    if destination_genome_builds is None:
        destination_genome_builds = GenomeBuild.builds_with_annotation()

    other_build_contigs_q_list = []
    other_builds = set()
    for genome_build in destination_genome_builds:
        if genome_build != inserted_genome_build:
            other_builds.add(genome_build)
            q = Q(variantallele__variant__locus__contig__in=genome_build.contigs)
            other_build_contigs_q_list.append(q)

    if not other_builds:
        return {}  # Nothing to do

    other_build_contigs_q = reduce(operator.or_, other_build_contigs_q_list)

    allele_qs = allele_source.get_allele_qs().select_related("clingen_allele")
    allele_builds = defaultdict(set)
    # We want to do a left outer join to variant allele
    qs = Allele.objects.filter(other_build_contigs_q, pk__in=allele_qs)
    for allele_id, genome_build_name in qs.values_list("pk", "variantallele__genome_build"):
        allele_builds[allele_id].add(genome_build_name)

    build_liftover_vcf_tuples = defaultdict(lambda: defaultdict(list))

    for allele in allele_qs:
        existing_builds = allele_builds[allele.pk]
        for genome_build in other_builds:
            if genome_build.pk in existing_builds:
                logging.info("%s already lifted over to %s", allele, genome_build)
                continue

            conversion_tool = None
            variant_id_or_coordinate = None
            try:
                conversion_tool, variant_id_or_coordinate = allele.get_liftover_tuple(genome_build)
            except (Contig.ContigNotInBuildError, GenomeFasta.ContigNotInFastaError):
                log_traceback()

            if variant_id_or_coordinate:
                # Converted ok - return VCF tuples in desired genome build
                if conversion_tool == AlleleConversionTool.SAME_CONTIG:
                    avt = (allele.pk, variant_id_or_coordinate)
                else:
                    chrom, position, ref, alt = variant_id_or_coordinate
                    avt = (chrom, position, allele.pk, ref, alt)
                build_liftover_vcf_tuples[genome_build][conversion_tool].append(avt)
            elif settings.LIFTOVER_NCBI_REMAP_ENABLED:
                if allele.liftovererror_set.filter(liftover__genome_build=genome_build,
                                                   liftover__conversion_tool=AlleleConversionTool.NCBI_REMAP).exists():
                    continue  # Skip as already failed NCBI liftover to desired build

                # Return VCF tuples in inserted genome build
                chrom, position, ref, alt = allele.variant_for_build(inserted_genome_build).as_tuple()
                if alt == Variant.REFERENCE_ALT:
                    alt = "."  # NCBI works with '.' but not repeating ref (ie ref = alt)
                avt = (chrom, position, allele.pk, ref, alt)
                build_liftover_vcf_tuples[genome_build][AlleleConversionTool.NCBI_REMAP].append(avt)

    return build_liftover_vcf_tuples


def liftover_alleles(allele_qs, user: User = None):
    """ Creates then runs (async) liftover pipelines for a queryset of alleles """
    if user is None:
        user = admin_bot()

    for genome_build in GenomeBuild.builds_with_annotation():
        allele_source = VariantAlleleCollectionSource.objects.create(genome_build=genome_build)

        records = []
        for variant_allele in VariantAllele.objects.filter(genome_build=genome_build, allele__in=allele_qs):
            records.append(VariantAlleleCollectionRecord(collection=allele_source,
                                                         variant_allele=variant_allele))
        if records:
            VariantAlleleCollectionRecord.objects.bulk_create(records)

            variants_qs = allele_source.get_variants_qs()
            populate_clingen_alleles_for_variants(genome_build, variants_qs)
            create_liftover_pipelines(user, allele_source, ImportSource.WEB,
                                      inserted_genome_build=genome_build)


def _liftover_using_same_contig(genome_build, av_tuples: List[Tuple[int, int]]):
    """ Special case of eg Mitochondria that has the same contig across multiple builds
        we just need to create a VariantAllele object - will already have annotation for both builds """

    print(f"_liftover_using_same_contig")
    print(av_tuples)

    variant_alleles = []
    for allele_id, variant_id in av_tuples:
        va = VariantAllele(variant_id=variant_id,
                           genome_build=genome_build,
                           allele_id=allele_id,
                           origin=AlleleOrigin.LIFTOVER,
                           conversion_tool=AlleleConversionTool.SAME_CONTIG)
        variant_alleles.append(va)

    VariantAllele.objects.bulk_create(variant_alleles, ignore_conflicts=True, batch_size=2000)
