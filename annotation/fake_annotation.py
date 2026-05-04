"""
Used for testing purposes
"""
import copy
import os
from uuid import uuid4

from django.conf import settings
from django.db.models.fields import IntegerField, TextField
from django.utils import timezone

from annotation.models import ClinVarReviewStatus, GeneAnnotationRelease
from annotation.models.models import VariantAnnotationVersion, ClinVarVersion, \
    HumanProteinAtlasAnnotationVersion, AnnotationVersion, ClinVar, ClinVarCitation, \
    ClinVarCitationsCollection, VariantAnnotation, AnnotationRun, AnnotationRangeLock, GeneAnnotationVersion
from annotation.models.models_citations import CitationIdNormalized, CitationSource
from genes.hgvs import HGVSMatcher
from genes.models import GeneAnnotationImport
from genes.models_enums import AnnotationConsortium
from ontology.tests.test_data_ontology import create_ontology_test_data, create_test_ontology_version
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf


def get_fake_annotation_settings_dict(columns_version: int) -> dict:
    TEST_IMPORT_PROCESSING_DIR = os.path.join(settings.PRIVATE_DATA_ROOT, 'import_processing',
                                              "test", str(uuid4()))

    TEST_ANNOTATION = copy.deepcopy(settings.ANNOTATION)
    # phastCons/phyloP custom tracks: v1-v3 fixtures were generated without the bigwig data,
    # so disable to keep importer bindings clean. v4 fixtures include them.
    if columns_version < 4:
        TEST_ANNOTATION[settings.BUILD_GRCH37]["vep_config"].update({
            "phastcons100way": None,
            "phastcons46way": None,
            "phylop100way": None,
            "phylop46way": None,
        })
        TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"].update({
            "phastcons100way": None,
            "phastcons30way": None,
            "phylop100way": None,
            "phylop30way": None,
        })

    # columns_version 4 fixtures were generated against gnomAD 4.1 (the FILTER column shifts
    # VEPConfig.gnomad4_minor_version, which gates which gnomAD CSQ fields apply); earlier
    # fixtures used gnomAD 4.0. Plugin/version-only data files (denovo_db, mave, etc.) are
    # gated through vep_columns column-defs (min_columns_version), so we don't need to null
    # them out per-version here.
    if columns_version >= 4:
        gnomad4_path = "annotation_data/GRCh38/gnomad4.1_GRCh38_contigs.vcf.gz"
    else:
        gnomad4_path = "annotation_data/GRCh38/gnomad4.0_GRCh38_combined_af.vcf.bgz"

    # Pin gnomAD so a developer's local override doesn't shift VEP CSQ fields and break fixture parsing.
    TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"]["gnomad4"] = gnomad4_path

    ANNOTATION_COLUMNS = copy.deepcopy(TEST_ANNOTATION)
    ANNOTATION_COLUMNS[settings.BUILD_GRCH37]["columns_version"] = columns_version
    ANNOTATION_COLUMNS[settings.BUILD_GRCH38]["columns_version"] = columns_version

    return {
        "IMPORT_PROCESSING_DIR": TEST_IMPORT_PROCESSING_DIR,
        "VARIANT_ZYGOSITY_GLOBAL_COLLECTION": "global",
        "ANNOTATION_VEP_FAKE_VERSION": True,
        "ANNOTATION": ANNOTATION_COLUMNS,
    }


def get_fake_vep_version(genome_build: GenomeBuild, annotation_consortium, columns_version: int):
    # We need to use a later assembly of GRCh37 as the 1st one didn't have MT in it
    if genome_build.name == "GRCh37":
        assembly = "GRCh37.p13"
    else:
        assembly = genome_build.name

    fake_version = {"id": None,
                    "genome_build": genome_build,
                    "assembly": assembly,
                    "annotation_consortium": annotation_consortium,
                    "columns_version": columns_version}
    for f in VariantAnnotationVersion._meta.fields:  # @UndefinedVariable
        if f.name in fake_version:
            continue  # already set
        if isinstance(f, IntegerField):
            value = -1
        elif isinstance(f, TextField):
            # Need proper gnomAD for get_classified_high_frequency_variants_qs
            if f.name == 'gnomad':
                if genome_build.name == 'GRCh37':
                    value = "2.1.1"
                else:
                    value = "3.1"
            elif f.name == 'dbnsfp':
                value = '4.0a'
            else:
                value = "fake"
        else:
            continue
        fake_version[f.name] = value
    return fake_version


def get_fake_annotation_version(genome_build: GenomeBuild):
    if not settings.UNIT_TEST:
        raise ValueError("Called get_fake_annotation_version while not in a test!")

    gene_annotation_import = GeneAnnotationImport.objects.get_or_create(genome_build=genome_build,
                                                                        annotation_consortium=AnnotationConsortium.ENSEMBL,
                                                                        url="fake")[0]
    gene_annotation_release = GeneAnnotationRelease.objects.get_or_create(version=42,
                                                                          genome_build=genome_build,
                                                                          annotation_consortium=AnnotationConsortium.ENSEMBL,
                                                                          defaults={
                                                                              "gene_annotation_import": gene_annotation_import,
                                                                          })[0]

    vav_kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.ENSEMBL, 2)
    vav_kwargs["gene_annotation_release"] = gene_annotation_release
    vav_kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
    variant_annotation_version = VariantAnnotationVersion.objects.create(**vav_kwargs)
    create_ontology_test_data()
    ontology_version = create_test_ontology_version()

    gene_annotation_version = GeneAnnotationVersion.objects.get_or_create(gene_annotation_release=gene_annotation_release,
                                                                          ontology_version=ontology_version,
                                                                          gnomad_import_date=timezone.now())[0]
    clinvar_version = ClinVarVersion.objects.get_or_create(filename="fake_clinvar.vcf",
                                                           sha256_hash="not_a_real_hash",
                                                           genome_build=genome_build)[0]
    human_protein_atlas_version = HumanProteinAtlasAnnotationVersion.objects.get_or_create(filename="fake_hpa",
                                                                                           sha256_hash="not_a_real_hash",
                                                                                           hpa_version=0.42)[0]

    av, _ = AnnotationVersion.objects.get_or_create(genome_build=genome_build,
                                                    variant_annotation_version=variant_annotation_version,
                                                    gene_annotation_version=gene_annotation_version,
                                                    clinvar_version=clinvar_version,
                                                    human_protein_atlas_version=human_protein_atlas_version,
                                                    ontology_version=ontology_version)
    return av


def create_fake_variants(genome_build: GenomeBuild):
    build_lc = genome_build.name.lower()
    vcf_filename = os.path.join(settings.BASE_DIR, f"annotation/tests/test_data/test_columns_version1_{build_lc}.vep_annotated.vcf")
    slowly_create_loci_and_variants_for_vcf(genome_build, vcf_filename, get_variant_id_from_info=True)


def create_fake_clinvar_data(clinvar_version: ClinVarVersion):
    create_fake_variants(clinvar_version.genome_build)
    variant = Variant.objects.filter(Variant.get_no_reference_q()).first()

    clinvar_variation_id = 42
    clinvar_allele_id = 42

    defaults = {
        "clinvar_variation_id": clinvar_variation_id,
        "clinvar_allele_id": clinvar_allele_id,
        "preferred_disease_name": "smelly feet",
        "disease_database_name": "blah",
        "review_status": ClinVarReviewStatus.CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS,
        "clinical_significance": "Pathogenic",
        "highest_pathogenicity": 5
    }
    ClinVar.objects.get_or_create(version=clinvar_version, variant=variant, defaults=defaults)
    citation = CitationIdNormalized.from_parts(source=CitationSource.PUBMED, index=20613862).for_bulk_create()
    citation.save()
    cvcc, _ = ClinVarCitationsCollection.objects.get_or_create(pk=1)

    ClinVarCitation.objects.get_or_create(clinvar_citations_collection=cvcc,
                                          clinvar_variation_id=clinvar_variation_id,
                                          clinvar_allele_id=clinvar_allele_id,
                                          citation=citation)


def create_fake_variant_annotation(variant, variant_annotation_version: VariantAnnotationVersion) -> VariantAnnotation:
    matcher = HGVSMatcher(variant_annotation_version.genome_build)
    defaults = {
        "hgvs_g": matcher.variant_to_g_hgvs(variant)
        # ??
    }
    annotation_range_lock, _ = AnnotationRangeLock.objects.get_or_create(version=variant_annotation_version,
                                                                         min_variant=variant,
                                                                         max_variant=variant,
                                                                         count=1)
    annotation_run, _ = AnnotationRun.objects.get_or_create(annotation_range_lock=annotation_range_lock)
    va, _ = VariantAnnotation.objects.get_or_create(variant=variant, version=variant_annotation_version,
                                                    annotation_run=annotation_run, defaults=defaults)
    return va
