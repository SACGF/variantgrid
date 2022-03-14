"""
Used for testing purposes
"""
import os

from django.conf import settings
from django.db.models.fields import IntegerField, TextField
from django.utils import timezone

from annotation.models import ClinVarReviewStatus, CitationSource, GeneAnnotationRelease
from annotation.models.models import VariantAnnotationVersion, ClinVarVersion, \
    HumanProteinAtlasAnnotationVersion, AnnotationVersion, ClinVar, Citation, ClinVarCitation, \
    ClinVarCitationsCollection, VariantAnnotation, AnnotationRun, AnnotationRangeLock, GeneAnnotationVersion
from genes.models import GeneAnnotationImport
from genes.models_enums import AnnotationConsortium
from ontology.models import OntologyImport
from ontology.tests.test_data_ontology import create_ontology_test_data
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf


def get_fake_vep_version(genome_build: GenomeBuild, annotation_consortium):
    fake_version = {"genome_build": genome_build,
                    "annotation_consortium": annotation_consortium}
    for f in VariantAnnotationVersion._meta.fields:  # @UndefinedVariable
        if f.name == "id":
            continue  # Don't set PK
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
    gene_annotation_import = GeneAnnotationImport.objects.get_or_create(genome_build=genome_build,
                                                                        annotation_consortium=AnnotationConsortium.ENSEMBL,
                                                                        url="fake")[0]
    gene_annotation_release = GeneAnnotationRelease.objects.get_or_create(version=42,
                                                                          genome_build=genome_build,
                                                                          annotation_consortium=AnnotationConsortium.ENSEMBL,
                                                                          defaults={
                                                                              "gene_annotation_import": gene_annotation_import,
                                                                          })[0]

    vav_kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.ENSEMBL)
    vav_kwargs["gene_annotation_release"] = gene_annotation_release
    variant_annotation_version = VariantAnnotationVersion.objects.create(**vav_kwargs)
    create_ontology_test_data()
    last_ontology_import = OntologyImport.objects.last()

    gene_annotation_version = GeneAnnotationVersion.objects.get_or_create(gene_annotation_release=gene_annotation_release,
                                                                          last_ontology_import=last_ontology_import,
                                                                          gnomad_import_date=timezone.now())[0]
    clinvar_version = ClinVarVersion.objects.get_or_create(filename="fake_clinvar.vcf",
                                                           md5_hash="not_a_real_hash",
                                                           genome_build=genome_build)[0]
    human_protein_atlas_version = HumanProteinAtlasAnnotationVersion.objects.get_or_create(filename="fake_hpa",
                                                                                           md5_hash="not_a_real_hash",
                                                                                           hpa_version=0.42)[0]
    av, _ = AnnotationVersion.objects.get_or_create(genome_build=genome_build,
                                                    variant_annotation_version=variant_annotation_version,
                                                    gene_annotation_version=gene_annotation_version,
                                                    clinvar_version=clinvar_version,
                                                    human_protein_atlas_version=human_protein_atlas_version)
    return av


def create_fake_variants(genome_build: GenomeBuild):
    build_lc = genome_build.name.lower()
    vcf_filename = os.path.join(settings.BASE_DIR, f"annotation/tests/test_data/test_{build_lc}.vep_annotated.vcf")
    slowly_create_loci_and_variants_for_vcf(genome_build, vcf_filename, get_variant_id_from_info=True)


def create_fake_clinvar_data(clinvar_version: ClinVarVersion):
    create_fake_variants(clinvar_version.genome_build)
    variant = Variant.objects.filter(Variant.get_no_reference_q()).first()

    clinvar_variation_id = 42
    clinvar_allele_id = 42

    defaults = {
        "clinvar_variation_id": clinvar_variation_id,
        "clinvar_allele_id": clinvar_allele_id,
        "clinvar_preferred_disease_name": "smelly feet",
        "clinvar_disease_database_name": "blah",
        "clinvar_review_status": ClinVarReviewStatus.CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS,
        "clinical_significance": "Pathogenic",
        "highest_pathogenicity": 5
    }
    ClinVar.objects.get_or_create(version=clinvar_version, variant=variant, defaults=defaults)
    citation, _ = Citation.objects.get_or_create(citation_source=CitationSource.PUBMED, citation_id=20613862)
    cvcc, _ = ClinVarCitationsCollection.objects.get_or_create(pk=1)

    ClinVarCitation.objects.get_or_create(clinvar_citations_collection=cvcc,
                                          clinvar_variation_id=clinvar_variation_id,
                                          clinvar_allele_id=clinvar_allele_id,
                                          citation=citation)


def create_fake_variant_annotation(variant, variant_annotation_version: VariantAnnotationVersion):
    defaults = {
        # ??
    }
    annotation_range_lock, _ = AnnotationRangeLock.objects.get_or_create(version=variant_annotation_version,
                                                                         min_variant=variant,
                                                                         max_variant=variant,
                                                                         count=1)
    annotation_run, _ = AnnotationRun.objects.get_or_create(annotation_range_lock=annotation_range_lock)
    vav, _ = VariantAnnotation.objects.get_or_create(variant=variant, version=variant_annotation_version,
                                                     annotation_run=annotation_run, defaults=defaults)
    return vav
