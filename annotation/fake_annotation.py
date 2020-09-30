"""
Used for testing purposes
"""
import os

from django.conf import settings
from django.db.models.fields import IntegerField, TextField

from annotation.models import ClinVarReviewStatus, CitationSource
from annotation.models.models import VariantAnnotationVersion, EnsemblGeneAnnotationVersion, ClinVarVersion, \
    HumanProteinAtlasAnnotationVersion, AnnotationVersion, ClinVar, Citation, ClinVarCitation, \
    ClinVarCitationsCollection, EnsemblGeneAnnotation, VariantAnnotation, AnnotationRun, AnnotationRangeLock
from genes.models_enums import AnnotationConsortium
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
            value = "fake"
        else:
            continue
        fake_version[f.name] = value
    return fake_version


def get_fake_annotation_version(genome_build: GenomeBuild):
    vav_kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.ENSEMBL)
    variant_annotation_version = VariantAnnotationVersion.objects.create(**vav_kwargs)
    ensembl_gene_annotation_version = EnsemblGeneAnnotationVersion.objects.get_or_create(filename="blah",
                                                                                         md5_hash="not_a_real_hash",
                                                                                         genome_build=genome_build,
                                                                                         ensembl_version=42)[0]
    clinvar_version = ClinVarVersion.objects.get_or_create(filename="fake_clinvar.vcf",
                                                           md5_hash="not_a_real_hash",
                                                           genome_build=genome_build)[0]
    human_protein_atlas_version = HumanProteinAtlasAnnotationVersion.objects.get_or_create(filename="fake_hpa",
                                                                                           md5_hash="not_a_real_hash",
                                                                                           hpa_version=0.42)[0]
    av, _ = AnnotationVersion.objects.get_or_create(genome_build=genome_build,
                                                    variant_annotation_version=variant_annotation_version,
                                                    ensembl_gene_annotation_version=ensembl_gene_annotation_version,
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


def create_fake_gene_annotation(ensembl_gene_annotation_version: EnsemblGeneAnnotationVersion):
    defaults = {
        'hgnc_symbol': 'RUNX1',
        'external_gene_name': 'RUNX1',
        'hgnc_symbol_lower': 'runx1',
        'hgnc_name': 'runt related transcription factor 1',
        'synonyms': 'PEBP2A2, AMLCR1',
        'previous_symbols': 'AML1, CBFA2',
        'hgnc_chromosome': '21q22.12',
        'gene_family_tag': None,
        'gene_family_description': None,
        'hgnc_id': 'HGNC:10471',
        'entrez_gene_id': 861,
        'uniprot_id': 'Q01196',
        'ucsc_id': 'uc002yuk.6',
        'omim_id': '151385',
        'enzyme_ids': None,
        'ccds_ids': 'CCDS13639, CCDS42922, CCDS46646',
        'rgd_id': 'RGD:2283',
        'mgi_id': 'MGI:99852',
        'rvis_percentile': '37.1137060628',
        'refseq_gene_summary': 'Core binding factor (CBF) is a heterodimeric transcription factor that binds to the core element of many enhancers and promoters. The protein encoded by this gene represents the alpha subunit of CBF and is thought to be involved in the development of normal hematopoiesis. Chromosomal translocations involving this gene are well-documented and have been associated with several types of leukemia. Three transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jul 2008]',
        'function_from_uniprotkb': "CBF binds to the core site, 5'-PYGPYGGT-3', of a number of enhancers and promoters, including murine leukemia virus, polyomavirus enhancer, T-cell receptor enhancers, LCK, IL-3 and GM-CSF promoters. The alpha subunit binds DNA and appears to have a role in the development of normal hematopoiesis. Isoform AML-1L interferes with the transactivation activity of RUNX1. Acts synergistically with ELF4 to transactivate the IL-3 promoter and with ELF2 to transactivate the mouse BLK promoter. Inhibits KAT6B- dependent transcriptional activation. Controls the anergy and suppressive function of regulatory T-cells (Treg) by associating with FOXP3. Activates the expression of IL2 and IFNG and down- regulates the expression of TNFRSF18, IL2RA and CTLA4, in conventional T-cells (PubMed:17377532). Positively regulates the expression of RORC in T-helper 17 cells (By similarity). {ECO:0000250|UniProtKB:Q03347, ECO:0000269|PubMed:10207087, ECO:0000269|PubMed:11965546, ECO:0000269|PubMed:14970218, ECO:0000269|PubMed:17377532, ECO:0000269|PubMed:17431401}.",
        'pathway_from_uniprotkb': None,
        'tissue_specificity_from_uniprotkb': 'Expressed in all tissues examined except brain and heart. Highest levels in thymus, bone marrow and peripheral blood.',
        'phenotypes_from_ensembl': 'Acute myeloid leukemia with t(8;21)(q22;q22) translocation :: Chronic myeloid leukemia :: Hereditary thrombocytopenia with normal platelets-hematological cancer predisposition syndrome :: Isolated delta-storage pool disease :: Precursor B-cell acute lymphoblastic leukemia',
        'omim_phenotypes': 'Platelet disorder, familial, with associated myeloid malignancy :: Leukemia, acute myeloid',
        'gene_biotype': 'protein_coding',
        'status': 'K',
        'chromosome_name': '21',
        'start_position': 36160098,
        'end_position': 37376965,
        'band': 'q22.12',
        'strand': '-',
        'percentage_gc_content': 41.54,
        'transcript_count': 19,
        'in_cancer_gene_census': True
    }
    ega, _ = EnsemblGeneAnnotation.objects.get_or_create(gene_id='ENSG00000159216',
                                                         version=ensembl_gene_annotation_version,
                                                         defaults=defaults)

    return ega
