"""
    Auto-populate a Classification object with evidence keys from
    info we can get from VariantGrid DB
"""
import socket
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.sites.models import Site

from annotation.models.models import AnnotationVersion
from classification.autopopulate_evidence_keys.evidence_from_sample_and_patient import (
    get_evidence_fields_for_sample_and_patient,
)
from classification.autopopulate_evidence_keys.evidence_from_variant import (
    AutopopulateData,
    get_evidence_fields_for_variant,
)
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationImport,
    ClassificationModification,
    EvidenceKey,
)
from classification.models.classification import COPY_SCOPES_ALL, COPY_SCOPES_GENE
from classification.models.classification_variant_fields_validation import (
    apply_somatic_tier_from_amp_level,
)
from classification.tasks.classification_import_process_variants_task import (
    liftover_classification_import,
)
from library.git import Git
from snpdb.clingen_allele_api import ClinGenAlleleRegistryAPI
from snpdb.models import GenomeBuild, ImportSource, Lab, Sample, Variant


def create_classification_for_sample_and_variant_objects(
        user: User,
        lab: Lab,
        sample: Optional[Sample],
        variant: Variant,
        genome_build: GenomeBuild,
        refseq_transcript_accession: str = None,
        ensembl_transcript_accession: str = None,
        annotation_version: str = None,
        clingen_api: ClinGenAlleleRegistryAPI = None):
    """ Create internally from existing variant - not used by API which may need to create variants """

    kwargs = {"user": user,
              "lab": lab,
              "variant": variant,
              "sample": sample,
              "populate_with_defaults": True}

    classification = Classification.create(**kwargs)
    classification_populate_from_variant(classification, genome_build,
                                         refseq_transcript_accession=refseq_transcript_accession,
                                         ensembl_transcript_accession=ensembl_transcript_accession,
                                         annotation_version=annotation_version,
                                         clingen_api=clingen_api)
    return classification


def classification_populate_from_variant(
        classification: Classification,
        genome_build: GenomeBuild,
        refseq_transcript_accession: str = None,
        ensembl_transcript_accession: str = None,
        annotation_version: str = None,
        clingen_api: ClinGenAlleleRegistryAPI = None):
    """ What a new record gets from its variant - evidence from annotation, its allele info and liftover """
    classification_auto_populate_fields(classification, genome_build,
                                        refseq_transcript_accession=refseq_transcript_accession,
                                        ensembl_transcript_accession=ensembl_transcript_accession,
                                        annotation_version=annotation_version,
                                        clingen_api=clingen_api)

    allele_info, allele_info_created = classification.ensure_allele_info_with_created()
    if allele_info and allele_info_created:
        vc_import = ClassificationImport.objects.create(user=classification.user, genome_build=genome_build)
        allele_info.set_variant_and_save(matched_variant=classification.variant)
        allele_info.classification_import = vc_import
        allele_info.save()
        liftover_classification_import(vc_import, ImportSource.WEB)

    # if the allele_info has already linked to an allele
    # call this to make sure the allele gets set
    classification.apply_allele_info_to_classification()


def classification_complete_web_create(
        classification: Classification,
        user: User,
        evidence: Optional[dict] = None,
        copy_from: Optional[ClassificationModification] = None,
        copy_gene_from: Optional[ClassificationModification] = None):
    """ After populating from the variant: the evidence the form sent, publish, then any record the curator
        chose to copy from """
    if evidence:
        classification.patch_value(
            patch=evidence,
            clear_all_fields=False,
            user=user,
            source=SubmissionSource.FORM,
            leave_existing_values=True,
            save=True,
            make_patch_fields_immutable=False)

    classification.publish_latest(user)

    # Allele level first so it beats the gene level copy, which only ever fills what is still empty
    for source, copy_scopes in [(copy_from, COPY_SCOPES_ALL), (copy_gene_from, COPY_SCOPES_GENE)]:
        if source:
            ClassificationConsensus(modification=source, copy_scopes=copy_scopes).apply_to(classification, user)

    apply_somatic_tier_from_amp_level(classification, user)


def generate_auto_populate_data(
        variant: Variant,
        genome_build: Optional[GenomeBuild] = None,
        annotation_version: Optional[AnnotationVersion] = None,
        refseq_transcript_accession: Optional[str] = None,
        ensembl_transcript_accession: Optional[str] = None,
        sample: Optional[Sample] = None,
        clingen_api: Optional[ClinGenAlleleRegistryAPI] = None) -> AutopopulateData:
    """
    Shows all complete auto-populate data for a would be classification.
    """

    if annotation_version is None:
        annotation_version = AnnotationVersion.latest(genome_build)

    evidence_keys_list = list(EvidenceKey.objects.all().select_related("variantgrid_column"))

    data = AutopopulateData("basic")
    data[SpecialEKeys.GENOME_BUILD] = genome_build.get_build_with_patch(annotation_version)
    data[SpecialEKeys.CURATION_SYSTEM] = get_curation_system()

    # Used to be a check if variant existed, but pretty sure it can be guaranteed to exist at this point
    data.annotation_version = annotation_version
    data.update(get_evidence_fields_for_variant(genome_build, variant,
                                                refseq_transcript_accession, ensembl_transcript_accession,
                                                evidence_keys_list, annotation_version, clingen_api=clingen_api))

    if sample:
        data.update(get_evidence_fields_for_sample_and_patient(variant, sample))

    data[SpecialEKeys.AUTOPOPULATE] = {"value": data.summary, "immutable": SubmissionSource.VARIANT_GRID}
    return data


def classification_auto_populate_fields(
        classification: Classification,
        genome_build: GenomeBuild,
        refseq_transcript_accession: str = None,
        ensembl_transcript_accession: str = None,
        leave_existing_values: bool = True,
        annotation_version: str = None,
        save: bool = True,
        clingen_api: ClinGenAlleleRegistryAPI = None):
    """
    Applies annotation data to the classification
    """

    auto_data = generate_auto_populate_data(
        variant=classification.variant,
        genome_build=genome_build,
        refseq_transcript_accession=refseq_transcript_accession,
        ensembl_transcript_accession=ensembl_transcript_accession,
        sample=classification.sample,
        annotation_version=annotation_version,
        clingen_api=clingen_api
    )
    classification.annotation_version = auto_data.annotation_version
    return classification.patch_value(auto_data.data,
                                      user=classification.user,
                                      source=SubmissionSource.VARIANT_GRID,
                                      leave_existing_values=leave_existing_values,
                                      save=save)


def get_curation_system():
    curation_system = "VariantGrid"

    details = {
        "site": Site.objects.get_current(),
        "hostname": socket.gethostname(),
        "git": Git(settings.BASE_DIR).hash,
    }
    explain = ""
    site_name = getattr(settings, "SITE_NAME", None)
    if site_name and site_name != curation_system:
        explain = f"{site_name}: "
    explain += ", ".join([f"{k}={v}" for k, v in details.items()])
    return {"value": curation_system,
            "explain": explain}
