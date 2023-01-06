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
from classification.autopopulate_evidence_keys.evidence_from_sample_and_patient import \
    get_evidence_fields_for_sample_and_patient
from classification.autopopulate_evidence_keys.evidence_from_variant import get_evidence_fields_for_variant, \
    AutopopulateData
from classification.enums import SubmissionSource, SpecialEKeys
from classification.models import EvidenceKey, Classification, ClassificationImport
from classification.tasks.classification_import_process_variants_task import liftover_classification_import
from library.git import Git
from snpdb.models import GenomeBuild, ImportSource, Sample, UserSettings, Variant


def create_classification_for_sample_and_variant_objects(
        user: User,
        sample: Sample,
        variant: Variant,
        genome_build: GenomeBuild,
        refseq_transcript_accession: str = None,
        ensembl_transcript_accession: str = None,
        annotation_version: str = None):
    """ Create internally from existing variant - not used by API which may need to create variants """
    user_settings = UserSettings.get_for_user(user)
    # TODO - if you have > 1 labs then redirect to pick page.
    lab = user_settings.get_lab()
    vc_import = ClassificationImport.objects.create(user=user, genome_build=genome_build)
    kwargs = {"user": user,
              "lab": lab,
              "variant": variant,
              "sample": sample,
              "classification_import": vc_import,
              "populate_with_defaults": True}

    classification = Classification.create(**kwargs)
    classification_auto_populate_fields(classification, genome_build,
                                        refseq_transcript_accession=refseq_transcript_accession,
                                        ensembl_transcript_accession=ensembl_transcript_accession,
                                        annotation_version=annotation_version)
    #classification.set_variant(variant)  # have to re-do this because we didn't have the transcript the 1st time around
    classification.ensure_allele_info()

    liftover_classification_import(vc_import, ImportSource.WEB)
    return classification


def generate_auto_populate_data(
        variant: Variant,
        genome_build: Optional[GenomeBuild] = None,
        annotation_version: Optional[AnnotationVersion] = None,
        refseq_transcript_accession: Optional[str] = None,
        ensembl_transcript_accession: Optional[str] = None,
        sample: Optional[Sample] = None) -> AutopopulateData:
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
                                                evidence_keys_list, annotation_version))

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
        save: bool = True):
    """
    Applies annotation data to the classification
    """

    auto_data = generate_auto_populate_data(
        variant=classification.variant,
        genome_build=genome_build,
        refseq_transcript_accession=refseq_transcript_accession,
        ensembl_transcript_accession=ensembl_transcript_accession,
        sample=classification.sample,
        annotation_version=annotation_version
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
