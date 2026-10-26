"""
Linking what upload imports to the sequencing records it was sent for. A file names its record by the path the
SeqAuto API registered (get_sequencing_record_for_path) or, for a run-level file, by its run
(get_sequencing_run_to_link); either way an upload from the web or API needs seqauto write access
(check_may_link_to_sequencing_data). An imported VCF becomes a BackendVCF (create_backend_vcf_links), and
link_samples_and_vcfs_to_sequencing copies what seqauto knows onto the VCF and its samples - run, enrichment kit,
extraction, BAMs, QC gene lists. upload asks for this through its signals (seqauto/signals/upload_signal_handlers.py),
and a re-sent sample sheet relinks through it (seqauto/sequencing_files/sample_sheet.py).
"""
import logging
import os

from django.conf import settings
from django.utils import timezone

from library.utils import get_single_element
from patients.enums import MatchStatus
from seqauto.models.models_seqauto import (
    DragenTSO500CombinedVariantOutput,
    JointCalledVCF,
    QCGeneList,
    SampleFromSequencingSample,
    SeqAutoRecord,
    SequencingRun,
    SingleSampleVCF,
    VCFFromSequencingRun,
)
from seqauto.models.models_upload import BackendVCF
from seqauto.seqauto_permissions import can_write_seqauto
from seqauto.signals.signals_list import backend_vcf_import_start_signal
from snpdb.enums import ImportSource, ImportStatus, VariantOrigin
from upload.models.models import FileUpload, SimpleVCFImportInfo, UploadedVCF, UploadUserError

MIXED_ENRICHMENT_KIT_MESSAGE = "Joint called VCF spans sequencing runs with different enrichment kits - " \
                               "enrichment kit settings were not applied to this VCF"


def check_may_link_to_sequencing_data(file_upload: FileUpload, sequencing_data: str):
    """ Linking an upload to a seqauto record writes onto sequencing data, so a file sent through the upload page
        or API needs the seqauto API's write access - command line and SeqAuto imports are run by the server.
        Raises UploadUserError naming sequencing_data where the upload may not link to it """
    if file_upload.import_source == ImportSource.WEB_UPLOAD and not can_write_seqauto(file_upload.user):
        raise UploadUserError(f"You don't have permission to link uploads to {sequencing_data}",
                              summary="No permission to link to sequencing data")


def get_sequencing_record_for_path(file_upload: FileUpload, *seqauto_record_classes: type[SeqAutoRecord]) -> SeqAutoRecord:
    """ The first record of seqauto_record_classes (in order) for the upload's path - the sequencing file the
        upload is linked to and writes onto. Raises UploadUserError where there is no such record or the upload
        may not link to it """
    for seqauto_record_class in seqauto_record_classes:
        if record := seqauto_record_class.objects.filter(path=file_upload.path).first():
            if record.sequencing_run_id:
                sequencing_data = f"sequencing run '{record.sequencing_run_id}'"
            else:
                sequencing_data = str(record)
            check_may_link_to_sequencing_data(file_upload, sequencing_data)
            return record
    raise UploadUserError(f"No sequencing data has the path '{file_upload.path}'",
                          summary="No sequencing data at path")


def get_sequencing_run_to_link(file_upload: FileUpload, sequencing_run_name: str) -> SequencingRun:
    """ The registered run a run-level file (TSO500 MetricsOutput / CombinedVariantOutput) names in its
        metadata - such a file has no path to link it by. Raises UploadUserError where the upload may not
        link to it or there is no such run """
    check_may_link_to_sequencing_data(file_upload, f"sequencing run '{sequencing_run_name}'")
    if sequencing_run := SequencingRun.objects.filter(name=sequencing_run_name).first():
        return sequencing_run
    raise UploadUserError(f"No sequencing run named '{sequencing_run_name}'",
                          summary="No such sequencing run")


def link_uploaded_vcf_to_sequencing(uploaded_vcf: UploadedVCF, upload_step=None):
    """ Also run when a pipeline is reloaded with its VCF already created, as the link is where a
        first attempt fails when the samples don't match the sheet """
    backend_vcf = create_backend_vcf_links(uploaded_vcf)
    if backend_vcf:
        logging.info("Handle backend VCF")
        link_samples_and_vcfs_to_sequencing(backend_vcf, upload_step=upload_step)
        backend_vcf_import_start_signal.send(sender=os.path.basename(__file__), backend_vcf=backend_vcf)


def create_backend_vcf_links(uploaded_vcf):
    """ returns backend_vcf if we can link to SeqAuto data (None if Web VCF) """

    backend_vcf = None
    if not settings.SEQAUTO_ENABLED:
        # Only SeqAuto deployments link uploads to backend sequencing VCFs by path.
        # Elsewhere 'path' is just an optional client hint (e.g. API dedup), so skip linking.
        return backend_vcf

    file_upload = uploaded_vcf.file_upload
    # APIFileUploadView (including the SeqAuto API client) comes through as WEB_UPLOAD
    sequencing_vcf_sources = {ImportSource.SEQAUTO, ImportSource.WEB_UPLOAD}
    if file_upload.path and file_upload.import_source in sequencing_vcf_sources:
        sequencing_vcf = get_sequencing_record_for_path(file_upload, JointCalledVCF, SingleSampleVCF)
        joint_called_vcf = sequencing_vcf if isinstance(sequencing_vcf, JointCalledVCF) else None
        single_sample_vcf = sequencing_vcf if isinstance(sequencing_vcf, SingleSampleVCF) else None
        defaults = {"uploaded_vcf": uploaded_vcf}
        backend_vcf, _ = BackendVCF.objects.update_or_create(single_sample_vcf=single_sample_vcf,
                                                             joint_called_vcf=joint_called_vcf,
                                                             defaults=defaults)
    return backend_vcf


def link_samples_and_vcfs_to_sequencing(backend_vcf, replace_existing=False, upload_step=None):
    if backend_vcf:
        logging.info("link_samples_and_vcfs_to_sequencing backend_vcf: %s", backend_vcf.pk)
        filesystem_vcf = backend_vcf.filesystem_vcf
        sequencing_run = backend_vcf.sample_sheet.sequencing_run
        vcf = backend_vcf.vcf
        if ek := filesystem_vcf.enrichment_kit:
            vcf.variant_zygosity_count = ek.variant_zygosity_count
            logging.info("Setting %s variant_zygosity_count to %s", vcf, vcf.variant_zygosity_count)
        elif filesystem_vcf.has_mixed_enrichment_kits:
            logging.warning("%s: %s", vcf, MIXED_ENRICHMENT_KIT_MESSAGE)
            if upload_step:
                SimpleVCFImportInfo.objects.create(upload_step=upload_step,
                                                   message_string=MIXED_ENRICHMENT_KIT_MESSAGE)
        vcf.fake_data = sequencing_run.fake_data
        vcf.save()

        samples_by_sequencing_sample = backend_vcf.get_samples_by_sequencing_sample()
        if not samples_by_sequencing_sample:
            expected = ", ".join(str(ss) for ss in filesystem_vcf.get_sequencing_samples())
            raise ValueError("Couldn't link VCF samples to sequencing samples. "
                             f"VCF samples must start with one of: {expected}")

        if backend_vcf.joint_called_vcf:
            VCFFromSequencingRun.objects.update_or_create(
                vcf=vcf,
                defaults={"sequencing_run": sequencing_run,
                          "variant_caller": backend_vcf.variant_caller},
            )
        elif backend_vcf.single_sample_vcf:
            # A re-analysis replaces the arm's earlier VCF from the same caller. A BAM can have a VCF per
            # caller (eg DRAGEN TSO 500's small variants, gene-level CNV and fusions), and those stay
            sequencing_sample = get_single_element(samples_by_sequencing_sample)
            existing = VCFFromSequencingRun.objects.filter(
                vcf__sample__in=sequencing_sample.samplefromsequencingsample_set.values_list("sample"),
                vcf__uploadedvcf__backendvcf__single_sample_vcf__variant_caller=backend_vcf.variant_caller,
            )
            if existing.exists():
                existing_vcfs = ", ".join([str(vfsr.vcf) for vfsr in existing])
                logging.info("Removing existing VCFs linked against run: %s / sequencing_sample=%s",
                             existing_vcfs, sequencing_sample)
                existing.delete()

            VCFFromSequencingRun.objects.create(sequencing_run=sequencing_run,
                                                variant_caller=backend_vcf.variant_caller,
                                                vcf=vcf)

        for sequencing_sample, sample in samples_by_sequencing_sample.items():
            modified_sample = False
            # Set sample.variants_type from EnrichmentKit (eg to Mixed somatic for somatic panels)
            if ek := sequencing_sample.enrichment_kit:
                if ek.sample_variants_type != VariantOrigin.UNKNOWN:
                    sample.variants_type = ek.sample_variants_type
                    modified_sample = True

            # One link call per sequencing sample reaches all of that arm's VCFs through here
            if sequencing_sample.extraction_id and not sample.extraction_id:
                sample.extraction = sequencing_sample.extraction
                sample.extraction_match_status = MatchStatus.MATCHED
                sample.extraction_match_date = timezone.now()
                modified_sample = True
            elif sequencing_sample.extraction_reference and not sample.extraction_reference:
                # Still parked upstream - carry the claim so one reconcile pass settles both rows
                sample.extraction_reference = sequencing_sample.extraction_reference
                sample.extraction_match_status = sequencing_sample.extraction_match_status
                sample.extraction_match_date = sequencing_sample.extraction_match_date
                modified_sample = True

            for alignment_file in sequencing_sample.alignmentfile_set.all():
                alignment_file.create_sample_file_path(sample)

            if modified_sample:
                sample.save()

            try:
                sfss = SampleFromSequencingSample.objects.get(sample=sample)
                if replace_existing or sfss.sample.import_status in ImportStatus.DELETION_STATES:
                    logging.info("Updating SampleFromSequencingSample")
                    sfss.sequencing_sample = sequencing_sample
                    sfss.sample = sample
                    sfss.save()
                else:
                    logging.warning("SS %s already linked to sample: %s (%s/%s)", sfss.sequencing_sample, sfss.sample,
                                    sfss.sample.pk, sfss.sample.import_status)
            except SampleFromSequencingSample.DoesNotExist:
                SampleFromSequencingSample.objects.create(sample=sample,
                                                          sequencing_sample=sequencing_sample)

            DragenTSO500CombinedVariantOutput.link_arm_sample(sample, sequencing_run)

            # Link any QCGeneLists
            for qcgl in QCGeneList.objects.filter(qc__sequencing_sample=sequencing_sample,
                                                  custom_text_gene_list__gene_list__isnull=False,
                                                  sample_gene_list__isnull=True).distinct():
                qcgl.create_and_assign_sample_gene_list(sample)
