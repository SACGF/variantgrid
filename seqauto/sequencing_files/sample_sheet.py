"""
    Model + signal operations on SampleSheet that are shared by the REST API and web UI.

    These were extracted from the old filesystem-scanning create_resource_models.py - see
    https://github.com/SACGF/variantgrid/issues/1643 - keeping only what the API-populated data needs.
"""
import logging
import os

from library.log_utils import log_traceback
from seqauto.models import JointCalledVCF, IlluminaFlowcellQC, SampleFromSequencingSample, \
    get_samples_by_sequencing_sample
from seqauto.signals.signals_list import sequencing_run_current_sample_sheet_changed_signal
from upload.models import BackendVCF
from upload.vcf.vcf_import import link_samples_and_vcfs_to_sequencing


def sample_sheet_meaningfully_changed(sample_sheet_1, sample_sheet_2):
    """ returns True, "string reason" if downstream stuff has to change """

    ss_dict_1 = sample_sheet_1.get_sequencing_samples_by_name()
    ss_dict_2 = sample_sheet_2.get_sequencing_samples_by_name()

    ss_1_set = set(ss_dict_1)
    ss_2_set = set(ss_dict_2)

    if ss_1_set != ss_2_set:
        reason = "SampleSheet SampleIDs don't match, samples not in common: %s" % str(ss_1_set ^ ss_2_set)
        return True, reason

    diff_barcodes = []
    for sequencing_sample_name, ss1 in ss_dict_1.items():
        ss2 = ss_dict_2[sequencing_sample_name]

        if ss1.barcode != ss2.barcode:
            diff_barcodes.append(f"{sequencing_sample_name} ({ss1.barcode}/{ss2.barcode})")

    if diff_barcodes:
        return True, ','.join(diff_barcodes)

    return False, None


def _reassign_sequencing_data_to_current_samples(old_samples_by_name, current_samples_by_name):
    """ Reassign the existing Fastq / UnalignedReads records onto the current sheet's SequencingSample
        with the same sample name (match by name). Used when a newer SampleSheet replaces an old one. """
    for sequencing_sample_name, old_sequencing_sample in old_samples_by_name.items():
        current_ss = current_samples_by_name.get(sequencing_sample_name)
        if current_ss is None:
            continue

        for fastq in old_sequencing_sample.fastq_set.all():
            logging.info("Updating FastQ with current sequencing sample")
            fastq.sequencing_sample = current_ss
            fastq.save()

        for unaligned_reads in old_sequencing_sample.unalignedreads_set.all():
            logging.info("Updating Unaligned reads with current sequencing sample")
            unaligned_reads.sequencing_sample = current_ss
            unaligned_reads.save()


def current_sample_sheet_changed(sequencing_run_current_sample_sheet, new_sample_sheet):
    logging.info("current_sample_sheet_changed - SS changed: %s", new_sample_sheet.path)

    old_sample_sheet = sequencing_run_current_sample_sheet.sample_sheet
    sequencing_run_current_sample_sheet.sample_sheet = new_sample_sheet
    sequencing_run_current_sample_sheet.save()

    sequencing_run = sequencing_run_current_sample_sheet.sequencing_run

    old_version = old_sample_sheet.get_version_string()
    new_version = new_sample_sheet.get_version_string()
    message = f"Updated CurrentSampleSheet: {old_sample_sheet.path}. Was: {old_version}, now: {new_version}"
    meaningfully_changed, reason = sample_sheet_meaningfully_changed(old_sample_sheet, new_sample_sheet)
    if meaningfully_changed:
        logging.warning("%s SampleSheets CHANGED: %s", message, reason)
    else:
        logging.info("%s SampleSheets not meaningfully changed", message)

    if not meaningfully_changed:
        # Can go through and update models etc...
        # Update samples from sequencing sample
        for joint_called_vcf in JointCalledVCF.objects.filter(sequencing_run=sequencing_run):
            if vcf := joint_called_vcf.vcf:
                # TODO: Raise some kind of manual task here to fix things???

                def _get_samples_by_sequencing_sample_id(sample_sheet):
                    return {ss.sample_id: s for ss, s in get_samples_by_sequencing_sample(sample_sheet, vcf).items()}

                old_samples_by_sequencing_sample_id = _get_samples_by_sequencing_sample_id(old_sample_sheet)
                new_samples_by_sequencing_sample_id = _get_samples_by_sequencing_sample_id(new_sample_sheet)

                for sequencing_sample_id, sample in old_samples_by_sequencing_sample_id.items():
                    new_sample = new_samples_by_sequencing_sample_id[sequencing_sample_id]
                    try:
                        new_ss = new_sample.samplefromsequencingsample.sequencing_sample
                    except SampleFromSequencingSample.DoesNotExist:
                        # Does this matter? maybe vcf wasn't imported etc...
                        logging.info("There was no SampleFromSequencingSample for sample %s", new_sample)
                        continue

                    try:
                        sfss = sample.samplefromsequencingsample
                        sfss.sequencing_sample = new_ss
                        sfss.save()
                        logging.info("Updated SampleFromSequencingSample with %s", new_ss)
                    except SampleFromSequencingSample.DoesNotExist:
                        # Does this matter? maybe vcf wasn't imported etc...
                        logging.info("There was no SampleFromSequencingSample for sample %s", sample)

        # Update downstream models - reassign Fastq / UnalignedReads to the current sequencing samples
        _reassign_sequencing_data_to_current_samples(old_sample_sheet.get_sequencing_samples_by_name(),
                                                     new_sample_sheet.get_sequencing_samples_by_name())

    sequencing_run_current_sample_sheet_changed_signal.send(sender=os.path.basename(__file__),
                                                            sequencing_run=sequencing_run,
                                                            meaningfully_changed=meaningfully_changed)

    sequencing_run.save()  # Re-validate ready based on errors


def assign_old_sample_sheet_data_to_current_sample_sheet(user, sequencing_run):
    """ Reassign the data the API already stored against older SampleSheets onto the current one.

        Operates purely on model data (the API sends records for the files; there is no filesystem
        to scan): relink old IlluminaFlowcellQC / JointCalledVCF and reassign the existing
        Fastq / UnalignedReads records by sample name. """
    current_sample_sheet = sequencing_run.get_current_sample_sheet()
    old_sample_sheets = sequencing_run.get_old_sample_sheets()

    # Relink old IlluminaFlowcellQC to the current sample sheet
    try:
        current_sample_sheet.illuminaflowcellqc
    except Exception:  # not linked to current
        illumina_qc = IlluminaFlowcellQC.objects.filter(sample_sheet__in=old_sample_sheets).order_by("-pk").first()
        if illumina_qc:
            logging.info("Assigning IlluminaQC to latest sample sheet")
            illumina_qc.sample_sheet = current_sample_sheet
            illumina_qc.save()

    # Relink old JointCalledVCFs to the current sample sheet
    old_joint_called_vcfs = JointCalledVCF.objects.filter(sample_sheet__in=old_sample_sheets)
    relink_samples = False
    for joint_called_vcf in old_joint_called_vcfs:
        try:
            logging.info("Assigning JointCalledVCF %s to latest sample sheet", joint_called_vcf)
            joint_called_vcf.sample_sheet = current_sample_sheet
            joint_called_vcf.save()
            relink_samples = True
        except Exception:
            log_traceback()

    if not relink_samples:
        missing_linked_sequencing_samples = current_sample_sheet.sequencingsample_set.filter(samplefromsequencingsample__isnull=True)
        relink_samples = missing_linked_sequencing_samples.exists()

    try:
        for joint_called_vcf in current_sample_sheet.jointcalledvcf_set.all():
            backend_vcf = joint_called_vcf.backendvcf

            if relink_samples or joint_called_vcf.needs_to_be_linked():
                link_samples_and_vcfs_to_sequencing(backend_vcf, replace_existing=True)
    except BackendVCF.DoesNotExist:
        pass

    # Reassign the existing Fastq / UnalignedReads records from the old sample sheets' sequencing
    # samples onto the current sheet's matching sequencing samples (match by sample name).
    current_samples_by_name = current_sample_sheet.get_sequencing_samples_by_name()
    for old_sample_sheet in old_sample_sheets:
        _reassign_sequencing_data_to_current_samples(old_sample_sheet.get_sequencing_samples_by_name(),
                                                     current_samples_by_name)

    logging.info("Assigned old data to current sample sheet %s", current_sample_sheet)

    sequencing_run.save()  # Re-validate "ready"
