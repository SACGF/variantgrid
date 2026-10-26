"""
    Model + signal operations on SampleSheet that are shared by the REST API and web UI.

    These were extracted from the old filesystem-scanning create_resource_models.py - see
    https://github.com/SACGF/variantgrid/issues/1643 - keeping only what the API-populated data needs.

    Entry points: current_sample_sheet_changed (a different sheet became current - move the run's data onto it)
    and current_sample_sheet_data_changed (a re-post changed the rows of the sheet already current). Both send
    sequencing_run_current_sample_sheet_changed_signal, which SA Pathology uses to match patients, and record
    what they did with api_request_note (#2349).
"""
import logging
import os
from collections import Counter

from eventlog.models import api_request_note
from library.django_utils.change_tracking import active_collector
from seqauto.models import IlluminaFlowcellQC, JointCalledVCF
from seqauto.signals.signals_list import (
    current_sample_sheet_data_moved_signal,
    sequencing_run_current_sample_sheet_changed_signal,
)


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
        return True, "Barcodes changed: " + ', '.join(diff_barcodes)

    return False, None


def _move_to_sequencing_sample(queryset, current_sequencing_sample) -> int:
    """ QuerySet.update sends no post_save, so the API log's change collector is told directly """
    if num_moved := queryset.update(sequencing_sample=current_sequencing_sample):
        if collector := active_collector():
            collector.record(queryset.model, updated=num_moved, fields=["sequencing_sample"])
    return num_moved


def _reassign_sequencing_data_to_current_samples(old_samples_by_name, current_samples_by_name) -> Counter:
    """ Reassign the existing Fastq / UnalignedReads / AlignmentFile / SingleSampleVCF / QC records, and the
        links to VCF samples, onto the current sheet's SequencingSample with the same sample name.
        Returns the number moved, by what they are """
    moved = Counter()
    for sequencing_sample_name, old_sequencing_sample in old_samples_by_name.items():
        current_ss = current_samples_by_name.get(sequencing_sample_name)
        if current_ss is None:
            continue

        moved["fastqs"] += _move_to_sequencing_sample(old_sequencing_sample.fastq_set.all(), current_ss)
        moved["unaligned reads"] += _move_to_sequencing_sample(old_sequencing_sample.unalignedreads_set.all(), current_ss)
        moved["alignment files"] += _move_to_sequencing_sample(old_sequencing_sample.alignmentfile_set.all(), current_ss)
        moved["single sample VCFs"] += _move_to_sequencing_sample(old_sequencing_sample.singlesamplevcf_set.all(), current_ss)
        moved["QC"] += _move_to_sequencing_sample(old_sequencing_sample.qc_set.all(), current_ss)
        moved["sample links"] += _move_to_sequencing_sample(old_sequencing_sample.samplefromsequencingsample_set.all(), current_ss)
    return moved


def _move_joint_called_vcf_members(joint_called_vcf, old_sample_sheet, new_sample_sheet) -> bool:
    """ Point a cross-run joint call's explicit members at the new sheet's matching rows (match by sample_id).
        Returns whether any members moved """
    old_members = list(joint_called_vcf.sequencing_samples.filter(sample_sheet=old_sample_sheet))
    if not old_members:
        return False

    moved = False
    new_by_sample_id = {ss.sample_id: ss for ss in new_sample_sheet.sequencingsample_set.all()}
    for old_member in old_members:
        if new_member := new_by_sample_id.get(old_member.sample_id):
            joint_called_vcf.sequencing_samples.remove(old_member)
            joint_called_vcf.sequencing_samples.add(new_member)
            logging.info("Moved %s member onto current sample sheet: %s", joint_called_vcf, new_member)
            moved = True
        else:
            logging.info("%s member %s has no match on %s", joint_called_vcf, old_member, new_sample_sheet)
    return moved


def current_sample_sheet_changed(sequencing_run_current_sample_sheet, new_sample_sheet):
    """ Points the run at new_sample_sheet, moves the data stored against its older sheets onto it, then sends
        sequencing_run_current_sample_sheet_changed_signal - after the move, so a receiver matching patients
        finds the new sheet's samples linked to theirs """
    old_sample_sheet = sequencing_run_current_sample_sheet.sample_sheet
    sequencing_run_current_sample_sheet.sample_sheet = new_sample_sheet
    sequencing_run_current_sample_sheet.save()
    sequencing_run = sequencing_run_current_sample_sheet.sequencing_run

    meaningfully_changed, reason = sample_sheet_meaningfully_changed(old_sample_sheet, new_sample_sheet)
    api_request_note(f"Current sample sheet changed {old_sample_sheet.hash} -> {new_sample_sheet.hash}: "
                     f"{reason or 'not meaningfully changed'}")
    new_sample_names = set(new_sample_sheet.get_sequencing_samples_by_name())
    if unmatched := sorted(set(old_sample_sheet.get_sequencing_samples_by_name()) - new_sample_names):
        api_request_note(f"Not on the new sheet, their data stays on {old_sample_sheet.hash}: {', '.join(unmatched)}")

    assign_old_sample_sheet_data_to_current_sample_sheet(sequencing_run)
    sequencing_run_current_sample_sheet_changed_signal.send(sender=os.path.basename(__file__),
                                                            sequencing_run=sequencing_run,
                                                            meaningfully_changed=meaningfully_changed)


def current_sample_sheet_data_changed(sequencing_run):
    """ A re-post changed the rows of the sheet that is already current - eg a SAPOrderNumber column that
        arrived late - so receivers re-match patients as they would for a new sheet """
    api_request_note("Sample sheet data changed, patient matching re-run")
    sequencing_run_current_sample_sheet_changed_signal.send(sender=os.path.basename(__file__),
                                                            sequencing_run=sequencing_run,
                                                            meaningfully_changed=False)


def _plural_counts(counts: Counter) -> str:
    return ", ".join(f"{count} {name}" for name, count in counts.items() if count)


def assign_old_sample_sheet_data_to_current_sample_sheet(sequencing_run):
    """ Reassign the data the API already stored against older SampleSheets onto the current one.

        Operates purely on model data (the API sends records for the files; there is no filesystem
        to scan): relink old IlluminaFlowcellQC / JointCalledVCF, reassign the existing
        Fastq / UnalignedReads / AlignmentFile / VCF / QC records by sample name, then send
        current_sample_sheet_data_moved_signal so the joint calls' samples are relinked. """
    current_sample_sheet = sequencing_run.get_current_sample_sheet()
    old_sample_sheets = sequencing_run.get_old_sample_sheets()
    moved = Counter()

    if not IlluminaFlowcellQC.objects.filter(sample_sheet=current_sample_sheet).exists():
        illumina_qc = IlluminaFlowcellQC.objects.filter(sample_sheet__in=old_sample_sheets).order_by("-pk").first()
        if illumina_qc:
            illumina_qc.sample_sheet = current_sample_sheet
            illumina_qc.save()
            moved["IlluminaFlowcellQC"] += 1

    for joint_called_vcf in JointCalledVCF.objects.filter(sample_sheet__in=old_sample_sheets):
        joint_called_vcf.sample_sheet = current_sample_sheet
        joint_called_vcf.save()
        moved["joint called VCFs"] += 1

    # Cross-run joint calls are owned by another run, so move their members over by sample_id
    member_joint_called_vcfs = JointCalledVCF.objects.filter(
        sequencing_samples__sample_sheet__in=old_sample_sheets).distinct()
    for joint_called_vcf in member_joint_called_vcfs:
        for old_sample_sheet in old_sample_sheets:
            if _move_joint_called_vcf_members(joint_called_vcf, old_sample_sheet, current_sample_sheet):
                moved["joint call memberships"] += 1

    current_samples_by_name = current_sample_sheet.get_sequencing_samples_by_name()
    for old_sample_sheet in old_sample_sheets:
        moved += _reassign_sequencing_data_to_current_samples(old_sample_sheet.get_sequencing_samples_by_name(),
                                                              current_samples_by_name)
    if moved:
        api_request_note(f"Moved onto sample sheet {current_sample_sheet.hash}: {_plural_counts(moved)}")

    relink_samples = bool(moved["joint called VCFs"] or moved["joint call memberships"])
    if not relink_samples:
        missing_linked_sequencing_samples = current_sample_sheet.sequencingsample_set.filter(samplefromsequencingsample__isnull=True)
        relink_samples = missing_linked_sequencing_samples.exists()
    # upload relinks the joint called VCFs it imported (upload/signals/signal_handlers.py)
    current_sample_sheet_data_moved_signal.send(sender=os.path.basename(__file__),
                                                sample_sheet=current_sample_sheet,
                                                relink_samples=relink_samples)

    logging.info("Assigned old data to current sample sheet %s", current_sample_sheet)
    sequencing_run.save()  # Re-validate "ready"
