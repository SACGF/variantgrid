import logging
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from typing import Iterable, Optional

import pandas as pd
from django.conf import settings
from django.db.models.query_utils import Q

from genes.canonical_transcripts.canonical_transcript_manager import CanonicalTranscriptManager
from genes.gene_matching import GeneSymbolMatcher
from genes.models import TranscriptVersion
from library.enums.log_level import LogLevel
from library.file_utils import name_from_filename, file_md5sum, file_to_array
from library.log_utils import get_traceback, log_traceback
from seqauto.illumina.run_parameters import get_run_parameters
from seqauto.illumina.samplesheet import convert_sheet_to_df, samplesheet_is_valid
from seqauto.models import Sequencer, SequencingRun, SequencingSample, SequencingSampleData, Fastq, SampleSheet, \
    UnalignedReads, BamFile, VCFFile, QC, SampleSheetCombinedVCFFile, IlluminaFlowcellQC, FastQC, Flagstats, \
    DontAutoLoadException, Experiment, SequencingRunCurrentSampleSheet, SampleFromSequencingSample, QCGeneList, \
    get_samples_by_sequencing_sample, QCGeneCoverage, SeqAutoMessage, SeqAutoRecord, get_variant_caller_from_vcf_file
from seqauto.models.models_enums import SequencingFileType
from seqauto.signals import sequencing_run_current_sample_sheet_changed_signal, sequencing_run_created_signal, \
    sequencing_run_sample_sheet_created_signal
from snpdb.models import DataState
from upload.models import BackendVCF
from upload.vcf.vcf_import import link_samples_and_vcfs_to_sequencing


class SeqAutoRunError(Exception):
    pass


def stripped_lines_set(lines):
    stripped_lines = set()
    for line in lines:
        stripped_lines.add(line.strip())
    return stripped_lines


def returning_existing_records_by_qs(qs, filter_kwargs=None):
    existing_records = {}
    if filter_kwargs:
        qs = qs.filter(**filter_kwargs)
    for record in qs:
        existing_records[record.path] = record
    return existing_records


def returning_existing_records_by_path(klass, filter_kwargs=None):
    qs = klass.objects.all()
    return returning_existing_records_by_qs(qs, filter_kwargs=filter_kwargs)


def create_samplesheet_samples(sample_sheet):
    try:
        date_on_file = None
        # Handle having prefixes added to SequencingRun dir - get off original path
        if original_sequencing_run := sample_sheet.sequencing_run.get_params().get("original_sequencing_run"):
            filename_parts = original_sequencing_run.split('_')
            date_on_file = filename_parts[0]

        df = convert_sheet_to_df(sample_sheet.path, date_on_file=date_on_file)
        sample_sheet_process_default, msg = samplesheet_is_valid(df)
    except Exception as e:
        logging.error("samplesheet_is_valid() died for %s", sample_sheet.path)
        exc_traceback = sys.exc_info()[2]
        logging.error(''.join(traceback.format_tb(exc_traceback)))
        raise e

    if not sample_sheet_process_default:
        logging.info("Not automatically processing samplesheet %s due to: %s", sample_sheet, msg)

    control_sample_pattern = None
    if settings.SEQAUTO_CONTROL_SAMPLE_REGEX:
        control_sample_pattern = re.compile(settings.SEQAUTO_CONTROL_SAMPLE_REGEX)

    sequencing_samples = []
    sample_number = 1
    for _, row in df.iterrows():
        sample_project = row["sample_project"]
        if pd.isnull(sample_project):
            sample_project = None

        sample_name = row["sample_name"]
        is_control = False
        if control_sample_pattern and control_sample_pattern.findall(sample_name):
            logging.info("%s = CONTROL", sample_name)
            is_control = True

        failed = row["failed"]
        automatically_process = sample_sheet_process_default and not (is_control or failed)

        ss = SequencingSample.objects.create(sample_sheet=sample_sheet,
                                             sample_number=sample_number,
                                             sample_id=row["sample_id"],
                                             sample_name=sample_name,
                                             sample_project=sample_project,
                                             lane=row["lane"],
                                             barcode=row["barcode"],
                                             is_control=is_control,
                                             failed=failed,
                                             automatically_process=automatically_process)
        sequencing_samples.append(ss)
        sample_number += 1

        # Extra columns
        for column in settings.SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS:
            if column in row:
                value = row[column]
                SequencingSampleData.objects.create(sequencing_sample=ss,
                                                    column=column,
                                                    value=value)

    return sequencing_samples


def load_from_file_if_complete(seqauto_run, record, **kwargs):
    if record.data_state == DataState.COMPLETE:
        try:
            logging.debug("Auto-loaded %s", record)
            record.load_from_file(seqauto_run, **kwargs)
        except DontAutoLoadException:
            logging.debug("Skipping auto-load of %s due to settings", record)
        except:
            error_exception = get_traceback()
            logging.error("Problem loading %s: '%s'", record, record.path)
            logging.error(error_exception)
            record.data_state = DataState.ERROR
            SeqAutoMessage.objects.update_or_create(seqauto_run=seqauto_run, record=record,
                                                    severity=LogLevel.ERROR, code="load_error",
                                                    defaults={"open": True, "message": error_exception})

        record.save()


class FlowcellChecker:

    def __init__(self):
        self.skip_patterns = []

        skip_patterns = getattr(settings, "SEQAUTO_SKIP_FLOWCELLS_PATTERNS")
        if skip_patterns:
            self.skip_patterns.extend(skip_patterns)

        skip_file = getattr(settings, "SEQAUTO_SKIP_FLOWCELLS_FILE")
        if skip_file:
            with open(skip_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("#"):
                        continue
                    self.skip_patterns.append(line)

        logging.info("Flowcell Skip patterns:")
        logging.info(self.skip_patterns)

        self.skip_flowcell_filename = getattr(settings, "SEQAUTO_SKIP_INDIVIDUAL_FLOWCELL_FILE")

    def skip(self, sequencing_run_dir):
        for p in self.skip_patterns:
            if re.match(p, sequencing_run_dir):
                return True
        if self.skip_flowcell_filename:
            skip_filename = os.path.join(sequencing_run_dir, self.skip_flowcell_filename)
            if os.path.exists(skip_filename):
                return True

        return False


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


def current_sample_sheet_changed(seqauto_run, sequencing_run_current_sample_sheet, new_sample_sheet):
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
        message += f"SampleSheets CHANGED: {reason}"
        SeqAutoMessage.objects.create(seqauto_run=seqauto_run, record=sequencing_run,
                                      message=message, severity=LogLevel.WARNING)
    else:
        message += "SampleSheets not meaningfully changed"
        SeqAutoMessage.objects.create(seqauto_run=seqauto_run, record=sequencing_run,
                                      message=message, severity=LogLevel.INFO)

    if not meaningfully_changed:
        # Can go through and update models etc...
        # Update samples from sequencing sample
        for combo_vcf in SampleSheetCombinedVCFFile.objects.filter(sequencing_run=sequencing_run):
            if vcf := combo_vcf.vcf:
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

        # Update downstream models
        old_samples_by_name = old_sample_sheet.get_sequencing_samples_by_name()
        new_samples_by_name = new_sample_sheet.get_sequencing_samples_by_name()

        for sequencing_sample_name, old_sequencing_sample in old_samples_by_name.items():
            new_ss = new_samples_by_name[sequencing_sample_name]

            for fastq in old_sequencing_sample.fastq_set.all():
                logging.info("Updating FastQ with current sequencing sample")
                fastq.sequencing_sample = new_ss
                fastq.save()

            for unaligned_reads in old_sequencing_sample.unalignedreads_set.all():
                logging.info("Updating Unaligned reads with current sequencing sample")
                unaligned_reads.sequencing_sample = new_ss
                unaligned_reads.save()

    sequencing_run_current_sample_sheet_changed_signal.send(sender=os.path.basename(__file__),
                                                            sequencing_run=sequencing_run,
                                                            meaningfully_changed=meaningfully_changed)

    sequencing_run.save()  # Re-validate ready based on errors


def assign_old_sample_sheet_data_to_current_sample_sheet(user, sequencing_run):
    current_sample_sheet = sequencing_run.get_current_sample_sheet()
    old_sample_sheets = sequencing_run.get_old_sample_sheets()

    try:
        current_sample_sheet.illuminaflowcellqc
    except:  # Exception - not linked to current
        try:
            illumina_qc = IlluminaFlowcellQC.objects.filter(sample_sheet__in=old_sample_sheets).order_by("-pk").first()
            logging.info("Assigning IlluminaQC to latest sample sheet")
            illumina_qc.sample_sheet = current_sample_sheet
            illumina_qc.save()
        except:
            pass

    old_sample_sheet_combined_vcf = SampleSheetCombinedVCFFile.objects.filter(sample_sheet__in=old_sample_sheets)
    relink_samples = False
    try:
        combo_vcf = old_sample_sheet_combined_vcf.get()
        logging.info("Assigning SampleSheetCombinedVCFFile to latest sample sheet")
        combo_vcf.sample_sheet = current_sample_sheet
        combo_vcf.save()

        relink_samples = True
    except:
        log_traceback()

    if not relink_samples:
        missing_linked_sequencing_samples = current_sample_sheet.sequencingsample_set.filter(samplefromsequencingsample__isnull=True)
        relink_samples = missing_linked_sequencing_samples.exists()

    try:
        for combo_vcf in current_sample_sheet.samplesheetcombinedvcffile_set.all():
            backend_vcf = combo_vcf.backendvcf

            if relink_samples or combo_vcf.needs_to_be_linked():
                link_samples_and_vcfs_to_sequencing(backend_vcf, replace_existing=True)
    except BackendVCF.DoesNotExist:
        pass

    # SeqAuto data
    old_unaligned_reads_qs = UnalignedReads.get_for_old_sample_sheets(sequencing_run)
    q_r1_or_r2 = Q(fastq_r1__in=old_unaligned_reads_qs) | Q(fastq_r2__in=old_unaligned_reads_qs)
    fastq_qs = Fastq.objects.filter(q_r1_or_r2)
    existing_fastq_records = returning_existing_records_by_qs(fastq_qs)

    existing_fastq_files = set()
    for fastq in existing_fastq_records.values():
        if os.path.exists(fastq.path):
            existing_fastq_files.add(fastq.path)

    fastq_pairs = defaultdict(dict)
    create_fastqs_for_sample_sheet(current_sample_sheet, existing_fastq_files, existing_fastq_records, fastq_pairs,
                                   assign_to_most_recent_sample_sheet=True)
    unaligned_reads_list_from_fastqs(fastq_pairs)

    # Delete anything left...
    UnalignedReads.get_for_old_sample_sheets(sequencing_run).delete()

    message = f"Assigning old data to current sample sheet {current_sample_sheet.hash} ({current_sample_sheet.last_modified_datetime})"
    SeqAutoMessage.objects.create(record=sequencing_run,
                                  severity=LogLevel.INFO,
                                  code="assign_old_date_to_current_sample_sheet",
                                  message=message)

    sequencing_run.save()  # Re-validate "ready"


def process_flowcells(seqauto_run, existing_files, results):
    existing_sequencing_runs = stripped_lines_set(existing_files)
    existing_sequencing_run_records = returning_existing_records_by_path(SequencingRun)
    sequencers = Sequencer.get_sequencers_dict_by_name()
    flowcell_checker = FlowcellChecker()

    sequencing_runs = {}
    # TODO: Only handling existing - need to handle deleted flowcells + set to deleted
    for sequencing_run_dir in existing_sequencing_runs:
        sequencing_run = existing_sequencing_run_records.get(sequencing_run_dir)
        if sequencing_run and sequencing_run.legacy:
            continue
        try:
            sequencing_run = process_sequencing_run(seqauto_run, sequencers, flowcell_checker, sequencing_run_dir, sequencing_run)
            if sequencing_run:
                # Process all files regardless of state
                sequencing_runs[sequencing_run.pk] = sequencing_run
        except Exception as e:
            raise SeqAutoRunError(f"Error processing {sequencing_run_dir}") from e

    return sequencing_runs


def process_sequencing_run(seqauto_run, sequencers, flowcell_checker, sequencing_run_dir, sequencing_run) -> Optional[SequencingRun]:
    # TODO: Handle skipped for this and spreadsheet missing below together
    if flowcell_checker.skip(sequencing_run_dir):
        message = f"Skipping Flowcell path '{sequencing_run_dir}'"
        logging.info(message)
        # sequencing_run may have been created pre-skipping
        if sequencing_run and sequencing_run.data_state != DataState.SKIPPED:
            sequencing_run.data_state = DataState.SKIPPED
            sequencing_run.save()

            SeqAutoMessage.objects.create(seqauto_run=seqauto_run, record=sequencing_run,
                                          message=message, severity=LogLevel.INFO)

        return None

    data_state = DataState.COMPLETE

    run_parameters_dir = os.path.join(sequencing_run_dir, settings.SEQAUTO_RUN_PARAMETERS_SUB_DIR)
    instrument_name, experiment_name = get_run_parameters(run_parameters_dir)
    if experiment_name:
        experiment, _ = Experiment.objects.get_or_create(name=experiment_name)
    else:
        experiment = None

    if sequencing_run is None:
        try:
            sequencer = Sequencer.get_or_create_sequencer(instrument_name, sequencers)
        except Exception as e:
            msg = f"Could not process sequencer from run parameters in {sequencing_run_dir}"
            raise ValueError(msg) from e

        name = os.path.basename(sequencing_run_dir)
        date = SequencingRun.get_date_from_name(name)
        sequencing_run = SequencingRun.objects.create(name=name,
                                                      date=date,
                                                      sequencer=sequencer,
                                                      path=sequencing_run_dir,
                                                      experiment=experiment,
                                                      data_state=data_state,
                                                      fake_data=seqauto_run.fake_data)
        sequencing_run.sequencing_run = sequencing_run  # Wat
        sequencing_run.save()

        sequencing_run_created_signal.send(sender=os.path.basename(__file__),
                                           sequencing_run=sequencing_run)
    else:
        if experiment and sequencing_run.experiment is None:
            logging.info("Setting Experiment on existing SequencingRun")
            sequencing_run.experiment = experiment

    sequencing_run.has_basecalls = sequencing_run.check_basecalls_dir()
    sequencing_run.has_interop = sequencing_run.check_interop_dir()

    # Set to skipped if no sample sheet
    samplesheet_path = SampleSheet.get_path_from_sequencing_run(sequencing_run)
    if not os.path.exists(samplesheet_path):
        message = f"Skipping Flowcell '{sequencing_run_dir}' as missing '{samplesheet_path}'"
        SeqAutoMessage.objects.create(seqauto_run=seqauto_run, record=sequencing_run,
                                      message=message, severity=LogLevel.ERROR)
        data_state = DataState.SKIPPED

    sequencing_run.data_state = data_state
    sequencing_run.save()

    if sequencing_run.data_state != DataState.SKIPPED:
        sample_sheet, created = SampleSheet.objects.get_or_create(sequencing_run=sequencing_run,
                                                                  hash=file_md5sum(samplesheet_path),
                                                                  path=samplesheet_path)
        if created:
            sample_sheet.file_last_modified = SeqAutoRecord.get_file_last_modified(samplesheet_path)
            sample_sheet.save()
            create_samplesheet_samples(sample_sheet)
            sequencing_run_sample_sheet_created_signal.send(sender=os.path.basename(__file__),
                                                            sample_sheet=sample_sheet)
        # Make sure SequencingRunCurrentSampleSheet is set to what we found on disk
        try:
            # Update existing
            current_ss = sequencing_run.sequencingruncurrentsamplesheet
            on_disk_not_current = current_ss.sample_sheet != sample_sheet
            if on_disk_not_current:
                current_sample_sheet_changed(seqauto_run, current_ss, sample_sheet)

        except SequencingRunCurrentSampleSheet.DoesNotExist:
            # Create new
            current_ss = SequencingRunCurrentSampleSheet.objects.create(sequencing_run=sequencing_run,
                                                                        sample_sheet=sample_sheet)
            logging.info("Created new SequencingRunCurrentSampleSheet: %s", current_ss)

    return sequencing_run


def get_data_state(input_data_state, file_exists):
    if input_data_state == DataState.SKIPPED:
        data_state = DataState.SKIPPED
    elif file_exists:
        data_state = DataState.COMPLETE
    elif input_data_state == DataState.DELETED:
        data_state = DataState.DELETED
    else:
        data_state = DataState.NON_EXISTENT
    return data_state


def process_illuminate_qc(seqauto_run, existing_files, results):
    existing_illumina_qcs = stripped_lines_set(existing_files)
    uses_current_samplesheet = {"sample_sheet__sequencingruncurrentsamplesheet__isnull": False}
    existing_illumina_qc_records = returning_existing_records_by_path(IlluminaFlowcellQC, uses_current_samplesheet)

    illumina_flowcell_qcs = []
    sequencing_runs = results[SequencingFileType.SAMPLE_SHEET]
    for sequencing_run in sequencing_runs.values():
        illuminate_path = IlluminaFlowcellQC.get_path_from_sequencing_run(sequencing_run)

        exists = illuminate_path in existing_illumina_qcs

        if sequencing_run.data_state == DataState.SKIPPED:
            data_state = DataState.SKIPPED
        elif exists:
            data_state = DataState.COMPLETE
        else:
            if sequencing_run.can_generate_qc():
                data_state = DataState.NON_EXISTENT
            elif sequencing_run.data_state == DataState.DELETED:
                data_state = DataState.DELETED
            else:
                logging.info("Can't run IlluminaQC on %s (data_state=%s)", sequencing_run.path, sequencing_run.data_state)
                data_state = DataState.ERROR

        illumina_qc = existing_illumina_qc_records.get(illuminate_path)
        if illumina_qc:
            if illumina_qc.data_state != data_state:
                logging.info("IlluminaQC %s data state changed from %s->%s", illumina_qc, illumina_qc.data_state, data_state)
                illumina_qc.data_state = data_state
                illumina_qc.save()
                load_from_file_if_complete(seqauto_run, illumina_qc)
        else:
            if DataState.should_create_new_record(data_state):
                logging.info("IlluminaQC: creating for %s", illuminate_path)
                sample_sheet = sequencing_run.get_current_sample_sheet()
                illumina_qc = IlluminaFlowcellQC.objects.create(sequencing_run=sequencing_run,
                                                                sample_sheet=sample_sheet,
                                                                path=illuminate_path,
                                                                data_state=data_state)
                load_from_file_if_complete(seqauto_run, illumina_qc)

        if illumina_qc:
            illumina_flowcell_qcs.append(illumina_qc)

    return illumina_flowcell_qcs


def create_fastqs_for_sample_sheet(sample_sheet, existing_fastq_files, existing_fastq_records, fastq_pairs, assign_to_most_recent_sample_sheet=False):
    fastqs = []
    for sequencing_sample in sample_sheet.sequencingsample_set.all():
        pair_paths = Fastq.get_pair_paths_from_sequencing_sample(sequencing_sample)
        # Try to match one of the path options.
        paths = []
        for pp in pair_paths:
            if pp[0] in existing_fastq_files or pp[1] in existing_fastq_files:
                paths = pp
                break

        data_state = get_data_state(sample_sheet.sequencing_run.data_state, paths)
        if data_state == DataState.NON_EXISTENT and not sample_sheet.sequencing_run.can_basecall():
            data_state = DataState.DELETED  # Can't do anything...

        # logging.info("paths=%s, data_state=%s", paths, data_state)

        if data_state == DataState.NON_EXISTENT:
            paths = pair_paths[0]  # 1st is what current pipeline will make

        # Why isn't R2 in this??
        for fastq_path in paths:
            # We want to be able to skip these sometimes
            if not sequencing_sample.automatically_process:
                # logging.info("Skipping sequencing sample %s", sequencing_sample)
                data_state = DataState.SKIPPED

            fastq = existing_fastq_records.get(fastq_path)
            if fastq:
                require_save = False
                if sequencing_sample != fastq.sequencing_sample:
                    if assign_to_most_recent_sample_sheet:
                        logging.info("Re-assigning to most recent sample sheet")
                        try:
                            unaligned_reads = fastq.fastq_r1.get()
                            unaligned_reads.sequencing_sample = sequencing_sample
                            unaligned_reads.save()
                        except:
                            pass

                        fastq.sequencing_sample = sequencing_sample
                    else:
                        logging.warning("FastQ %s on disk set for old samplesheet", fastq)
                        fastq.data_state = DataState.ERROR
                    require_save = True

                if fastq.data_state != data_state:
                    logging.info("FastQ %s data state changed from %s->%s", fastq, fastq.data_state, data_state)
                    fastq.data_state = data_state
                    require_save = True

                if require_save:
                    fastq.save()
            else:
                if DataState.should_create_new_record(data_state):
                    name = name_from_filename(fastq_path)
                    fastq = Fastq.objects.create(sequencing_run=sequencing_sample.sample_sheet.sequencing_run,
                                                 sequencing_sample=sequencing_sample,
                                                 path=fastq_path,
                                                 name=name,
                                                 data_state=data_state)

            if fastq and fastq.data_state != DataState.SKIPPED:
                fastqs.append(fastq)

                pair_name, read = fastq.get_common_pair_name_and_read()
                fastq_pairs[pair_name][read] = fastq

    return fastqs


def unaligned_reads_list_from_fastqs(fastq_pairs):
    unaligned_reads_list = []

    for reads in fastq_pairs.values():
        fastq_r1 = reads.get("R1")
        fastq_r2 = reads.get("R2")
        if fastq_r1 and fastq_r2:
            sequencing_sample = fastq_r1.sequencing_sample
            try:
                unaligned_reads = UnalignedReads.objects.get(sequencing_sample=sequencing_sample)
                check_and_update_unaligned_reads(unaligned_reads, fastq_r1, fastq_r2)
            except UnalignedReads.DoesNotExist:
                unaligned_reads = UnalignedReads.objects.create(sequencing_sample=sequencing_sample,
                                                                fastq_r1=fastq_r1,
                                                                fastq_r2=fastq_r2)
            unaligned_reads_list.append(unaligned_reads)
    return unaligned_reads_list


def process_fastq(seqauto_run, existing_files, results):
    existing_fastq_files = stripped_lines_set(existing_files)
    existing_fastq_records = returning_existing_records_by_path(Fastq)

    sequencing_runs = results[SequencingFileType.SAMPLE_SHEET]
    fastq_pairs = defaultdict(dict)
    fastqs = []

    for sequencing_run in sequencing_runs.values():
        try:
            sample_sheet = sequencing_run.get_current_sample_sheet()
        except:
            logging.info("process_fastq is skipping flowcell %s", sequencing_run)
            continue

        fastqs.extend(create_fastqs_for_sample_sheet(sample_sheet, existing_fastq_files, existing_fastq_records, fastq_pairs))

    print_data_state_stats(fastqs)

    return unaligned_reads_list_from_fastqs(fastq_pairs)


# Unaligned reads may have been created with NON_EXISTANT data state, but now we've found FastQs but with a different path.
def check_and_update_unaligned_reads(unaligned_reads, fastq_r1, fastq_r2):

    def ok_to_replace(old_fastq, new_fastq):
        return old_fastq.data_state in [DataState.NON_EXISTENT, DataState.DELETED] and new_fastq.data_state in [DataState.COMPLETE, DataState.NON_EXISTENT]

    def bad_fastq_replace(old_fastq, new_fastq):
        params = (unaligned_reads.sequencing_sample, old_fastq.read, old_fastq.path, old_fastq, old_fastq.data_state, new_fastq, new_fastq.path, new_fastq.data_state)
        msg = "Existing UnalignedRead for sequencing_sample %s has read %s '%s' %s (state: %s) path not equal to %s '%s' (state: %s)" % params
        raise Exception(msg)

    replaced_fastqs = []

    if unaligned_reads.fastq_r1 != fastq_r1:
        if ok_to_replace(unaligned_reads.fastq_r1, fastq_r1):
            replaced_fastqs.append(unaligned_reads.fastq_r1)
            unaligned_reads.fastq_r1 = fastq_r1
        else:
            bad_fastq_replace(unaligned_reads.fastq_r1, fastq_r1)

    if unaligned_reads.fastq_r2 != fastq_r2:
        if ok_to_replace(unaligned_reads.fastq_r2, fastq_r2):
            replaced_fastqs.append(unaligned_reads.fastq_r2)
            unaligned_reads.fastq_r2 = fastq_r2
        else:
            bad_fastq_replace(unaligned_reads.fastq_r2, fastq_r2)

    if replaced_fastqs:
        logging.info("Replacing FastQs for UnalignedReads %s and deleting: %s", unaligned_reads, ','.join(map(str, replaced_fastqs)))
        unaligned_reads.save()

        for fastq in replaced_fastqs:
            fastq.delete()


def process_fastqc(seqauto_run, existing_files, results):
    existing_fastqcs = stripped_lines_set(existing_files)
    existing_fastqc_records = returning_existing_records_by_path(FastQC)

    fastqcs = []
    unaligned_reads_list = results[SequencingFileType.FASTQ]
    for unaligned_reads in unaligned_reads_list:
        for fastq in [unaligned_reads.fastq_r1, unaligned_reads.fastq_r2]:
            fastqc_path = FastQC.get_path_from_fastq(fastq)

            exists = fastqc_path in existing_fastqcs
            data_state = get_data_state(fastq.data_state, exists)
            fastqc = existing_fastqc_records.get(fastqc_path)
            if fastqc:
                if fastqc.data_state != data_state:
                    logging.info("FastQC %s data state changed from %s->%s", fastqc, fastqc.data_state, data_state)
                    fastqc.data_state = data_state
                    fastqc.save()
            else:
                if DataState.should_create_new_record(data_state):
                    fastqc = FastQC.objects.create(sequencing_run=fastq.sequencing_run,
                                                   fastq=fastq,
                                                   path=fastqc_path,
                                                   data_state=data_state)
                    load_from_file_if_complete(seqauto_run, fastqc)

            if fastqc:
                fastqcs.append(fastqc)

    return fastqcs


def get_expected_bams_and_unaligned_reads(unaligned_reads_list):
    expected_bams_and_unaligned_reads = {}
    for unaligned_reads in unaligned_reads_list:
        bam_path = BamFile.get_path_from_unaligned_reads(unaligned_reads)
        expected_bams_and_unaligned_reads[bam_path] = unaligned_reads

    return expected_bams_and_unaligned_reads


def process_bam(seqauto_run, existing_files, results):
    existing_bam_files = stripped_lines_set(existing_files)
    existing_bam_records = returning_existing_records_by_path(BamFile)
    return process_bam_files_and_records(seqauto_run, existing_bam_files, existing_bam_records, results)


def process_bam_files_and_records(seqauto_run, existing_bam_files, existing_bam_records, results):
    unaligned_reads_list = results[SequencingFileType.FASTQ]
    expected_bams_and_unaligned_reads = get_expected_bams_and_unaligned_reads(unaligned_reads_list)
    bam_files = set()

    for bam_path, unaligned_reads in expected_bams_and_unaligned_reads.items():
        exists = bam_path in existing_bam_files
        data_state = get_data_state(unaligned_reads.data_state, exists)

        if bam := existing_bam_records.get(bam_path):
            if bam.data_state != data_state:
                logging.info("Bam %s data state changed from %s->%s", bam, bam.data_state, data_state)
                bam.data_state = data_state
                bam.save()
        else:
            if DataState.should_create_new_record(data_state):
                name = name_from_filename(bam_path)
                aligner = BamFile.get_aligner_from_bam_file(bam_path)
                sequencing_run = unaligned_reads.fastq_r1.sequencing_run
                bam = BamFile.objects.create(sequencing_run=sequencing_run,
                                             unaligned_reads=unaligned_reads,
                                             path=bam_path,
                                             name=name,
                                             data_state=data_state,
                                             aligner=aligner)
        if bam:
            bam_files.add(bam)

    return bam_files


def process_flagstats(seqauto_run, existing_files, results):
    existing_flagstats_files = stripped_lines_set(existing_files)
    existing_flagstats_records = returning_existing_records_by_path(Flagstats)

    bams = results[SequencingFileType.BAM]
    flagstats_list = []
    for bam_file in bams:
        flagstats_path = Flagstats.get_path_from_bam_file(bam_file)

        exists = flagstats_path in existing_flagstats_files
        data_state = get_data_state(bam_file.data_state, exists)

        flagstats = existing_flagstats_records.get(flagstats_path)
        if flagstats:
            if flagstats.data_state != data_state:
                logging.info("Flagstats %s data state changed from %s->%s", flagstats, flagstats.data_state, data_state)
                flagstats.data_state = data_state
                flagstats.save()
        else:
            if DataState.should_create_new_record(data_state):
                flagstats = Flagstats.objects.create(sequencing_run=bam_file.sequencing_run,
                                                     bam_file=bam_file,
                                                     path=flagstats_path,
                                                     data_state=data_state)
                load_from_file_if_complete(seqauto_run, flagstats)

        if flagstats:
            flagstats_list.append(flagstats)

    return flagstats_list


def get_expected_vcf_and_bams(bams):
    expected_vcf_and_bams = {}
    for bam in bams:
        vcf_path = VCFFile.get_path_from_bam(bam)
        expected_vcf_and_bams[vcf_path] = bam

    return expected_vcf_and_bams


def process_single_sample_vcfs(seqauto_run, existing_files, results):
    existing_vcf_files = stripped_lines_set(existing_files)
    existing_vcf_records = returning_existing_records_by_path(VCFFile)

    bams = results[SequencingFileType.BAM]
    expected_vcf_and_bams = get_expected_vcf_and_bams(bams)
    vcf_files = set()

    for vcf_path, bam_file in expected_vcf_and_bams.items():
        exists = vcf_path in existing_vcf_files
        data_state = get_data_state(bam_file.data_state, exists)

        vcf_file = existing_vcf_records.get(vcf_path)
        if vcf_file:
            if vcf_file.data_state != data_state:
                logging.info("vcf_file %s data state changed from %s->%s", vcf_file, vcf_file.data_state, data_state)
                vcf_file.data_state = data_state
                vcf_file.save()
        else:
            if DataState.should_create_new_record(data_state):
                variant_caller = get_variant_caller_from_vcf_file(vcf_path)
                vcf_file = VCFFile.objects.create(sequencing_run=bam_file.sequencing_run,
                                                  bam_file=bam_file,
                                                  path=vcf_path,
                                                  data_state=data_state,
                                                  variant_caller=variant_caller)
                load_from_file_if_complete(seqauto_run, vcf_file)

        if vcf_file:
            vcf_files.add(vcf_file)

    return vcf_files


def process_combo_vcfs(seqauto_run, existing_files, results):
    existing_combo_vcf_files = stripped_lines_set(existing_files)
    existing_combo_vcf_records = returning_existing_records_by_path(SampleSheetCombinedVCFFile)

    combined_vcf_files = set()
    sequencing_runs = results[SequencingFileType.SAMPLE_SHEET]
    for sequencing_run in sequencing_runs.values():
        try:
            sample_sheet = sequencing_run.get_current_sample_sheet()
        except:
            continue

        for vcf_path in SampleSheetCombinedVCFFile.get_paths_from_sample_sheet(sample_sheet):
            exists = vcf_path in existing_combo_vcf_files
            data_state = get_data_state(sequencing_run.data_state, exists)
            combined_vcf_file = existing_combo_vcf_records.get(vcf_path)

            if combined_vcf_file:
                require_save = False
                if sample_sheet != combined_vcf_file.sample_sheet:
                    logging.info("Combo VCF %s on disk set for old sample_sheet, setting to latest", combined_vcf_file)
                    combined_vcf_file.sample_sheet = sample_sheet
                    require_save = True

                if combined_vcf_file.data_state != data_state:
                    logging.info("combined_vcf_file %s data state changed from %s->%s", combined_vcf_file,
                                 combined_vcf_file.data_state, data_state)
                    combined_vcf_file.data_state = data_state
                    require_save = True

                if require_save:
                    combined_vcf_file.save()
            else:
                if DataState.should_create_new_record(data_state):
                    variant_caller = get_variant_caller_from_vcf_file(vcf_path)
                    combined_vcf_file = SampleSheetCombinedVCFFile.objects.create(sequencing_run=sequencing_run,
                                                                                  sample_sheet=sample_sheet,
                                                                                  path=vcf_path,
                                                                                  data_state=data_state,
                                                                                  variant_caller=variant_caller)
            if combined_vcf_file:
                combined_vcf_files.add(combined_vcf_file)
                try:
                    combined_vcf_file.backendvcf  # If there, VCF is already loaded
                except:
                    load_from_file_if_complete(seqauto_run, combined_vcf_file)
    return combined_vcf_files


def process_vcf(seqauto_run, existing_files, results):
    return {"vcf": process_single_sample_vcfs(seqauto_run, existing_files, results),
            "combined_vcf": process_combo_vcfs(seqauto_run, existing_files, results)}


def get_expected_qc_and_vcfs(vcfs):
    expected_qc_and_vcfs = {}
    for vcf in vcfs:
        qc_path = QC.get_path_from_vcf(vcf)
        expected_qc_and_vcfs[qc_path] = vcf
    return expected_qc_and_vcfs


def process_qc(seqauto_run, existing_files, results):
    existing_qcs = stripped_lines_set(existing_files)
    existing_qc_records = returning_existing_records_by_path(QC)

    data = results[SequencingFileType.VCF]
    vcfs = data['vcf']  # Single sample VCFs
    expected_qc_and_vcfs = get_expected_qc_and_vcfs(vcfs)
    qc_set = set()

    for qc_path, vcf_file in expected_qc_and_vcfs.items():
        exists = qc_path in existing_qcs
        data_state = get_data_state(vcf_file.data_state, exists)

        qc = existing_qc_records.get(qc_path)
        if qc:
            if qc.data_state == DataState.ERROR:  # Leave errored ones alone.
                continue

            if qc.data_state != data_state:
                logging.info("QC %s data state changed from %s->%s", qc, qc.data_state, data_state)
                qc.data_state = data_state
                qc.save()
                load_from_file_if_complete(seqauto_run, qc)
        else:
            if DataState.should_create_new_record(data_state):
                qc = QC.objects.create(path=qc_path,
                                       sequencing_run=vcf_file.sequencing_run,
                                       bam_file=vcf_file.bam_file,
                                       vcf_file=vcf_file,
                                       data_state=data_state)
                load_from_file_if_complete(seqauto_run, qc)

        if qc:
            qc_set.add(qc)

    process_other_qc(seqauto_run, existing_qcs, qc_set)
    return qc_set


def process_other_qc(seqauto_run, existing_qcs, qc_set):
    gene_matcher = GeneSymbolMatcher()
    canonical_transcript_manager = CanonicalTranscriptManager()
    transcript_versions_by_id = TranscriptVersion.transcript_versions_by_id()  # Don't know build - just get all

    process_other_qc_class(seqauto_run, gene_matcher, canonical_transcript_manager, transcript_versions_by_id,
                           existing_qcs, qc_set, QCGeneList)
    if settings.SEQAUTO_LOAD_GENE_COVERAGE:
        process_other_qc_class(seqauto_run, gene_matcher, canonical_transcript_manager, transcript_versions_by_id,
                               existing_qcs, qc_set, QCGeneCoverage, delete_out_of_date_records=True)


def process_other_qc_class(seqauto_run, gene_matcher, canonical_transcript_manager, transcript_versions_by_id,
                           existing_qcs, qc_set, klass, delete_out_of_date_records=False):
    logging.info("process_other_qc_class: %s", klass)

    existing_other_qc_records = returning_existing_records_by_path(klass)
    created = 0
    updated = 0
    reloaded = 0
    missing = 0
    unchanged = 0

    for qc in qc_set:
        other_qc_path = klass.get_path_from_qc(qc)
        exists = other_qc_path in existing_qcs
        data_state = get_data_state(qc.data_state, exists)

        if exists:
            file_last_modified = SeqAutoRecord.get_file_last_modified(other_qc_path)
        else:
            file_last_modified = 0.0

        def _create_new_record():
            return klass.objects.create(path=other_qc_path,
                                        sequencing_run=qc.sequencing_run,
                                        qc=qc,
                                        data_state=data_state,
                                        file_last_modified=file_last_modified)

        record = existing_other_qc_records.get(other_qc_path)
        if record:
            if record.data_state == DataState.ERROR:  # Leave errored ones alone.
                continue

            if record.data_state != data_state:
                updated += 1
                logging.info("QC %s data state changed from %s->%s", record, record.data_state, data_state)
                if data_state == DataState.COMPLETE:
                    record.file_last_modified = file_last_modified

                record.data_state = data_state
                record.save()
                load_from_file_if_complete(seqauto_run,
                                           record,
                                           gene_matcher=gene_matcher,
                                           canonical_transcript_manager=canonical_transcript_manager,
                                           transcript_versions_by_id=transcript_versions_by_id)
            else:
                # logging.info("%s - previous last modified: %f - now %f",
                #             record.path, record.file_last_modified, file_last_modified)
                # Compare last modified at second level - otherwise we get float equality issues
                if exists and int(file_last_modified) > int(record.file_last_modified):
                    logging.info("**** RECORD last modified has changed: old: %s, new: %s ***",
                                 record.file_last_modified, file_last_modified)

                    if delete_out_of_date_records:
                        # Can only have 1 QCGeneCoverage for QC
                        record.delete()

                    record = _create_new_record()
                    load_from_file_if_complete(seqauto_run,
                                               record,
                                               gene_matcher=gene_matcher,
                                               canonical_transcript_manager=canonical_transcript_manager,
                                               transcript_versions_by_id=transcript_versions_by_id)
                    reloaded += 1
                else:
                    unchanged += 1
        else:
            if DataState.should_create_new_record(data_state):
                created += 1
                record = _create_new_record()
                load_from_file_if_complete(seqauto_run,
                                           record,
                                           gene_matcher=gene_matcher,
                                           canonical_transcript_manager=canonical_transcript_manager,
                                           transcript_versions_by_id=transcript_versions_by_id)
            else:
                missing += 1

    logging.info("%s: created: %d, updated: %d, reloaded: %d, missing: %d, unchanged: %d",
                 klass, created, updated, reloaded, missing, unchanged)


def print_data_state_stats(data):
    if isinstance(data, dict):
        data = data.values()

    c = Counter()
    for d in data:
        if isinstance(d, Iterable):
            for sd in d:
                c[sd.data_state] += 1
        else:
            c[d.data_state] += 1
    logging.info(c)


def create_resource_models(seqauto_run, seqauto_file_types_and_scripts):
    PROCESSORS = {
        SequencingFileType.SAMPLE_SHEET: process_flowcells,
        SequencingFileType.ILLUMINA_FLOWCELL_QC: process_illuminate_qc,
        SequencingFileType.FASTQ: process_fastq,
        SequencingFileType.FASTQC: process_fastqc,
        SequencingFileType.BAM: process_bam,
        SequencingFileType.FLAGSTATS: process_flagstats,
        SequencingFileType.VCF: process_vcf,
        SequencingFileType.QC: process_qc
    }
    results = {}

    for file_type, script_name in seqauto_file_types_and_scripts:
        output_filename = os.path.join(seqauto_run.scan_resources_dir, "%s.txt" % name_from_filename(script_name))

        logging.info("create_resource_models: %s", SequencingFileType(file_type).label)

        lines = file_to_array(output_filename)
        processor = PROCESSORS[file_type]

        data = processor(seqauto_run, lines, results)
        results[file_type] = data

        print_data_state_stats(data)
        logging.info("-" * 20)
