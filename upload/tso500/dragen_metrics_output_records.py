"""
The database side of a TSO 500 run's MetricsOutput - the LibraryQC rows it writes. The file itself
is read by upload.tso500.dragen_metrics_output_parser and the import step is
upload.tasks.import_dragen_tso500_metrics_output_task.

Entry point is write_library_qc: one row per (run, pair, QC category). The file is written per
sequencing run and names that run nowhere, so the run comes off the upload's 'sequencing_run'
metadata - it is what makes a re-sequenced pair a second row rather than an overwrite, and a
re-analysis of the same run replaces its rows. The SequencingRun is linked by name where seqauto has
registered it, and each row's arm to its SequencingSample where the run's current sheet carries the
TSO500 Pair_ID / Sample_Type columns (seqauto.models.sequencing_sample_for_pair) - both left null
otherwise, for reconcile_pending_extractions to fill in once the sheet is registered.

A column is a Pair ID, not a sample, so it names no extraction. It claims the Specimen whose
accession it ends in (PAIR_ID_ACCESSION_PATTERN) and never creates one - accessioning a case is the
CombinedVariantOutput's job - so QC that arrives first parks its claim and
patients.tasks.extraction_matching_tasks.reconcile_pending_extractions attaches it once the specimen
lands.

The lab writes the Pair ID either as the whole pair sample name or as the patient's code alone, so a
column often carries no accession. The run is what supplies it then: the sheet's SequencingSample
names carry the code and the accession together ('1_TSO_DNAHRD_C17817_2517114977C_B4'), so the
accession comes off the linked arm's name, or - on a sheet posted without Pair_ID data - off the row
the column's patient code (settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX, the CVO's reader) finds. A
column nothing on the run is named for - a control, another assay sharing the flowcell - still gets
its rows, with the claim parked saying so.
"""
import logging
import re
from typing import Optional

from django.contrib.auth.models import User

from patients.external_references import ExternalReference, resolve_reference
from patients.models import Specimen
from seqauto.models import (
    LibraryQC,
    SequencingRun,
    SequencingSample,
    sequencing_sample_for_pair,
)
from upload.tso500.dragen_combined_variant_output_parser import (
    OUTPUT_DATE,
    OUTPUT_TIME,
)
from upload.tso500.dragen_combined_variant_output_records import (
    SAMPLE_ID_ACCESSION_PATTERN,
    CombinedVariantOutputIdentityError,
    measured_date,
    parse_patient_code,
)
from upload.tso500.dragen_metrics_output_parser import (
    CATEGORY_NUCLEIC_ACID,
    HEADER,
    LibraryQCMetrics,
    get_workflow_version,
)

METHOD = "DRAGEN TSO500 MetricsOutput"

# A Pair ID ends in the specimen's ten-digit lab accession - '5_C0000001_FCUP_2600000001'. Unlike a
# sample ID it carries no container suffix, so it names the specimen rather than either extraction
PAIR_ID_ACCESSION_PATTERN = re.compile(r"(?P<specimen>\d{10})$")


def metrics_method(sections) -> str:
    """ The tool and version that judged the library, as SpecimenMeasure.method words it """
    version = get_workflow_version(sections)
    return f"{METHOD} {version}" if version else METHOD


def metrics_measured_date(sections):
    """ When the module wrote the file - '[Header]' spells Output Date / Time as the CVO does """
    if section := sections.get(HEADER):
        values = section.values
        return measured_date({OUTPUT_DATE: values.get(OUTPUT_DATE),
                              OUTPUT_TIME: values.get(OUTPUT_TIME)})
    return None


def specimen_reference(pair_id: str) -> Optional[str]:
    """ The specimen the pair was taken from - the accession DRAGEN carries through into its Pair IDs """
    if m := PAIR_ID_ACCESSION_PATTERN.search(pair_id):
        return m.group("specimen")
    return None


def specimen_reference_from_sequencing_samples(sequencing_samples: list[SequencingSample]) -> Optional[str]:
    """ The accession off the linked arms' names - the sheet's sample names carry the code and the
        accession together """
    for sequencing_sample in sequencing_samples:
        if m := SAMPLE_ID_ACCESSION_PATTERN.search(sequencing_sample.sample_name):
            return m.group("specimen")
    return None


def specimen_reference_from_run(sequencing_run: Optional[SequencingRun],
                               pair_id: str) -> Optional[str]:
    """ The accession a run's sample sheet holds for a pair named only by its patient code, on a sheet
        posted without Pair_ID data: the sheet row whose name carries the code """
    if sequencing_run is None:
        return None
    try:
        patient_code = parse_patient_code(pair_id)
    except CombinedVariantOutputIdentityError:
        return None
    if not patient_code:
        return None
    sequencing_samples = SequencingSample.get_current().filter(
        sample_sheet__sequencing_run=sequencing_run, sample_name__contains=patient_code)
    for sequencing_sample in sequencing_samples.order_by("pk"):
        if m := SAMPLE_ID_ACCESSION_PATTERN.search(sequencing_sample.sample_name):
            return m.group("specimen")
    return None


def write_library_qc(libraries: list[LibraryQCMetrics], user: User, sequencing_run_name: str,
                     method: str, date=None, file_upload=None) -> list[LibraryQC]:
    """ Every category each pair has QC for. Both arms are in the one column, so a pair sequenced on
        both gets the five DNA categories and RNA - an arm it does not have is all NA, which is no QC
        rather than a failure, and is not written """
    sequencing_run = SequencingRun.objects.filter(name=sequencing_run_name).first()
    rows = []
    for library in libraries:
        arm_samples = {arm: sequencing_sample_for_pair(sequencing_run, library.pair_id, arm)
                       for arm in set(CATEGORY_NUCLEIC_ACID.values())}
        reference_id = specimen_reference(library.pair_id) or \
            specimen_reference_from_sequencing_samples([ss for ss in arm_samples.values() if ss]) or \
            specimen_reference_from_run(sequencing_run, library.pair_id)
        parked_error = None
        if reference_id is None:
            parked_error = (f"'{library.pair_id}' carries no lab accession, and nothing on "
                            f"sequencing run '{sequencing_run_name}' is named for it")
            logging.warning("%s - its library QC is parked", parked_error)
        else:
            resolved = resolve_reference(Specimen, ExternalReference(reference_id=reference_id), user)
        for category, metrics in library.categories.items():
            passed = library.category_passed(category)
            if passed is None:
                continue  # an arm this pair does not have - every value in the section is NA
            library_qc, _ = LibraryQC.objects.update_or_create(
                sequencing_run_name=sequencing_run_name, pair_id=library.pair_id,
                category=category,
                defaults={
                    "sequencing_run": sequencing_run,
                    "sequencing_sample": arm_samples[CATEGORY_NUCLEIC_ACID[category]],
                    "specimen_reference": reference_id or "",
                    "nucleic_acid": CATEGORY_NUCLEIC_ACID[category],
                    "passed": passed,
                    "completed": library.completed,
                    "metrics": metrics,
                    "method": method,
                    "measured_date": date,
                    "file_upload": file_upload,
                    "user": user,
                })
            if parked_error:
                library_qc.park_specimen_claim(parked_error)
            else:
                library_qc.apply_specimen_match(resolved)
            rows.append(library_qc)
    return rows
