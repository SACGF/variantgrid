"""
The patient records CSV import (PatientColumns are the columns):

    import_patient_records:
        - Read the CSV with pandas, then process_record each row in its own transaction, so a row that
          fails leaves nothing behind but its invalid PatientRecord and the other rows still import
        - A value that won't parse is a validation message on the row's PatientRecord; a row that can't be
          imported at all raises InvalidPatientRecord and becomes an invalid PatientRecord
        - Patients whose phenotype text changed are matched in bulk at the end
        - One pass: each row is parsed, matched and written straight away - there is no review step
        - When any row created an extraction, reconcile_pending_extractions runs so VCF samples parked on an
          extraction reference (eg 2600000001C) can find it

    process_record:
        - Match or create the patient, update its deceased state and the fields not used to match,
          then the specimen (a blank column leaves the stored value alone), its extraction and the sample
        - One row per extraction: the EXTRACTION_* columns name the arm (Specimen.get_or_create_extraction
          has the matching order), so the DNA and RNA arms of one block are two rows repeating the specimen
"""

import logging
from collections import Counter

import pandas as pd
from dateutil import parser
from django.db import transaction
from django.utils import timezone

from library.django_utils import ensure_timezone_aware
from library.guardian_utils import assign_permission_to_user_and_groups
from library.log_utils import report_exc_info
from library.pandas_utils import df_nan_to_none
from patients.enums import (
    NucleicAcid,
    PatientRecordMatchType,
    PatientRecordOriginType,
    Sex,
    TissueStatus,
)
from patients.models import (
    AmbiguousExtraction,
    Extraction,
    Patient,
    PatientColumns,
    PatientModification,
    PatientRecord,
    Specimen,
)
from patients.phenotype_matching import bulk_patient_phenotype_matching
from patients.tasks.extraction_matching_tasks import reconcile_pending_extractions
from snpdb.models import Sample

UNKNOWN_STRING = 'UNKNOWN'  # Upper
LEGACY_COLUMNS = {'Specimen Nucleic acid source (DNA/RNA)': PatientColumns.EXTRACTION_NUCLEIC_ACID_SOURCE}
OPTIONAL_COLUMNS = [PatientColumns.EXTRACTION_REFERENCE_ID, PatientColumns.EXTRACTION_DATE]


class InvalidPatientRecord(ValueError):
    """ A row that can't be imported - the message becomes its PatientRecord's validation message """


def assign_patient_to_sample(patient_import, user, sample, patient, origin):
    """ Creates patient modification record """

    if sample.patient == patient:
        return

    old_patient = sample.patient

    sample.patient = patient
    sample.save()

    description = f"Assigned patient to {sample=}"
    if old_patient:
        description += f" (previously patient was: {old_patient})"

    PatientModification.objects.create(patient=patient,
                                       user=user,
                                       description=description,
                                       origin=origin,
                                       patient_import=patient_import)


def assign_extraction_to_sample(patient_import, user, sample, extraction, origin):
    """ Creates patient modification record """

    if sample.extraction == extraction:
        return

    old_extraction = sample.extraction

    sample.extraction = extraction
    sample.save()

    description = f"Assigned extraction to {sample=}"
    if old_extraction:
        description += f" (previously extraction was: {old_extraction})"

    PatientModification.objects.create(patient=extraction.specimen.patient,
                                       user=user,
                                       description=description,
                                       origin=origin,
                                       patient_import=patient_import)


def parse_date(row, column, validation_messages):
    date_string = row[column]
    d = None
    if date_string and not pd.isnull(date_string):
        if date_string.upper() != UNKNOWN_STRING:
            try:
                d = parser.isoparse(date_string)
            except ValueError:
                try:
                    d = parser.parse(date_string, dayfirst=True)  # Column headers say DD-MM-YYYY
                except Exception:
                    message = f"{column}: Could not parse date '{date_string}'"
                    validation_messages.append(message)

    return d


def parse_int(row, column, validation_messages):
    value = row[column]
    if not value:
        return None

    try:
        return int(value)
    except ValueError:
        validation_messages.append(f"{column}: Could not parse whole number from '{value}'")
        return None


def parse_boolean(row, column, validation_messages, nullable=True):
    value = row[column]
    if value:
        lowercase_value = value.lower()
        TRUE_VALUES = ['true', 'y']
        FALSE_VALUES = ['false', 'n']

        if lowercase_value in TRUE_VALUES:
            value = True
        elif lowercase_value in FALSE_VALUES:
            value = False
        else:
            message = f"{column}: couldn't interpret boolean from '{value}'"
            validation_messages.append(message)
            value = None
    else:
        if not nullable:
            message = f"{column}: Non-nullable boolean field was None"
            validation_messages.append(message)

    return value


def parse_choice(choices, row, column, validation_messages):
    """ Key, label, or the label's first word when that is unambiguous (any case) - so TissueStatus
        takes 'A', 'Affected / lesional', or the 'Affected' its column header offers """
    choice_string = row[column]
    if choice_string is None:
        return None

    choice_string = choice_string.upper()

    choice_dict = dict(choices)
    if choice_string in choice_dict:
        return choice_string

    reverse_choice_dict = {label.upper(): key for key, label in choices}
    first_words = Counter(label.upper().split()[0] for label in choice_dict.values())
    for key, label in choices:
        first_word = label.upper().split()[0]
        if first_words[first_word] == 1:
            reverse_choice_dict.setdefault(first_word, key)
    value = reverse_choice_dict.get(choice_string)

    if value is None:
        valid = ','.join(list(choice_dict.keys()) + list(reverse_choice_dict.keys()))
        message = f"{column}: Could not parse choice '{choice_string}' (valid: {valid})"
        validation_messages.append(message)

    return value


def match_sample(user, sample_id, sample_name, validation_messages):
    sample = None
    if sample_id or sample_name:
        msg = None
        samples_qs = Sample.filter_for_user(user, has_write_permission=True)

        if sample_id:
            kwargs = {"id": sample_id}
        else:
            kwargs = {"name": sample_name}

        try:
            sample = samples_qs.get(**kwargs)
        except Sample.MultipleObjectsReturned:
            msg = "Matched multiple records!"
        except Sample.DoesNotExist:
            msg = f"Couldn't load sample '{sample_name}'"
            if sample_id:
                sample_without_permission = Sample.objects.filter(pk=sample_id)
            else:
                sample_without_permission = Sample.objects.filter(name=sample_name)

            if sample_without_permission.exists():
                msg += ". Samples exists with that name, but you don't have write access."

        if msg:
            logging.warning(msg)
            validation_messages.append(f"Match Sample({kwargs}): {msg}")

    return sample


def create_patient(patient_import, first_name, last_name, sex, date_of_birth, user):
    if first_name:
        first_name = first_name.upper()
    if last_name:
        last_name = last_name.upper()

    patient = Patient.objects.create(last_name=last_name,
                                     first_name=first_name,
                                     date_of_birth=date_of_birth,
                                     sex=sex or Sex.UNKNOWN)
    assign_permission_to_user_and_groups(user, patient)

    description = "Imported record"
    PatientModification.objects.create(patient=patient,
                                       user=user,
                                       description=description,
                                       origin=PatientRecordOriginType.UPLOADED_CSV,
                                       patient_import=patient_import)

    return patient


def _patient_details(patient, include_pk=True) -> str:
    """ The two patients in a clash often share name and sex, so add what else Patient.match compared """
    dob = patient.date_of_birth.strftime("%Y-%m-%d") if patient.date_of_birth else "blank"
    details = f"DOB: {dob}"
    if include_pk:
        details = f"Patient:{patient.pk}, {details}"
    return f"{patient} [{details}]"


def specimen_patient_clash_message(specimen, patient, patient_match_type, user) -> str:
    """ A specimen belongs to one patient. The specimen lookup isn't restricted to the user's patients,
        so the owner is only described when they can view it """
    owner = specimen.patient
    if owner.can_view(user):
        owner_description = _patient_details(owner)
    else:
        owner_description = f"a patient you don't have access to (Patient:{owner.pk})"

    if patient_match_type == PatientRecordMatchType.CREATED:
        # The failed row is rolled back, so this patient's pk won't exist
        row_patient = "matched no existing patient on name, sex and DOB, so would have created new patient " \
                      f"{_patient_details(patient, include_pk=False)}"
    else:
        row_patient = f"matched patient {_patient_details(patient)}"
    return f"Specimen '{specimen.reference_id}' belongs to {owner_description}. This row {row_patient}. " \
           "A specimen can only belong to one patient."


def set_fields_from_row(obj, field_values) -> bool:
    """ A blank CSV column is 'no answer' rather than 'clear it', so it leaves the stored value alone
        and a value replaces what is stored. Returns whether anything changed """
    changed = False
    for k, v in field_values.items():
        if v is not None and getattr(obj, k) != v:
            changed = True
            setattr(obj, k, v)

    return changed


def process_record(patient_records, record_id, row):
    patient_to_check_for_phenotype_match = None
    user = patient_records.user

    logging.info("row:")
    logging.info(row)

    validation_messages = []
    family_code = row[PatientColumns.PATIENT_FAMILY_CODE]
    patient_code = row[PatientColumns.PATIENT_CODE]
    first_name = row[PatientColumns.PATIENT_FIRST_NAME]
    last_name = row[PatientColumns.PATIENT_LAST_NAME]
    if not last_name:
        raise InvalidPatientRecord(f"{PatientColumns.PATIENT_LAST_NAME}: is blank, it is needed to match or create a patient")
    date_of_birth = parse_date(row, PatientColumns.DATE_OF_BIRTH, validation_messages)
    date_of_death = parse_date(row, PatientColumns.DATE_OF_DEATH, validation_messages)
    sex = parse_choice(Sex.choices, row, PatientColumns.SEX, validation_messages)
    affected = parse_boolean(row, PatientColumns.AFFECTED, validation_messages)
    consanguineous = parse_boolean(row, PatientColumns.CONSANGUINEOUS, validation_messages)
    patient_phenotype = row[PatientColumns.PATIENT_PHENOTYPE]
    deceased = parse_boolean(row, PatientColumns.DECEASED, validation_messages)

    # only set patient_deceased if deceased flag set but date_of_death not set
    if isinstance(deceased, bool) and date_of_death is None:
        patient_deceased = deceased
    else:
        patient_deceased = None

    specimen_reference_id = Specimen.normalize_reference_id(row[PatientColumns.SPECIMEN_REFERENCE_ID])
    specimen_description = row[PatientColumns.SPECIMEN_DESCRIPTION]
    specimen_collected_by = row[PatientColumns.SPECIMEN_COLLECTED_BY]
    specimen_collection_date = parse_date(row, PatientColumns.SPECIMEN_COLLECTION_DATE, validation_messages)
    specimen_received_date = parse_date(row, PatientColumns.SPECIMEN_RECEIVED_DATE, validation_messages)
    specimen_tissue_status = parse_choice(TissueStatus.choices, row, PatientColumns.SPECIMEN_TISSUE_STATUS, validation_messages)
    specimen_age_at_collection = parse_int(row, PatientColumns.SPECIMEN_AGE_AT_COLLECTION_DATE, validation_messages)
    extraction_reference_id = Extraction.normalize_reference_id(row[PatientColumns.EXTRACTION_REFERENCE_ID])
    extraction_nucleic_acid_source = parse_choice(NucleicAcid.choices, row, PatientColumns.EXTRACTION_NUCLEIC_ACID_SOURCE, validation_messages)
    extraction_date = parse_date(row, PatientColumns.EXTRACTION_DATE, validation_messages)

    sample_id = parse_int(row, PatientColumns.SAMPLE_ID, validation_messages)
    sample_name = row[PatientColumns.SAMPLE_NAME]

    sample = match_sample(user, sample_id, sample_name, validation_messages)
    try:
        patient, patient_match_type = Patient.match(first_name, last_name, sex, date_of_birth, user=user)
    except Patient.MultipleObjectsReturned as e:
        raise InvalidPatientRecord(str(e)) from e
    if patient is None:
        patient = create_patient(patient_records.patient_import, first_name, last_name, sex, date_of_birth, user)
        patient_match_type = PatientRecordMatchType.CREATED

    # update date of death if not already set
    # set _deceased back to False if DOD is being set
    # print('Deceased = %s AND Date of Death = %s' % (patient_deceased, date_of_death))
    description = None
    if date_of_death is not None:
        if patient.date_of_death != date_of_death:
            patient._deceased = None
            patient.date_of_death = date_of_death
            patient.save(check_patient_text_phenotype=False)
            description = "Updated patient date of death"

    else:
        # only set deceased if no date of death
        # update patient_deceased if not already set
        if patient_deceased:
            if patient._deceased != patient_deceased:
                patient.date_of_death = None
                patient._deceased = patient_deceased
                patient.save(check_patient_text_phenotype=False)
                description = "Updated patient as deceased = True"
        elif patient_deceased is False:
            patient.date_of_death = None
            patient._deceased = patient_deceased
            patient.save(check_patient_text_phenotype=False)
            description = "Updated patient as deceased = False"

    # if patient has been modified, create a PatientModification record
    if description is not None:
        PatientModification.objects.create(patient=patient,
                                           user=user,
                                           description=description,
                                           origin=PatientRecordOriginType.UPLOADED_CSV,
                                           patient_import=patient_records.patient_import)

    # Fields that can change (ie not used to match)
    PATIENT_FIELDS = {"family_code": family_code,
                      "patient_code": patient_code,
                      "affected": affected,
                      "consanguineous": consanguineous}

    patient_modified = False
    for patient_field, field_value in PATIENT_FIELDS.items():
        if field_value is not None:
            old_value = getattr(patient, patient_field)
            setattr(patient, patient_field, field_value)
            description = f"Set {patient_field} to {field_value}"
            if old_value != field_value:
                PatientModification.objects.create(patient=patient,
                                                   user=user,
                                                   description=description,
                                                   origin=PatientRecordOriginType.UPLOADED_CSV,
                                                   patient_import=patient_records.patient_import)
                patient_modified = True

    PHENOTYPE_FIELDS = {"phenotype": patient_phenotype}

    for phenotype_field, phenotype_value in PHENOTYPE_FIELDS.items():
        if phenotype_value:
            updated_phenotype = None
            existing_value = getattr(patient, phenotype_field)
            if existing_value:
                if phenotype_value not in existing_value:  # not already in there
                    # TODO: be more sophisticated?
                    phenotype_lines = [existing_value,
                                       f"-- From Imported CSV on {timezone.now()}:",
                                       phenotype_value]
                    updated_phenotype = "\n".join(phenotype_lines)
            else:
                updated_phenotype = phenotype_value

            if updated_phenotype:
                patient_to_check_for_phenotype_match = patient
                setattr(patient, phenotype_field, updated_phenotype)
                patient_modified = True

    if patient_modified:
        patient.save(check_patient_text_phenotype=False)  # Will do bulk at the end

    specimen = None
    extraction = None
    extraction_match_type = None
    if specimen_reference_id:
        specimen = Specimen.objects.filter(patient=patient, reference_id=specimen_reference_id).first()
        if specimen:
            specimen_match_type = PatientRecordMatchType.EXACT
        else:
            other_patients_specimen = Specimen.objects.filter(reference_id=specimen_reference_id).first()
            if other_patients_specimen:
                msg = specimen_patient_clash_message(other_patients_specimen, patient, patient_match_type, user)
                raise InvalidPatientRecord(msg)
            specimen = Specimen.objects.create(reference_id=specimen_reference_id,
                                               patient=patient)
            specimen_match_type = PatientRecordMatchType.CREATED

        field_values = {
            "description": specimen_description,
            "collected_by": specimen_collected_by,
            "collection_date": ensure_timezone_aware(specimen_collection_date),
            "received_date": ensure_timezone_aware(specimen_received_date),
            "tissue_status": specimen_tissue_status,
            "_age_at_collection_date": specimen_age_at_collection
        }
        if set_fields_from_row(specimen, field_values):
            specimen.save()

        existing_extraction_ids = set(specimen.extraction_set.values_list("pk", flat=True))
        try:
            extraction = specimen.get_or_create_extraction(extraction_nucleic_acid_source, sample=sample,
                                                           reference_id=extraction_reference_id,
                                                           extraction_date=ensure_timezone_aware(extraction_date))
        except AmbiguousExtraction as e:
            validation_messages.append(str(e))
        else:
            if extraction.pk in existing_extraction_ids:
                extraction_match_type = PatientRecordMatchType.EXACT
            else:
                extraction_match_type = PatientRecordMatchType.CREATED
    else:
        specimen_match_type = None
        if extraction_reference_id:
            validation_messages.append(f"{PatientColumns.EXTRACTION_REFERENCE_ID} '{extraction_reference_id}' "
                                       f"needs a {PatientColumns.SPECIMEN_REFERENCE_ID}")

    if sample:
        origin = PatientRecordOriginType.UPLOADED_CSV
        assign_patient_to_sample(patient_records.patient_import, user, sample, patient, origin)

        if extraction:
            assign_extraction_to_sample(patient_records.patient_import, user, sample, extraction, origin)

    validation_message = '\n'.join(validation_messages)
    valid = not validation_messages
    PatientRecord.objects.create(patient_records=patient_records,
                                 record_id=record_id,
                                 valid=valid,
                                 validation_message=validation_message,
                                 sample=sample,
                                 patient=patient,
                                 patient_match=patient_match_type,
                                 specimen=specimen,
                                 specimen_match=specimen_match_type,
                                 extraction=extraction,
                                 extraction_match=extraction_match_type,
                                 sample_identifier=sample_id,
                                 sample_name=sample_name,
                                 patient_family_code=family_code,
                                 patient_code=patient_code,
                                 patient_first_name=first_name,
                                 patient_last_name=last_name,
                                 date_of_birth=date_of_birth,
                                 date_of_death=date_of_death,
                                 sex=sex,
                                 affected=affected,
                                 consanguineous=consanguineous,
                                 _deceased=patient_deceased,
                                 patient_phenotype=patient_phenotype,
                                 specimen_reference_id=specimen_reference_id,
                                 specimen_description=specimen_description,
                                 specimen_collected_by=specimen_collected_by,
                                 specimen_collection_date=specimen_collection_date,
                                 specimen_received_date=specimen_received_date,
                                 specimen_tissue_status=specimen_tissue_status,
                                 specimen_age_at_collection_date=specimen_age_at_collection,
                                 extraction_reference_id=extraction_reference_id,
                                 extraction_nucleic_acid_source=extraction_nucleic_acid_source,
                                 extraction_date=row[PatientColumns.EXTRACTION_DATE])

    return patient_to_check_for_phenotype_match, extraction_match_type == PatientRecordMatchType.CREATED


def create_failed_patient_record(patient_records, record_id, row, exception):
    """ The row as text only, as any typed column may be what failed """
    PatientRecord.objects.create(patient_records=patient_records,
                                 record_id=record_id,
                                 valid=False,
                                 validation_message=str(exception),
                                 sample_name=row[PatientColumns.SAMPLE_NAME],
                                 patient_family_code=row[PatientColumns.PATIENT_FAMILY_CODE],
                                 patient_code=row[PatientColumns.PATIENT_CODE],
                                 patient_first_name=row[PatientColumns.PATIENT_FIRST_NAME],
                                 patient_last_name=row[PatientColumns.PATIENT_LAST_NAME] or "",
                                 patient_phenotype=row[PatientColumns.PATIENT_PHENOTYPE],
                                 specimen_reference_id=row[PatientColumns.SPECIMEN_REFERENCE_ID],
                                 specimen_description=row[PatientColumns.SPECIMEN_DESCRIPTION],
                                 specimen_collected_by=row[PatientColumns.SPECIMEN_COLLECTED_BY],
                                 extraction_reference_id=row[PatientColumns.EXTRACTION_REFERENCE_ID])


def pandas_read_encoded_csv(*args, **kwargs):
    """ Try opening as UTF8, then Windows code page 1252 (Western Excel) if that fails """
    try:
        df = pd.read_csv(*args, **kwargs)
    except UnicodeDecodeError as ude:
        if "encoding" in kwargs:  # Already there explicity - just fail...
            raise ude

        logging.warning("Reading CSV failed %s, retrying with windows code page", str(args))
        kwargs["encoding"] = 'cp1252'
        try:
            df = pd.read_csv(*args, **kwargs)
        except Exception as exc:
            raise ude from exc

    return df


def get_patient_record_imports_dataframe(f):
    """ A file written before #1870 (one extraction per specimen) has the nucleic acid column under its old
        Specimen name and no extraction reference / date - blank ones match the sample's extraction as before """
    df = pandas_read_encoded_csv(f, index_col=None, dtype=str)
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    df = df.rename(columns=LEGACY_COLUMNS)
    for column in OPTIONAL_COLUMNS:
        if column not in df.columns:
            df[column] = None
    return df_nan_to_none(df)


def import_patient_records(patient_records):
    file_upload = patient_records.file_upload
    filename = file_upload.get_filename()
    df = get_patient_record_imports_dataframe(filename)
    missing_columns = set(PatientColumns.COLUMNS) - set(df.columns)
    if missing_columns:
        expected = len(PatientColumns.COLUMNS)
        found_columns = set(PatientColumns.COLUMNS) & set(df.columns)
        found = len(found_columns)
        missing_str = ','.join([f'"{s}"' for s in missing_columns])
        msg = f"Invalid Patient Records Import file. Only {found} of {expected} columns supplied, missing: {missing_str}"
        raise ValueError(msg)

    items_processed = 0
    created_extraction = False
    patients_to_check_for_phenotype_matches = []
    # iterrows() would build each row as a str Series, turning df_nan_to_none's None back into NaN
    for i, row in enumerate(df.to_dict("records")):
        logging.info("import_patient_records, process_record: %s", i)
        # A row that fails leaves nothing behind but its invalid PatientRecord, and the other rows still import
        try:
            with transaction.atomic():
                patient_with_phenotype, row_created_extraction = process_record(patient_records, i, row)
            created_extraction |= row_created_extraction
        except Exception as e:
            if not isinstance(e, InvalidPatientRecord):
                report_exc_info(extra_data={"patient_records": patient_records.pk, "record_id": i})
            create_failed_patient_record(patient_records, i, row, e)
            patient_with_phenotype = None
        if patient_with_phenotype:
            patients_to_check_for_phenotype_matches.append(patient_with_phenotype)
        items_processed += 1

    if patients_to_check_for_phenotype_matches:
        bulk_patient_phenotype_matching(patients_to_check_for_phenotype_matches)

    if created_extraction:
        # A VCF imported before this file may have parked a claim on one of the new extractions
        reconcile_pending_extractions.delay()

    return items_processed
