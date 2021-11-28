"""
Imports work by:

    import_patient_records:
        - Read the CSV with pandas, and create PatientRecord entries
        - Set whether a record is valid or needs manual intervention
        - Display what's going to happen with the records (ie good, bad etc)
        - Then, click SUBMIT button after review

    process_patient_records:
        - Do the actual conversion from PatientRecord into the various samples etc etc.

"""

import logging

import pandas as pd
from dateutil import parser
from django.utils import timezone
from guardian.shortcuts import get_objects_for_user

from annotation.phenotype_matching import bulk_patient_phenotype_matching
from library.guardian_utils import assign_permission_to_user_and_groups
from library.pandas_utils import df_nan_to_none
from patients.models import PatientColumns, PatientRecord, Specimen, Patient, \
    PatientModification, PatientRecordOriginType
from patients.models_enums import Sex, NucleicAcid, Mutation
from snpdb.models import Sample

UNKNOWN_STRING = 'UNKNOWN'  # Upper


def assign_patient_to_sample(patient_import, user, sample, patient, description, origin):
    """ Creates patient modification record """

    if sample.patient == patient:
        return

    old_patient = sample.patient

    sample.patient = patient
    sample.save()

    if old_patient:
        description += f" (previously patient was: {old_patient})"

    PatientModification.objects.create(patient=patient,
                                       user=user,
                                       description=description,
                                       origin=origin,
                                       patient_import=patient_import)


def assign_specimen_to_sample(patient_import, user, sample, specimen, description, origin):
    """ Creates patient modification record """

    if sample.specimen == specimen:
        return

    old_specimen = sample.specimen

    sample.specimen = specimen
    sample.save()

    if old_specimen:
        description += f" (previously specimen was: {old_specimen})"

    PatientModification.objects.create(patient=specimen.patient,
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
                d = parser.parse(date_string)
            except:
                message = f"{column}: Could not parse date '{date_string}'"
                validation_messages.append(message)

    return d


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
    else:
        if not nullable:
            message = f"{column}: Non-nullable boolean field was None"
            validation_messages.append(message)

    return value


def parse_choice(choices, row, column, validation_messages):
    """ Can be either the key or values in a choice (of any case) """
    choice_string = row[column]
    if choice_string is None:
        return None

    choice_string = choice_string.upper()

    choice_dict = dict(choices)
    if choice_string in choice_dict:
        return choice_string

    reverse_choice_dict = {b.upper(): a for a, b in choices}
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
        samples_qs = get_objects_for_user(user, 'snpdb.change_sample')

        if sample_id:
            kwargs = {"id": sample_id}
        elif sample_name:
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


def set_fields_if_blank(obj, field_values):
    """ returns true if change """
    changed = False
    for k, v in field_values.items():
        existing_value = getattr(obj, k)
        if existing_value:
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
    first_name = row[PatientColumns.PATIENT_FIRST_NAME]
    last_name = row[PatientColumns.PATIENT_LAST_NAME]
    date_of_birth = parse_date(row, PatientColumns.DATE_OF_BIRTH, validation_messages)
    date_of_death = parse_date(row, PatientColumns.DATE_OF_DEATH, validation_messages)
    sex = parse_choice(Sex.choices, row, PatientColumns.SEX, validation_messages)
    affected = parse_boolean(row, PatientColumns.AFFECTED, validation_messages)
    consanguineous = parse_boolean(row, PatientColumns.CONSANGUINEOUS, validation_messages)
    patient_phenotype = row[PatientColumns.PATIENT_PHENOTYPE]
    deceased = row[PatientColumns.DECEASED]

    # only set patient_deceased if deceased flag set but date_of_death not set
    if (deceased == 'Y' and date_of_death is None):
        patient_deceased = True
    elif (deceased == 'N' and date_of_death is None):
        patient_deceased = False
    else:
        patient_deceased = None

    specimen_reference_id = row[PatientColumns.SPECIMEN_REFERENCE_ID]
    specimen_description = row[PatientColumns.SPECIMEN_DESCRIPTION]
    specimen_collected_by = row[PatientColumns.SPECIMEN_COLLECTED_BY]
    specimen_collection_date = parse_date(row, PatientColumns.SPECIMEN_COLLECTION_DATE, validation_messages)
    specimen_received_date = parse_date(row, PatientColumns.SPECIMEN_RECEIVED_DATE, validation_messages)
    specimen_mutation_type = parse_choice(Mutation.choices, row, PatientColumns.SPECIMEN_MUTATION_TYPE, validation_messages)
    specimen_nucleic_acid_source = parse_choice(NucleicAcid.choices, row, PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE, validation_messages)
    specimen_age_at_collection = row[PatientColumns.SPECIMEN_AGE_AT_COLLECTION_DATE]

    sample_id = row[PatientColumns.SAMPLE_ID] or None
    sample_name = row[PatientColumns.SAMPLE_NAME]

    matched_sample = match_sample(user, sample_id, sample_name, validation_messages)
    if matched_sample:
        matched_sample_id = int(matched_sample.pk)
    else:
        matched_sample_id = None

    matched_patient = Patient.match(first_name, last_name, sex, date_of_birth, user=user)
    if matched_patient:
        created_patient = None
    else:
        created_patient = create_patient(patient_records.patient_import, first_name, last_name, sex, date_of_birth, user)

    patient = matched_patient or created_patient
    # update date of death if not already set
    # set _deceased back to False if DOD is being set
    # print('Deceased = %s AND Date of Death = %s' % (patient_deceased, date_of_death))
    description = None
    if date_of_death is not None:
        if patient.date_of_death != date_of_death:
            patient._deceased = None
            patient.date_of_death = date_of_death
            patient.save()
            description = "Updated patient date of death"

    else:
        # only set deceased if no date of death
        # update patient_deceased if not already set
        if patient_deceased:
            if patient._deceased != patient_deceased:
                patient.date_of_death = None
                patient._deceased = patient_deceased
                patient.save()
                description = "Updated patient as deceased = True"
        elif not patient_deceased:
            patient.date_of_death = None
            patient._deceased = patient_deceased
            patient.save()
            description = "Updated patient as deceased = False"
        else:
            # only clear date_of_death and _deceased if they are not currently NULL
            if patient.date_of_death is not None or patient._deceased is not None:
                patient.date_of_death = None
                patient._deceased = None
                patient.save()
                description = "Updated patient deceased and date_of_death to NULL values"

    # if patient has been modified, create a PatientModification record
    if description is not None:
        PatientModification.objects.create(patient=patient,
                                           user=user,
                                           description=description,
                                           origin=PatientRecordOriginType.UPLOADED_CSV,
                                           patient_import=patient_records.patient_import)

    # Fields that can change (ie not used to match)
    PATIENT_FIELDS = {"family_code": family_code,
                      "affected": affected,
                      "consanguineous": consanguineous}

    patient_modified = False
    for patient_field, field_value in PATIENT_FIELDS.items():
        if field_value:
            setattr(patient, patient_field, field_value)
            description = f"Set {patient_field} to {field_value}"
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
                                       "-- From Imported CSV on %s:" % timezone.now(),
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
    matched_specimen = None
    created_specimen = None
    # print("specimen_reference_id=%s" %(specimen_reference_id))
    if specimen_reference_id:
        # print("process specimen id=%s" %(specimen_reference_id))
        try:
            specimen = Specimen.objects.get(reference_id=specimen_reference_id)
            if specimen.patient != patient:

                msg = f"{specimen} had patient {patient}, tried to assign to patient {specimen.patient}"
                raise ValueError(msg)
            matched_specimen = specimen
        except Specimen.DoesNotExist:
            specimen = Specimen.objects.create(reference_id=specimen_reference_id,
                                               patient=patient)
            created_specimen = specimen

        field_values = {"reference_id": specimen_reference_id,
                        "description": specimen_description,
                        "collected_by": specimen_collected_by,
                        "patient": patient,
                        #tissue=tissue,
                        "collection_date": specimen_collection_date,
                        "received_date": specimen_received_date,
                        "mutation_type": specimen_mutation_type,
                        "nucleic_acid_source": specimen_nucleic_acid_source,
                        "_age_at_collection_date": specimen_age_at_collection}
        changed = set_fields_if_blank(specimen, field_values)
        if changed:
            # print("save specimen id=%s" %(specimen_reference_id))
            specimen.description = specimen_description
            specimen.collected_by = specimen_collected_by
            specimen.patient = patient
            specimen.collection_date = specimen_collection_date
            specimen.received_date = specimen_received_date
            specimen.mutation_type = specimen_mutation_type
            specimen.nucleic_acid_source = specimen_nucleic_acid_source
            specimen.age_at_collection = specimen_age_at_collection

            specimen.save()

    else:
        created_specimen = None

    if matched_sample:
        if specimen:
            matched_sample.specimen = specimen

        description = "Set during patient records import"
        assign_patient_to_sample(patient_records.patient_import, user, matched_sample, patient, description, origin=PatientRecordOriginType.UPLOADED_CSV)

    validation_message = '\n'.join(validation_messages)
    print(validation_message)

    PatientRecord.objects.create(patient_records=patient_records,
                                 record_id=record_id,
                                 validation_message=validation_message,
                                 matched_sample_id=matched_sample_id,
                                 matched_patient=matched_patient,
                                 matched_specimen=matched_specimen,
                                 created_patient=created_patient,
                                 created_specimen=created_specimen,
                                 sample_id=sample_id,
                                 sample_name=sample_name,
                                 patient_family_code=family_code,
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
                                 specimen_mutation_type=specimen_mutation_type,
                                 specimen_nucleic_acid_source=specimen_nucleic_acid_source,
                                 specimen_age_at_collection_date=specimen_age_at_collection)

    return patient_to_check_for_phenotype_match


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
        except:
            raise ude

    return df


def get_patient_record_imports_dataframe(f):
    df = pandas_read_encoded_csv(f, index_col=None, dtype=str)
    df = df_nan_to_none(df)
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def import_patient_records(patient_records):
    uploaded_file = patient_records.uploaded_file
    filename = uploaded_file.get_filename()
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
    patients_to_check_for_phenotype_matches = []
    for i, row in df.iterrows():
        logging.info("import_patient_records, process_record: %s", i)
        patient_with_phenotype = process_record(patient_records, i, row)
        if patient_with_phenotype:
            patients_to_check_for_phenotype_matches.append(patient_with_phenotype)
        items_processed += 1

    if patients_to_check_for_phenotype_matches:
        bulk_patient_phenotype_matching(patients_to_check_for_phenotype_matches)

    return items_processed
