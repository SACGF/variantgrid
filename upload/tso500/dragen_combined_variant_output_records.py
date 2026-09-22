"""
The database side of a TSO 500 pair's CombinedVariantOutput - everything in the file that is not a
variant. The file itself is read by upload.tso500.dragen_combined_variant_output_parser and its
splice calls become Variants in upload.tasks.import_dragen_tso500_combined_variant_output_task.

'[Analysis Details]' names the whole chain. 'Pair ID' is the pair's sample name, whose second
underscore-separated field is the patient's code (the lab's C-number) - the sequencing sample ID
leading it changes when the patient is re-sequenced, so the code is what one patient comes back
under and the whole pair ID is not (settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX reads it, so a lab
naming pairs some other way says so there). The ten-digit accession inside each sample ID is the
specimen, and the container suffix on it names that arm's extraction. Each level is resolved against what is
already there and created when absent, so a CVO arriving before anything has been accessioned
leaves a stub Patient holding only its code, a Specimen and two named Extractions - the patients
API keeps what the file lacks (name, DOB, sex, tissue, dates) and fills the stub in. In practice
the patient is pushed before sequencing starts, so the stub is the exception this tolerates rather
than the path it is built for.

The same two sample IDs are the join to everything else: DRAGEN writes them from the SampleSheet's
Sample_ID, so each is exactly a Sample.vcf_sample_name and exactly a SequencingSample.sample_name.
That is what links both arms' samples to their extraction and the splice VCF to its sequencing run,
in place of the filename matching seqauto does for a VCF it found on disk.

Nothing here fails an import: a chain that cannot be made leaves the variants and says why on the
import page (@see the caller), because a splice call is worth having whether or not the pair's
patient has been accessioned yet.
"""
import logging
import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import NamedTuple, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction
from django.db.models import Model
from django.utils import timezone
from django.utils.dateparse import parse_datetime

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.external_references import (
    ExternalReference,
    ResolvedReference,
    resolve_reference,
)
from patients.models import Extraction, Patient, Specimen, SpecimenMeasure
from patients.models_enums import MatchStatus, NucleicAcid, SpecimenMeasureType
from patients.serializers import upsert_specimen_measure
from seqauto.models import (
    SampleFromSequencingSample,
    SequencingRun,
    SequencingSample,
    VCFFromSequencingRun,
)
from snpdb.models import VCF, Sample
from upload.tso500.dragen_combined_variant_output_parser import (
    DNA_SAMPLE_ID,
    GENOMIC_INSTABILITY_SCORE,
    GIS,
    MSI,
    OUTPUT_DATE,
    OUTPUT_TIME,
    PAIR_ID,
    PERCENT_UNSTABLE_MSI_SITES,
    PLOIDY,
    RNA_SAMPLE_ID,
    TMB,
    TOTAL_TMB,
    TUMOR_FRACTION,
    USABLE_MSI_SITES,
    CombinedVariantOutputSection,
    get_section_values,
)

# The lab accession DRAGEN carries through into the sample IDs: ten digits identifying the specimen
# and a container suffix naming the extraction taken off it ('..._2600000001C', 'SA-C23755-2535115161C-D').
# The extraction half is what settings.PATIENT_EXTRACTION_SAMPLE_NAME_REGEX picks out of a VCF
# sample name on deployments with nothing upstream to quote an identifier
SAMPLE_ID_ACCESSION_PATTERN = re.compile(r"(?P<specimen>\d{10})(?P<container>[A-Za-z])")


class CombinedVariantOutputIdentityError(ValueError):
    """ The file names a chain we cannot make - a pair ID carrying no patient code, one arm's
        accession disagreeing with the other's, or a specimen already held by a different patient """


@dataclass(frozen=True)
class ArmIdentifiers:
    """ One nucleic acid arm of the pair, as '[Analysis Details]' names it """
    sample_id: str
    extraction_reference: str
    nucleic_acid: str


@dataclass(frozen=True)
class PairIdentifiers:
    pair_id: str
    patient_code: str
    specimen_reference: str
    dna: Optional[ArmIdentifiers]
    rna: Optional[ArmIdentifiers]

    @property
    def arms(self) -> list[ArmIdentifiers]:
        return [arm for arm in (self.dna, self.rna) if arm]


@dataclass(frozen=True)
class ResolvedPair:
    patient: Patient
    specimen: Specimen
    extractions: dict[str, Extraction]  # by sample id, so a Sample is reached by name alone

    def arm_extraction(self, arm: Optional[ArmIdentifiers]) -> Optional[Extraction]:
        if arm:
            return self.extractions.get(arm.sample_id)
        return None


class MeasureCall(NamedTuple):
    """ The lab's call on one measure and the policy that produced it. The call is None where the
        policy is set but the numbers cannot answer it - too few usable MSI sites - so the threshold
        that was applied is still recorded against the measure """
    call: Optional[str]
    threshold: str
    threshold_source: str


@dataclass(frozen=True)
class MeasureSource:
    """ Where one SpecimenMeasure is written in the file, and what turns it into a call """
    measure_type: str
    section: str
    key: str
    unit: Optional[str] = None
    call: Optional[Callable[[dict], Optional[MeasureCall]]] = None


def _value(values: dict, key: str) -> Optional[float]:
    try:
        return float(values[key])
    except (KeyError, TypeError, ValueError):
        return None


def band_call(value: float, bands: list) -> Optional[str]:
    """ The call whose lower bound the value reaches, off a setting's [(lower bound, call), ...] -
        SA Path's MSI is MSI-High >= 30%, MSI-Low >= 10%, MSS below that """
    for lower_bound, call in sorted(bands, key=lambda band: band[0], reverse=True):
        if value >= lower_bound:
            return call
    return None


def describe_bands(bands: list, unit: str) -> str:
    """ The policy in words, as the measure's threshold and the build form show it -
        'MSI-High >= 30%, MSI-Low >= 10%, MSS < 10%' """
    ordered = sorted(bands, key=lambda band: band[0], reverse=True)
    parts = [f"{call} >= {lower_bound:g}{unit}" for lower_bound, call in ordered[:-1]]
    if ordered:
        lowest = ordered[-1]
        parts.append(f"{lowest[1]} < {ordered[-2][0]:g}{unit}" if len(ordered) > 1 else lowest[1])
    return ", ".join(parts)


def msi_call(section_values: dict) -> Optional[MeasureCall]:
    """ The lab's MSI category off the percent of unstable sites, where the pair has enough usable
        sites for the percentage to mean anything. Both are the lab's policy, not vendor output, so
        an installation that has not set them gets the value and no call """
    min_usable_sites = settings.TSO500_MSI_MIN_USABLE_SITES
    bands = settings.TSO500_MSI_CALL_BANDS
    if min_usable_sites is None or not bands:
        return None
    percent = _value(section_values, PERCENT_UNSTABLE_MSI_SITES)
    if percent is None:
        return None
    threshold = f"{describe_bands(bands, '%')} unstable sites, needs >= {min_usable_sites} usable sites"
    source = "settings.TSO500_MSI_CALL_BANDS / settings.TSO500_MSI_MIN_USABLE_SITES"
    usable_sites = _value(section_values, USABLE_MSI_SITES)
    if usable_sites is None or usable_sites < min_usable_sites:
        return MeasureCall(None, threshold, source)
    return MeasureCall(band_call(percent, bands), threshold, source)


def tmb_call(section_values: dict) -> Optional[MeasureCall]:
    """ High / Low off the mutations per megabase, against the lab's bands """
    bands = settings.TSO500_TMB_CALL_BANDS
    if not bands:
        return None
    value = _value(section_values, TOTAL_TMB)
    if value is None:
        return None
    return MeasureCall(band_call(value, bands), describe_bands(bands, " mut/Mb"),
                       "settings.TSO500_TMB_CALL_BANDS")


# The five scalars the pair-level sections carry. '[GIS]' holds three: the score, and the tumour
# fraction and ploidy the caller estimated it from. MSI and TMB are the two the lab has a policy
# for - the rest are a number the report quotes and a pathologist reads
MEASURE_SOURCES = (
    MeasureSource(SpecimenMeasureType.TMB, TMB, TOTAL_TMB, "mut/Mb", call=tmb_call),
    MeasureSource(SpecimenMeasureType.MSI, MSI, PERCENT_UNSTABLE_MSI_SITES, "%", call=msi_call),
    MeasureSource(SpecimenMeasureType.GIS, GIS, GENOMIC_INSTABILITY_SCORE),
    MeasureSource(SpecimenMeasureType.TUMOUR_FRACTION, GIS, TUMOR_FRACTION),
    MeasureSource(SpecimenMeasureType.PLOIDY, GIS, PLOIDY),
)


def parse_patient_code(pair_id: str) -> str:
    """ The patient's code out of the pair's sample name, as the lab's naming writes it (@see the
        setting). A pair ID the regex does not read names a patient we cannot identify, and the
        whole pair ID is not it - one per sequencing of the patient would leave a patient per run """
    pattern = settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX
    if not pattern:
        return pair_id
    m = re.match(pattern, pair_id)
    if m is None:
        raise CombinedVariantOutputIdentityError(
            f"Pair ID '{pair_id}' does not match settings.TSO500_PAIR_ID_PATIENT_CODE_REGEX "
            f"('{pattern}') - no patient code to accession the pair against")
    return m.group("patient_code")


def parse_pair_identifiers(analysis_details: dict) -> Optional[PairIdentifiers]:
    """ The chain the file names, or None where it names none - an older module version writing no
        pair ID, or sample IDs carrying no accession """
    pair_id = analysis_details.get(PAIR_ID)
    arms = {}
    specimen_references = set()
    for key, nucleic_acid in ((DNA_SAMPLE_ID, NucleicAcid.DNA), (RNA_SAMPLE_ID, NucleicAcid.RNA)):
        sample_id = analysis_details.get(key)
        if not sample_id:
            continue
        m = SAMPLE_ID_ACCESSION_PATTERN.search(sample_id)
        if m is None:
            logging.info("'%s' carries no lab accession - no extraction named for it", sample_id)
            continue
        specimen_references.add(m.group("specimen"))
        arms[nucleic_acid] = ArmIdentifiers(sample_id=sample_id, extraction_reference=m.group(0),
                                            nucleic_acid=nucleic_acid)

    if not (pair_id and specimen_references):
        return None
    if len(specimen_references) > 1:
        raise CombinedVariantOutputIdentityError(
            f"'{pair_id}' pairs arms from different specimens: "
            f"{', '.join(sorted(specimen_references))}")

    return PairIdentifiers(pair_id=pair_id, patient_code=parse_patient_code(pair_id),
                           specimen_reference=specimen_references.pop(),
                           dna=arms.get(NucleicAcid.DNA), rna=arms.get(NucleicAcid.RNA))


def _resolve(model: type[Model], reference_id: str, user: User) -> Optional[Model]:
    """ What the identifier already names, or None to create it. Ambiguity is a human's to settle -
        the file is one lab's output, so guessing between two rows would attach a run to the wrong one """
    resolved = resolve_reference(model, ExternalReference(reference_id=reference_id), user)
    if resolved.status == MatchStatus.NEEDS_ATTENTION:
        raise CombinedVariantOutputIdentityError(resolved.error)
    return resolved.obj


@transaction.atomic
def resolve_pair(identifiers: PairIdentifiers, user: User) -> ResolvedPair:
    """ The pair's Patient / Specimen / Extraction rows, creating whichever are not there yet """
    patient = _resolve(Patient, identifiers.patient_code, user)
    if patient is None:
        patient = Patient.objects.create(patient_code=identifiers.patient_code)
        assign_permission_to_user_and_groups(user, patient)
        logging.info("Created %s from the patient code in Pair ID '%s'", patient,
                     identifiers.pair_id)

    specimen = _resolve(Specimen, identifiers.specimen_reference, user)
    if specimen is None:
        specimen = Specimen.objects.create(patient=patient,
                                           reference_id=identifiers.specimen_reference)
    elif specimen.patient != patient:
        raise CombinedVariantOutputIdentityError(
            f"Specimen '{identifiers.specimen_reference}' belongs to {specimen.patient}, "
            f"which is not the pair's patient {patient}")

    extractions = {}
    for arm in identifiers.arms:
        extraction = _resolve(Extraction, arm.extraction_reference, user)
        if extraction is None:
            extraction = Extraction.objects.create(specimen=specimen,
                                                   reference_id=arm.extraction_reference,
                                                   nucleic_acid_source=arm.nucleic_acid)
        elif extraction.specimen != specimen:
            raise CombinedVariantOutputIdentityError(
                f"Extraction '{arm.extraction_reference}' belongs to specimen "
                f"{extraction.specimen}, not {specimen}")
        extractions[arm.sample_id] = extraction

    return ResolvedPair(patient=patient, specimen=specimen, extractions=extractions)


def link_samples_to_extractions(resolved: ResolvedPair, user: User) -> int:
    """ Both arms' samples, wherever they were imported from - the file's sample IDs are the VCFs'
        sample names, so this settles a claim the arm files could only park """
    linked = 0
    for sample_id, extraction in resolved.extractions.items():
        reference = ExternalReference(reference_id=extraction.reference_id)
        match = ResolvedReference(reference, MatchStatus.MATCHED, obj=extraction)
        samples = Sample.filter_for_user(user, has_write_permission=True)
        for sample in samples.filter(vcf_sample_name=sample_id):
            if sample.apply_extraction_match(match):
                linked += 1
    return linked


def link_to_sequencing_run(vcf: VCF, sample: Sample, sample_id: str) -> Optional[SequencingRun]:
    """ The rows seqauto writes for a VCF it matched by path (@see
        upload.vcf.vcf_import.link_samples_and_vcfs_to_sequencing), against the arm the calls came
        off. A sample sheet that has not registered this sample yet means no link rows """
    sequencing_sample = SequencingSample.get_current().filter(sample_name=sample_id).order_by("-pk").first()
    if sequencing_sample is None:
        return None

    sequencing_run = sequencing_sample.sample_sheet.sequencing_run
    VCFFromSequencingRun.objects.update_or_create(vcf=vcf, defaults={"sequencing_run": sequencing_run})
    SampleFromSequencingSample.objects.update_or_create(sample=sample,
                                                        defaults={"sequencing_sample": sequencing_sample})
    return sequencing_run


def measured_date(analysis_details: dict):
    """ When the module wrote the file - the closest thing it carries to when the pair was measured """
    output_date = analysis_details.get(OUTPUT_DATE)
    if not output_date:
        return None
    output_time = analysis_details.get(OUTPUT_TIME) or "00:00:00"
    value = parse_datetime(f"{output_date} {output_time}")
    if value and timezone.is_naive(value):
        value = timezone.make_aware(value)
    return value


def write_specimen_measures(sections: dict[str, CombinedVariantOutputSection],
                            resolved: ResolvedPair, identifiers: PairIdentifiers,
                            user: User, method: str, date=None) -> list[SpecimenMeasure]:
    """ The pair's TMB, MSI, GIS, tumour fraction and ploidy. They describe the specimen, and the DNA
        arm is what produced them, so each is written against both. A re-analysis of the same
        specimen replaces the values rather than accumulating them (@see SpecimenMeasure.Meta) """
    extraction = resolved.arm_extraction(identifiers.dna)
    measures = []
    for measure_source in MEASURE_SOURCES:
        section_values = get_section_values(sections, measure_source.section)
        value = _value(section_values, measure_source.key)
        if value is None:
            continue
        measure_call = measure_source.call(section_values) if measure_source.call else None
        data = {
            "extraction": extraction,
            "measure_type": measure_source.measure_type,
            "value": value,
            "unit": measure_source.unit,
            "call": measure_call.call if measure_call else None,
            "threshold": measure_call.threshold if measure_call else None,
            "threshold_source": measure_call.threshold_source if measure_call else None,
            # The whole section, so which numbers this came off stays answerable
            "source_payload": section_values,
            "method": method,
            "measured_date": date,
        }
        measures.append(upsert_specimen_measure(resolved.specimen, data, user))
    return measures
