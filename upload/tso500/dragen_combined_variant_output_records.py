"""
The database side of a TSO 500 pair's CombinedVariantOutput - everything in the file that is not a
variant. The file itself is read by upload.tso500.dragen_combined_variant_output_parser and the import
step is upload.tasks.import_dragen_tso500_combined_variant_output_task.

The analysis itself - '[Analysis Details]' and the '[TMB]', '[MSI]' and '[GIS]' scalars - is one
seqauto.models.DragenTSO500CombinedVariantOutput row per (run, pair), written by
write_combined_variant_output. The run is the upload's 'sequencing_run' metadata, falling back to the
run whose current sample sheet names one of the pair's sample IDs.

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

Nothing here fails an import: a chain that cannot be made still writes the analysis, with its specimen
claim parked, and says why on the import page (@see the caller).
"""
import logging
import re
from dataclasses import dataclass
from typing import Optional

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
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import MatchStatus, NucleicAcid
from seqauto.models import (
    DragenTSO500CombinedVariantOutput,
    SampleFromSequencingSample,
    SequencingRun,
    SequencingSample,
    VCFFromSequencingRun,
)
from snpdb.models import VCF, Sample
from upload.tso500.dragen_combined_variant_output_parser import (
    CODING_REGION_SIZE,
    DNA_SAMPLE_ID,
    GENOMIC_INSTABILITY_SCORE,
    GIS,
    MODULE_VERSION,
    MSI,
    OUTPUT_DATE,
    OUTPUT_TIME,
    PAIR_ID,
    PASSING_ELIGIBLE_VARIANTS,
    PERCENT_UNSTABLE_MSI_SITES,
    PIPELINE_VERSION,
    PLOIDY,
    RNA_SAMPLE_ID,
    TMB,
    TOTAL_MSI_SITES_UNSTABLE,
    TOTAL_TMB,
    TUMOR_FRACTION,
    USABLE_MSI_SITES,
    CombinedVariantOutputSection,
    get_analysis_details,
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


def parse_output_datetime(analysis_details: dict):
    """ When the module wrote the file - the closest thing it carries to when the pair was measured """
    output_date = analysis_details.get(OUTPUT_DATE)
    if not output_date:
        return None
    output_time = analysis_details.get(OUTPUT_TIME) or "00:00:00"
    value = parse_datetime(f"{output_date} {output_time}")
    if value and timezone.is_naive(value):
        value = timezone.make_aware(value)
    return value


def sequencing_run_for_sample_ids(sample_ids: list[str]) -> Optional[SequencingRun]:
    """ The run whose current sample sheet names one of the pair's sample IDs - the run a CVO sent without
        'sequencing_run' metadata came off """
    for sample_id in sample_ids:
        sequencing_sample = SequencingSample.get_current().filter(sample_name=sample_id).order_by("-pk").first()
        if sequencing_sample:
            return sequencing_sample.sample_sheet.sequencing_run
    return None


def _float(values: dict, key: str) -> Optional[float]:
    try:
        return float(values[key])
    except (KeyError, TypeError, ValueError):
        return None


def _int(values: dict, key: str) -> Optional[int]:
    value = _float(values, key)
    return None if value is None else int(value)


def combined_variant_output_values(sections: dict[str, CombinedVariantOutputSection]) -> dict:
    """ The row's columns off the file, bar its key and links """
    analysis_details = get_analysis_details(sections)
    tmb = get_section_values(sections, TMB)
    msi = get_section_values(sections, MSI)
    gis = get_section_values(sections, GIS)
    return {
        "dna_sample_name": analysis_details.get(DNA_SAMPLE_ID) or "",
        "rna_sample_name": analysis_details.get(RNA_SAMPLE_ID) or "",
        "output_datetime": parse_output_datetime(analysis_details),
        "module_version": analysis_details.get(MODULE_VERSION) or "",
        "pipeline_version": analysis_details.get(PIPELINE_VERSION) or "",
        "total_tmb": _float(tmb, TOTAL_TMB),
        "coding_region_size_mb": _float(tmb, CODING_REGION_SIZE),
        "passing_eligible_variants": _int(tmb, PASSING_ELIGIBLE_VARIANTS),
        "usable_msi_sites": _int(msi, USABLE_MSI_SITES),
        "total_msi_sites_unstable": _int(msi, TOTAL_MSI_SITES_UNSTABLE),
        "percent_unstable_msi_sites": _float(msi, PERCENT_UNSTABLE_MSI_SITES),
        "genomic_instability_score": _float(gis, GENOMIC_INSTABILITY_SCORE),
        "tumor_fraction": _float(gis, TUMOR_FRACTION),
        "ploidy": _float(gis, PLOIDY),
    }


def write_combined_variant_output(sections: dict[str, CombinedVariantOutputSection], pair_id: str,
                                  sequencing_run: Optional[SequencingRun], sequencing_run_name: str,
                                  user: User, resolved: Optional[ResolvedPair] = None,
                                  specimen_reference: Optional[str] = None,
                                  parked_error: Optional[str] = None,
                                  file_upload=None) -> DragenTSO500CombinedVariantOutput:
    """ The pair's analysis on this run - a re-analysis of the run replaces it. It claims the specimen
        resolve_pair made, or parks the claim with why there is none """
    cvo, _ = DragenTSO500CombinedVariantOutput.objects.update_or_create(
        sequencing_run_name=sequencing_run_name, pair_id=pair_id,
        defaults={
            **combined_variant_output_values(sections),
            "sequencing_run": sequencing_run,
            "specimen_reference": specimen_reference or "",
            "file_upload": file_upload,
            "user": user,
        })
    cvo.link_sequencing_samples()
    cvo.link_samples()
    if resolved:
        reference = ExternalReference(reference_id=resolved.specimen.reference_id)
        cvo.apply_specimen_match(ResolvedReference(reference, MatchStatus.MATCHED, obj=resolved.specimen),
                                 save=False)
    else:
        cvo.park_specimen_claim(parked_error or f"'{pair_id}' names no specimen", save=False)
    cvo.save()
    return cvo
