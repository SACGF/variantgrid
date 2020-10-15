from typing import Optional

import logging

from patients.models_enums import Zygosity, Sex
from snpdb.models import Sample, VCFFilter, Specimen, SampleGenotype
from snpdb.models.models_variant import Variant
from classification.autopopulate_evidence_keys.evidence_from_variant import AutopopulateData
from classification.enums import SpecialEKeys, ValidationCode
from classification.models import VCDataDict


def get_zygosity(sample_genotype: SampleGenotype) -> Optional[str]:
    """ returns Homozygous/het/hem """
    if sample_genotype.zygosity in (Zygosity.HOM_ALT, Zygosity.HOM_REF):
        return "Homozygous"
    elif sample_genotype.zygosity == Zygosity.HET:
        has_chr = sample_genotype.sample.genome_build.reference_fasta_has_chr
        CHR_X = "chrX" if has_chr else 'X'
        if sample_genotype.variant.locus.chrom == CHR_X:
            patient = sample_genotype.sample.patient
            if patient and patient.sex == Sex.MALE:
                return "Hemizygous"
        return "Heterozygous"
    return None


def negative_as_none(value):
    if value is not None:
        if value < 0:
            return None
    return value


def get_evidence_fields_for_sample_and_patient(variant: Variant, sample: Sample) -> AutopopulateData:
    data = AutopopulateData("sample")
    data[SpecialEKeys.SAMPLE_ID] = sample.name

    patient = sample.patient
    if patient:
        data[SpecialEKeys.PATIENT_ID] = patient.external_pk or patient.pk
        if patient.family_code:
            data[SpecialEKeys.FAMILY_ID] = patient.family_code

    if sample.specimen:
        data.update(get_evidence_fields_for_specimen(sample.specimen))

    if variant and sample:
        if sample_genotype := sample.get_genotype(variant):
            data[SpecialEKeys.ZYGOSITY] = get_zygosity(sample_genotype)
            # SampleGenotype has unknown as -1
            data[SpecialEKeys.ALLELE_DEPTH] = negative_as_none(sample_genotype.allele_depth)
            data[SpecialEKeys.READ_DEPTH] = negative_as_none(sample_genotype.read_depth)
            data[SpecialEKeys.PHRED_LIKELIHOOD] = negative_as_none(sample_genotype.phred_likelihood)
            data[SpecialEKeys.GENOTYPE_QUALITY] = negative_as_none(sample_genotype.genotype_quality)
            data[SpecialEKeys.ALLELE_FREQUENCY] = negative_as_none(sample_genotype.allele_frequency)
            data[SpecialEKeys.VCF_FILTER] = ", ".join(sample_genotype.get_vcf_filters())
        else:
            logging.error("%s does not have ObservedVariant record for %s", sample, variant)

    try:
        sfss = sample.samplefromsequencingsample
        data[SpecialEKeys.CAPTURE_METHOD] = sfss.sequencing_sample.enrichment_kit.name  # version?

        sequencing_model = sfss.sequencing_run.sequencer.sequencing_model
        if sequencing_model.manufacturer.name == 'Illumina':
            model_lower = sequencing_model.model.lower()
            if "hiseq" in model_lower:
                sequencing_platform = "Illumina_HiSeq"
            elif "miseq" in model_lower:
                sequencing_platform = "Illumina_MiSeq"
            elif "nextseq" in model_lower:
                sequencing_platform = "Illumina_NextSeq"
            else:
                # couldn't match to a specific platform, include the raw text
                sequencing_platform = model_lower

            data[SpecialEKeys.SEQUENCING_PLATFORM] = sequencing_platform
    except:
        pass

    if sample.patient:
        data.update(get_evidence_fields_for_patient(sample.patient))

    return data


def get_evidence_fields_for_specimen(specimen: Specimen) -> AutopopulateData:
    data = AutopopulateData("specimen")
    data[SpecialEKeys.SPECIMEN_ID] = specimen.reference_id
    data[SpecialEKeys.ALLELE_ORIGIN] = specimen.get_mutation_type_display()
    data[SpecialEKeys.SAMPLE_TYPE] = str(specimen.tissue)
    data[SpecialEKeys.AGE] = specimen.age_at_collection_date
    data[SpecialEKeys.NUCLEIC_ACID_SOURCE] = specimen.get_nucleic_acid_source_display()

    return data


def get_evidence_fields_for_patient(patient) -> AutopopulateData:
    data = AutopopulateData("patient")
    sex = None
    if patient.sex == Sex.MALE:
        sex = "male"
    elif patient.sex == Sex.FEMALE:
        sex = 'female'

    data[SpecialEKeys.SEX] = sex
    ancestry = list(patient.patientpopulation_set.all().values_list("population", flat=True))
    if ancestry:
        data[SpecialEKeys.ANCESTRY] = ancestry
    return data
