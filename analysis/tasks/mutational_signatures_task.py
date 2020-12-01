from django.conf import settings
import celery
import logging

from analysis.models.enums import MinimisationResultType
from analysis.models.mutational_signatures import MutationalSignatureCalculator, \
    MutationalSignature, MutationalSignatureMinimisationResult, MutationalSignatureMutationCount
from library.genomics.calculate_cancer_mutation_signatures import MutationSignatures, invert_muttype
from snpdb.models import Sample, Variant
from snpdb.models.models_enums import ImportStatus
import numpy as np


@celery.task
def calculate_mutational_signature(sample_id):
    sample = Sample.objects.get(pk=sample_id)

    ms = MutationSignatures(sigdatafile=settings.MUTATIONAL_SIGNATURE_DATA_FILE,
                            siginfofile=settings.MUTATIONAL_SIGNATURE_INFO_FILE,
                            iterations=settings.MUTATIONAL_SIGNATURE_ITERATIONS,
                            sampling=settings.MUTATIONAL_SIGNATURE_SAMPLING_FRACTION)

    calculator = get_or_create_calculator_for_mutation_signatures(ms)
    # Check if exists and delete
    kwargs = {"sample": sample,
              "calculator": calculator}
    MutationalSignature.objects.filter(**kwargs).delete()
    mutational_signature = MutationalSignature.objects.create(**kwargs)

    # Extract data from VCF
    genome_build = sample.vcf.genome_build
    raw_mut_index_list = get_raw_mut_index_list_for_sample(genome_build.genome_fasta.fasta, ms.mut_index, sample)
    logging.info("Returned: %d muts", len(raw_mut_index_list))

    ms.UpdateMutIndexList(raw_mut_index_list)

    ms.calculate_signatures()

    # Store results to DB
    store_results_to_db(mutational_signature, ms.Results.full, MinimisationResultType.FULL)
    store_results_to_db(mutational_signature, ms.Results.bootstrapped, MinimisationResultType.BOOTSTRAPPED)

    # Mean
    all_data = np.array([msr.x for msr in ms.Results.full])
    mean = all_data.mean(axis=0).tolist()

    mutational_signature.import_status = ImportStatus.SUCCESS
    mutational_signature.mean = mean
    mutational_signature.summary = get_mutational_signature_summary(mean, threshold=0.02)
    mutational_signature.num_snps = ms.ObservedRawFrequency.sum()
    mutational_signature.save()

    labels = ms.get_mutation_type_labels()
    counts = ms.ObservedRawFrequency
    for mutation_type, count in zip(labels, counts):
        MutationalSignatureMutationCount.objects.create(mutational_signature=mutational_signature,
                                                        mutation_type=mutation_type,
                                                        count=count)


def get_or_create_calculator_for_mutation_signatures(ms):
    calculator, _ = MutationalSignatureCalculator.objects.get_or_create(name=settings.MUTATIONAL_SIGNATURE_CALCULATOR,
                                                                        version=ms.version,
                                                                        num_iterations=ms.ITERATIONS,
                                                                        sampling_fraction=ms.SAMPLING,
                                                                        signature_data_filename=settings.MUTATIONAL_SIGNATURE_DATA_FILE,
                                                                        minimisation_strategy=ms.minimization)
    return calculator


def store_results_to_db(mutational_signature, ms_results, result_type):
    for i, msr in enumerate(ms_results):
        MutationalSignatureMinimisationResult.objects.create(mutational_signature=mutational_signature,
                                                             result_type=result_type,
                                                             iteration=i,
                                                             solution_array=msr.x.tolist(),
                                                             fit_data=msr.fitdata.tolist(),
                                                             diff_data=msr.diffdata.tolist(),
                                                             ls_sum_diff=msr.ls_sum_diff,
                                                             la_sum_diff=msr.la_sum_diff,)


def get_raw_mut_index_list_for_sample(fasta_reference, mut_index, sample):
    qs = sample.get_variant_qs()
    qs = qs.filter(Variant.get_no_reference_q())
    snps_qs = qs.filter(locus__ref__length=1, alt__length=1)

    raw_mut_index_list = []
    for chrom, position, alt in snps_qs.values_list("locus__contig__name", "locus__position", "alt__seq"):

        start_ = position - 1
        end_ = position + 1
        context = fasta_reference[chrom][start_ - 1:end_]
        mutkey = f"{context}_{alt}"
        if mutkey not in mut_index:
            invkey = invert_muttype(mutkey)
            mutidx = mut_index[invkey]
        else:
            mutidx = mut_index[mutkey]
        raw_mut_index_list.append(mutidx)

    return raw_mut_index_list


def get_mutational_signature_summary(mean, threshold=0.02):
    """ threshold = minimum to display (default=0.02 = 2%)"""
    sorted_index = np.flip(np.argsort(mean), axis=0)
    signatures = []
    for i in sorted_index:
        val = mean[i]
        if val >= threshold:
            signatures.append("Sig %d: %.2f%%" % (i + 1, 100 * val))

    if signatures:
        summary_str = ", ".join(signatures)
    else:
        summary_str = "No clear signature detected."
    return summary_str
