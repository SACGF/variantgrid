from django.core.exceptions import ObjectDoesNotExist
from django.db import transaction
import itertools
import logging

from analysis.models.nodes.node_utils import reload_analysis_nodes
from snpdb.tasks.cohort_genotype_tasks import create_cohort_genotype_and_launch_task


@transaction.atomic
def vcf_replace_data(old_vcf, new_vcf):
    if new_vcf == old_vcf:
        raise ValueError("VCFs are the same!")

    # Find reasons why we can't do it
    old_vcf_samples = {s.vcf_sample_name: s for s in old_vcf.sample_set.all()}
    new_vcf_samples = {s.vcf_sample_name: s for s in new_vcf.sample_set.all()}

    missing_samples = set(old_vcf_samples) - set(new_vcf_samples)
    if missing_samples:
        raise ValueError(f"New VCF is missing the following samples: {missing_samples}")

    new_for_old = {}
    for vcf_sample_name, old_sample in old_vcf_samples.items():
        new_for_old[old_sample] = new_vcf_samples[vcf_sample_name]

    def get_new_samples(old_samples):
        return [new_for_old[s] for s in old_samples]

    # can't be FK or 1-to-1 fields on BOTH old and new
    ONE_TO_ONE_FIELDS = ["pathologytestordersample"]  # Won't move, just prevent moving
    MIGRATING_FIELDS = ["patient"]

    for old_sample, new_sample in new_for_old.items():
        for mf in ONE_TO_ONE_FIELDS + MIGRATING_FIELDS:
            old_val = getattr(old_sample, mf, None)
            if old_val:
                try:
                    new_val = getattr(new_sample, mf)
                    if new_val and new_val != old_val:
                        msg = f"Error: old '{old_sample}' and new '{new_sample}' samples both have field '{mf}' " \
                            f"old value: '{old_val}', new_value: '{new_val}'"
                        raise ValueError(msg)
                except ObjectDoesNotExist:
                    old_val.sample = new_sample
                    old_val.save()

    ## VCF stuff
    old_vcf.vcftag_set.update(vcf=new_vcf)

    # Core cohort (based off VCF)
    logging.info("VCF Cohort")
    num_old = len(old_vcf_samples)
    num_new = len(new_vcf_samples)

    if num_new > num_old:
        logging.info("Creating sub-cohort of %d old out of %d new samples", num_old, num_new)
        # make sub cohort representing old cohort
        sample_list = get_new_samples(old_vcf.cohort.get_samples())
        new_vcf_cohort = new_vcf.cohort.create_sub_cohort(new_vcf.user, sample_list)
    else:
        new_vcf_cohort = new_vcf.cohort

    old_new_cohorts = [(old_vcf.cohort, new_vcf_cohort)]

    # Handle sub-cohorts
    for old_sub_cohort in old_vcf.cohort.sub_cohort_set.all():
        sample_list = get_new_samples(old_sub_cohort.get_samples())
        new_sub_cohort = new_vcf.cohort.create_sub_cohort(new_vcf.user, sample_list)
        logging.info("Creating copy of sub cohort: %s", old_sub_cohort)
        old_new_cohorts.append((old_sub_cohort, new_sub_cohort))

    # Replace cohort in related objects
    logging.info("Replacing Cohort Related objects")
    for old_cohort, new_cohort in old_new_cohorts:
        new_for_old_cohort_sample = {}
        for old_cs in old_cohort.cohortsample_set.all():
            new_sample = new_for_old[old_cs.sample]
            new_cs = new_cohort.cohortsample_set.get(sample=new_sample)
            new_for_old_cohort_sample[old_cs] = new_cs

        # Do save instead of update so it clears internal stuff
        for node in old_cohort.cohortnode_set.all():
            node.cohort = new_cohort
            node.save()

        for cnzfc in old_cohort.cohortnodezygosityfilterscollection_set.all():
            cnzfc.cohort = new_cohort
            cnzfc.save()

            for cnzf in cnzfc.cohortnodezygosityfilter_set.all():
                cnzf.cohort_sample = new_for_old_cohort_sample[cnzf.cohort_sample]
                cnzf.save()

        for pedigree in old_cohort.pedigree_set.all():
            pedigree.cohort = new_cohort
            pedigree.save()

            for cspfr in pedigree.cohortsamplepedfilerecord_set.all():
                cspfr.cohort_sample = new_for_old_cohort_sample[cnzf.cohort_sample]
                cspfr.save()

        for trio in old_cohort.trio_set.all():
            trio.cohort = new_cohort
            trio.mother = new_for_old_cohort_sample[trio.mother]
            trio.father = new_for_old_cohort_sample[trio.father]
            trio.proband = new_for_old_cohort_sample[trio.proband]
            trio.save()

    # Now any further cohorts (not VCF or sub)
    # This swaps out CohortSample - so won't need to change
    old_vcf_based_cohorts = [on[0] for on in old_new_cohorts]
    sample_cohorts = set()

    logging.info("Non VCF Cohorts")
    for old_sample in old_vcf.sample_set.all():
        new_sample = new_for_old[old_sample]
        for cs in old_sample.cohortsample_set.exclude(cohort__in=old_vcf_based_cohorts):
            logging.info("%s %s => %s", cs.cohort.pk, old_sample, new_sample)
            cs.sample = new_sample
            cs.save()
            sample_cohorts.add(cs.cohort)

    for cohort in sample_cohorts:
        logging.info("Need to reload count tasks for %s", cohort)
        for cgc in cohort.cohortgenotypecollection_set.all():
            cgc.delete_related_objects()
            cgc.delete()
        create_cohort_genotype_and_launch_task(cohort, run_async=False)

    analyses_to_reload = set()
    new_vcf_based_cohorts = [on[1] for on in old_new_cohorts]  # Related objects now linked to this
    for cohort in itertools.chain(new_vcf_based_cohorts, sample_cohorts):
        for node in cohort.cohortnode_set.all():
            analyses_to_reload.add(node.analysis)

        for trio in cohort.trio_set.all():
            for node in trio.trionode_set.all():
                analyses_to_reload.add(node.analysis)

        for pedigree in cohort.pedigree_set.all():
            for node in pedigree.pedigreenode_set.all():
                analyses_to_reload.add(node.analysis)

    logging.info("Changing references to samples")
    for old_sample, new_sample in new_for_old.items():
        for mf in MIGRATING_FIELDS:
            v = getattr(old_sample, mf, None)
            if v:
                logging.info("Copying field %s (%s)", mf, v)
                setattr(new_sample, mf, v)
                new_sample.save()

        old_sample.sampletag_set.update(sample=new_sample)
        # Don't need to do anything else here
        old_sample.variantclassification_set.update(sample=new_sample)

        for related_field in ["genelistnode_set", "samplenode_set", "zygositynode_set"]:
            related = getattr(old_sample, related_field)
            for node in related.all():
                analyses_to_reload.add(node.analysis)
                node.sample = new_sample
                node.save()

    for analysis in analyses_to_reload:
        logging.info("Reloading: %s", analysis)
        reload_analysis_nodes(analysis.pk)

    logging.info("Done!")
