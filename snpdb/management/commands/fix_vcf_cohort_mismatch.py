import logging
import os

import cyvcf2
from django.core.management import BaseCommand

from snpdb.models import CohortGenotypeCollection


class Command(BaseCommand):

    def handle(self, *args, **options):
        for cgc in CohortGenotypeCollection.objects.filter(cohort__vcf__isnull=False):
            cohort = cgc.cohort
            vcf = cohort.vcf
            if not vcf.has_genotype:
                continue  # Will only have 1 sample
            filename = vcf.uploadedvcf.uploaded_file.get_filename()
            if os.path.exists(filename):
                cohort_samples_by_vcf_sample_name = {cs.sample.vcf_sample_name: cs for cs in
                                                     cohort.cohortsample_set.all()}
                reader = cyvcf2.VCF(filename)
                num_vcf_samples = len(reader.samples)
                if cg := cgc.cohortgenotype_set.first():
                    num_cgc_samples = len(cg.samples_zygosity)
                    if num_cgc_samples != num_vcf_samples:
                        logging.warning("VCF %d, VCF samples: %d, CohortGenotype samples: %d",
                                        vcf.pk, num_vcf_samples, num_cgc_samples)

                for i, sample_name in enumerate(reader.samples):
                    if cs := cohort_samples_by_vcf_sample_name.get(sample_name):
                        if cs.cohort_genotype_packed_field_index != i:
                            print(f"{sample_name} {cs.cohort_genotype_packed_field_index} -> {i}")
                    else:
                        logging.warning("VCF %d, sample_name '%s' deleted", vcf.pk, sample_name)
            else:
                logging.warning("VCF %d, filename '%s' not found", vcf.pk, filename)
