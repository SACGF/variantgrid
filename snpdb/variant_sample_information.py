from collections import Counter
from collections import defaultdict
from django.contrib.postgres.aggregates.general import StringAgg
from lazy import lazy

from annotation.models.models_phenotype_match import PATIENT_TPM_PATH, PATIENT_ONTOLOGY_TERM_PATH
from patients.models import Patient
from patients.models_enums import Zygosity
from snpdb.models import Variant, Sample, Locus, CohortGenotypeCollection
import pandas as pd


class VariantSampleInformation:

    def __init__(self, user, variant):
        locus_qs = Locus.objects.filter(pk=variant.locus.pk)
        values_qs = self._get_sample_values_for_variant_via_cohort_genotype(locus_qs)
        values_qs = self._cohort_genotype_to_sample_genotypes(values_qs)

        self.variant = variant
        self.num_samples = Sample.objects.count()
        self.user_sample_ids = set(Sample.filter_for_user(user).values_list("pk", flat=True))
        self.num_user_samples = len(self.user_sample_ids)

        locus_counter = defaultdict(Counter)
        locus_patients = defaultdict(set)
        num_observations = 0
        self.visible_rows = []
        for row in values_qs:
            pk = row["variant"]
            zygosity = row["zygosity"]
            if zygosity:
                locus_counter[pk][zygosity] += 1

            patient = row["sample__patient"]
            if patient:
                locus_patients[pk].add(patient)

            if variant.pk == pk:
                sample_id = row["sample"]
                if sample_id:
                    num_observations += 1

                    if sample_id in self.user_sample_ids:
                        self.visible_rows.append(row)

        has_hidden_samples = self.num_samples > self.num_user_samples
        self.hidden_samples_details = {}
        if has_hidden_samples:
            num_visible_observations = len(self.visible_rows)
            num_invisible_observations = num_observations - num_visible_observations
            self.hidden_samples_details = {"num_observations": num_observations,
                                           "num_visible_observations": num_visible_observations,
                                           "num_invisible_observations": num_invisible_observations}

        self.has_observations = num_observations > 0
        if self.has_observations:
            self.locus_counts_df = self._get_locus_counts(variant.pk, locus_counter)
        else:
            self.locus_counts_df = None

        self.patient_ids = locus_patients[variant.pk]

    @property
    def has_phenotype_match_graphs(self):
        patients_qs = Patient.objects.filter(pk__in=self.patient_ids)
        return patients_qs.filter(**{PATIENT_TPM_PATH + "__isnull": False}).exists()

    @lazy
    def has_unknown_zygosity(self):
        return 'Unknown' in self.locus_counts_df.columns

    @staticmethod
    def _get_sample_values_for_variant_via_cohort_genotype(locus_qs):
        """ This is the new, preferred way - as it gets all of the samples for a VCF at once """
        COHORT_PATH = "variant__cohortgenotype__collection__cohort"
        qs = locus_qs.filter(**{COHORT_PATH + "__vcf__isnull": False})  # Only VCFs
        return qs.values("variant",
                         "variant__alt",
                         "variant__cohortgenotype__collection",
                         "variant__cohortgenotype__collection__cohort__vcf__name",
                         "variant__cohortgenotype__samples_zygosity",
                         "variant__cohortgenotype__samples_allele_frequency",
                         "variant__cohortgenotype__samples_allele_depth",
                         "variant__cohortgenotype__samples_read_depth",
                         "variant__cohortgenotype__samples_phred_likelihood")

    @staticmethod
    def _cohort_genotype_to_sample_genotypes(values_qs):
        """ We're now joining to CohortGenotype - break up into samples so old code works """
        SAMPLE_ENRICHMENT_KIT_PATH = "samplefromsequencingsample__sequencing_sample__enrichment_kit"

        # turn a row of CohortGenotype into multiple rows of ObservedVariant values
        for cg_values in values_qs:
            variant = cg_values["variant"]
            cgc_id = cg_values["variant__cohortgenotype__collection"]
            vcf_name = cg_values["variant__cohortgenotype__collection__cohort__vcf__name"]
            samples_zygosity = cg_values["variant__cohortgenotype__samples_zygosity"]
            num_samples = len(samples_zygosity)
            empty = [-1] * num_samples
            samples_allele_frequency = cg_values["variant__cohortgenotype__samples_allele_frequency"] or empty
            samples_allele_depth = cg_values["variant__cohortgenotype__samples_allele_depth"] or empty
            samples_read_depth = cg_values["variant__cohortgenotype__samples_read_depth"] or empty
            samples_phred_likelihood = cg_values["variant__cohortgenotype__samples_phred_likelihood"] or empty

            cgc = CohortGenotypeCollection.objects.get(pk=cgc_id)
            samples_qs = cgc.cohort.get_samples()

            annotation_kwargs = {"ontology_terms": StringAgg("patient__" + PATIENT_ONTOLOGY_TERM_PATH + "__name", '|', distinct=True)}
            samples_qs = samples_qs.annotate(**annotation_kwargs)

            COPY_SAMPLE_FIELDS = ["id", "name", "patient", SAMPLE_ENRICHMENT_KIT_PATH]
            sample_values = samples_qs.values(*COPY_SAMPLE_FIELDS, *list(annotation_kwargs.keys()))

            for s_values in sample_values:
                sample_id = s_values["id"]
                i = cgc.get_array_index_for_sample_id(sample_id)
                zygosity = samples_zygosity[i]
                allele_frequency = samples_allele_frequency[i]
                allele_depth = samples_allele_depth[i]
                read_depth = samples_read_depth[i]
                phred_likelihood = samples_phred_likelihood[i]

                if zygosity == Zygosity.MISSING:
                    continue

                sample_genotype_values = {
                    "variant": variant,
                    "zygosity": zygosity,
                    "allele_frequency": allele_frequency,
                    "allele_depth": allele_depth,
                    "read_depth": read_depth,
                    "phred_likelihood": phred_likelihood,
                    "sample": sample_id,
                    "sample__vcf__name": vcf_name,
                }
                for k, v in s_values.items():
                    if k in COPY_SAMPLE_FIELDS:
                        k = "sample__" + k
                    sample_genotype_values[k] = v
                yield sample_genotype_values

    @staticmethod
    def _get_locus_counts(this_variant_id, locus_counter):
        variants_qs = Variant.objects.filter(pk__in=locus_counter.keys())
        variant_by_id = {v.pk: v for v in variants_qs}
        zygosity_display = dict(Zygosity.CHOICES)

        rows = []
        for variant_id, zygosity_counts in locus_counter.items():
            row = {"variant_id": variant_id,
                   "variant": variant_by_id[variant_id],
                   "description": '',
                   "Total": 0}

            for z in zygosity_display.values():
                row[z] = 0

            if variant_id == this_variant_id:
                row["sort_order"] = 1
                row["description"] = 'This variant'
            else:
                v = variant_by_id[variant_id]
                if v.is_reference:
                    row["sort_order"] = 0
                else:
                    row["sort_order"] = sum(map(ord, v.alt.seq))

            for z, count in zygosity_counts.items():
                zyg_display = zygosity_display[z]
                row[zyg_display] = row.get(zyg_display, 0) + count
                row["Total"] += count

            rows.append(row)

        df = pd.DataFrame.from_records(rows).sort_values("sort_order")
        return df
