import logging
from collections import Counter, defaultdict
from functools import cached_property
from typing import Optional

from django.conf import settings
from django.contrib.postgres.aggregates.general import StringAgg
from django.db.models import Q, TextField
from django.urls import reverse

from annotation.models.models import VariantAnnotation
from annotation.models.models_phenotype_match import PATIENT_ONTOLOGY_TERM_PATH
from library.log_utils import log_traceback
from library.unit_percent import format_af
from ontology.models import OntologyService
from patients.models_enums import Zygosity
from snpdb.models import Variant, Sample, Locus, CohortGenotypeCollection, CohortGenotype, VCFFilter
from upload.models import ModifiedImportedVariant

SAMPLE_ENRICHMENT_KIT_PATH = "samplefromsequencingsample__sequencing_sample__enrichment_kit__name"
# Zygosity.CHOICES uses '?' for unknown, which doesn't make a great JSON key or column header
LOCUS_COUNT_ZYGOSITY_LABELS = [(Zygosity.HOM_REF, "REF"), (Zygosity.HET, "HET"),
                               (Zygosity.HOM_ALT, "HOM_ALT"), (Zygosity.UNKNOWN_ZYGOSITY, "Unknown")]


class VariantZygosityCounts:
    """ Permission-scoped per-variant sample/zygosity lookup - the reusable core extracted
        from VariantSampleGenotypes (no presentation concerns). Shared by the
        Variantopedia sample-information API (subclass below) and the Beacon
        g_variants endpoint (§5.2 stage 2). Gates by Sample.filter_for_user (anonymous ->
        public group), so presence/counts are only ever asserted from samples the requester
        may read. """

    def __init__(self, user, variant):
        locus_qs = Locus.objects.filter(pk=variant.locus.pk)
        values_qs = self._get_sample_values_for_variant_via_cohort_genotype(locus_qs)
        values_qs = self._cohort_genotype_to_sample_genotypes(values_qs)

        self.variant = variant
        # Builds can share contigs (GRCh37/38 share MT and some unplaced scaffolds) in which case both builds
        # have the same variant. Scope to every build the variant is in - scoping to one would count the other
        # build's samples as observations without ever making them visible rows, reporting them as hidden.
        self.genome_builds = sorted(variant.genome_builds, key=lambda gb: gb.name)
        self.num_samples = Sample.objects.filter(vcf__genome_build__in=self.genome_builds).count()
        visible_samples_qs = Sample.filter_for_user(user).filter(vcf__genome_build__in=self.genome_builds)
        self.user_sample_ids = set(visible_samples_qs.values_list("pk", flat=True))
        self.num_user_samples = len(self.user_sample_ids)

        self.visible_zygosity_counter = Counter()  # raw {zygosity: count} for visible samples
        self.locus_counter = defaultdict(Counter)
        self.num_observations = 0
        self.visible_rows = []
        for row in values_qs:
            pk = row["variant"]
            zygosity = row.get("zygosity", Zygosity.UNKNOWN_ZYGOSITY)
            self.locus_counter[pk][zygosity] += 1

            if variant.pk == pk:
                if sample_id := row.get("sample"):
                    self.num_observations += 1

                    if sample_id in self.user_sample_ids:
                        self.visible_rows.append(row)
                        self.visible_zygosity_counter[zygosity] += 1

    @cached_property
    def genome_builds_str(self) -> str:
        return ", ".join(str(gb) for gb in self.genome_builds)

    @property
    def has_hidden_samples(self) -> bool:
        """ Are there samples in these builds the user can't see (regardless of this variant) """
        return self.num_samples > self.num_user_samples

    @property
    def num_visible_observations(self) -> int:
        return len(self.visible_rows)

    @property
    def num_invisible_observations(self) -> int:
        return self.num_observations - self.num_visible_observations

    @property
    def num_visible_het(self) -> int:
        return self.visible_zygosity_counter[Zygosity.HET]

    @property
    def num_visible_hom_alt(self) -> int:
        return self.visible_zygosity_counter[Zygosity.HOM_ALT]

    @property
    def num_visible_alt(self) -> int:
        """ Visible het + hom_alt observations - the presence/count a Beacon reports. """
        return self.num_visible_het + self.num_visible_hom_alt

    @property
    def has_visible_observations(self) -> bool:
        return self.num_visible_alt > 0

    @property
    def has_observations(self) -> bool:
        return self.num_observations > 0

    @staticmethod
    def _get_sample_values_for_variant_via_cohort_genotype(locus_qs):
        """ This is the new, preferred way - as it gets all of the samples for a VCF at once.

            Archived VCFs naturally fall out: archive deletes the CohortGenotypeCollection rows
            (and partitions), so the cohortgenotype join yields no rows for archived VCFs. """
        no_cohort = Q(variant__cohortgenotype__isnull=True)
        vcf_cohort = Q(variant__cohortgenotype__collection__cohort__vcf__isnull=False)
        qs = locus_qs.filter(no_cohort | vcf_cohort)  # Only VCF CohortGenotypes (not generated cohorts)
        return qs.values("variant",
                         "variant__alt",
                         "variant__cohortgenotype__collection",
                         "variant__cohortgenotype__collection__cohort__vcf__name",
                         "variant__cohortgenotype__collection__cohort__vcf__genome_build",
                         "variant__cohortgenotype__filters",
                         "variant__cohortgenotype__samples_filters",
                         "variant__cohortgenotype__samples_zygosity",
                         "variant__cohortgenotype__samples_allele_frequency",
                         "variant__cohortgenotype__samples_allele_depth",
                         "variant__cohortgenotype__samples_read_depth",
                         "variant__cohortgenotype__samples_phred_likelihood")

    @staticmethod
    def _cohort_genotype_to_sample_genotypes(values_qs):
        """ We're now joining to CohortGenotype - break up into samples so old code works """

        filter_code_lookup_by_vcf = {}  # Filter codes are per-VCF, so can't share a lookup between them

        # turn a row of CohortGenotype into multiple rows of ObservedVariant values
        for cg_values in values_qs:
            variant = cg_values["variant"]
            cgc_id = cg_values["variant__cohortgenotype__collection"]
            if cgc_id is None:
                yield {
                    "variant": variant,
                }
                continue

            vcf_name = cg_values["variant__cohortgenotype__collection__cohort__vcf__name"]
            genome_build = cg_values["variant__cohortgenotype__collection__cohort__vcf__genome_build"]
            samples_zygosity = cg_values["variant__cohortgenotype__samples_zygosity"]
            num_samples = len(samples_zygosity)
            empty = [-1] * num_samples
            samples_allele_frequency = cg_values["variant__cohortgenotype__samples_allele_frequency"] or empty
            samples_allele_depth = cg_values["variant__cohortgenotype__samples_allele_depth"] or empty
            samples_read_depth = cg_values["variant__cohortgenotype__samples_read_depth"] or empty
            samples_phred_likelihood = cg_values["variant__cohortgenotype__samples_phred_likelihood"] or empty
            samples_filters = cg_values["variant__cohortgenotype__samples_filters"]

            cgc = CohortGenotypeCollection.objects.get(pk=cgc_id)
            if cgc.cohort_version != cgc.cohort.version:
                logging.warning("CohortGenotypeCollection (without VCF) %s v.%d out of date with Cohort v.%d",
                    cgc.pk, cgc.cohort_version, cgc.cohort.version)

                cohort_can_change = cgc.cohort.vcf is None
                if cohort_can_change:
                    continue  # Otherwise not really out of date, just a bug - see issue #1053

            vcf = cgc.cohort.vcf
            if (filter_code_lookup := filter_code_lookup_by_vcf.get(vcf.pk)) is None:
                filter_code_lookup = VCFFilter.get_code_lookup(vcf)
                filter_code_lookup_by_vcf[vcf.pk] = filter_code_lookup
            filters = VCFFilter.format_filter_codes(filter_code_lookup,
                                                    cg_values["variant__cohortgenotype__filters"])

            samples_qs = cgc.cohort.get_samples()

            ontology_path = f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__name"
            q_hpo = Q(**{f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.HPO})
            q_omim = Q(**{f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.OMIM})
            q_mondo = Q(**{f"patient__{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.MONDO})
            annotation_kwargs = {"patient_hpo": StringAgg(ontology_path, '|', filter=q_hpo,
                                                          distinct=True, output_field=TextField()),
                                 "patient_omim": StringAgg(ontology_path, '|', filter=q_omim,
                                                           distinct=True, output_field=TextField()),
                                 "patient_mondo": StringAgg(ontology_path, '|', filter=q_mondo,
                                                            distinct=True, output_field=TextField())}
            samples_qs = samples_qs.annotate(**annotation_kwargs)

            COPY_SAMPLE_FIELDS = ["id", "name", "patient", "patient__patient_code",
                                  SAMPLE_ENRICHMENT_KIT_PATH]
            sample_values = samples_qs.values("no_dna_control", "vcf__allele_frequency_percent",
                                              *COPY_SAMPLE_FIELDS, *list(annotation_kwargs.keys()))

            for s_values in sample_values:
                sample_id = s_values["id"]
                i = cgc.get_array_index_for_sample_id(sample_id)
                zygosity = samples_zygosity[i]
                if zygosity == Zygosity.MISSING or s_values["no_dna_control"]:
                    continue

                # The displayed AF may be percent formatted text, so keep the unit value for graphing
                source_in_percent = s_values["vcf__allele_frequency_percent"]
                allele_frequency_unit = format_af(samples_allele_frequency[i],
                                                  source_in_percent=source_in_percent,
                                                  dest_in_percent=False,
                                                  missing_value=CohortGenotype.MISSING_NUMBER_VALUE)
                allele_frequency = format_af(samples_allele_frequency[i],
                                             source_in_percent=source_in_percent,
                                             dest_in_percent=settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT,
                                             missing_value=CohortGenotype.MISSING_NUMBER_VALUE)
                allele_depth = samples_allele_depth[i]
                read_depth = samples_read_depth[i]
                phred_likelihood = samples_phred_likelihood[i]

                sample_filters = None  # VCF had no FT field
                if samples_filters:
                    if (ft := samples_filters[i]) != CohortGenotype.MISSING_FT_VALUE:
                        sample_filters = VCFFilter.format_filter_codes(filter_code_lookup, ft)

                sample_genotype_values = {
                    "variant": variant,
                    "genome_build": genome_build,
                    "zygosity": zygosity,
                    "allele_frequency": allele_frequency,
                    "allele_frequency_unit": allele_frequency_unit,
                    "allele_depth": allele_depth,
                    "read_depth": read_depth,
                    "phred_likelihood": phred_likelihood,
                    "filters": filters,
                    "sample_filters": sample_filters,
                    "sample": sample_id,
                    "sample__vcf__name": vcf_name,
                }
                for k, v in s_values.items():
                    if k in COPY_SAMPLE_FIELDS:
                        k = "sample__" + k
                    sample_genotype_values[k] = v
                yield sample_genotype_values


class VariantSampleGenotypes(VariantZygosityCounts):
    """ Builds the JSON payload for the variant details samples grid - one response per variant.
        The grid is client side, so everything it draws (rows, locus counts, multi-allelic
        variants, sample/observation summary) comes from here. """

    def to_json(self) -> dict:
        return {
            "genome_builds": [gb.name for gb in self.genome_builds],
            "variant_id": self.variant.pk,
            "variant": str(self.variant),
            "variant_url": reverse("view_variant", kwargs={"variant_id": self.variant.pk}),
            "g_hgvs": self._get_g_hgvs(),
            "samples": {
                "total": self.num_samples,
                "visible": self.num_user_samples,
            },
            "observations": {
                "total": self.num_observations,
                "visible": self.num_visible_observations,
                "invisible": self.num_invisible_observations,
            },
            "zygosity_counts": dict(self.visible_zygosity_counter),
            "locus_counts": self._get_locus_counts(),
            "multiallelic": self._get_multiallelic(),
            "rows": [self._row_to_json(row) for row in self.visible_rows],
        }

    def _get_g_hgvs(self) -> Optional[str]:
        """ Only used to name the CSV export - reference variants have to generate it, which needs
            annotation the deployment may not have, and that's not worth losing the samples over """
        try:
            return VariantAnnotation.get_hgvs_g(self.variant)
        except Exception:  # pylint: disable=broad-except
            log_traceback()
            return None

    @staticmethod
    def _row_to_json(row: dict) -> dict:
        sample_id = row["sample"]
        allele_frequency_unit = row["allele_frequency_unit"]
        if allele_frequency_unit == CohortGenotype.MISSING_NUMBER_VALUE:
            allele_frequency_unit = None
        return {
            "genome_build": row["genome_build"],
            "sample": sample_id,
            "sample_name": row["sample__name"],
            "sample_url": reverse("view_sample", kwargs={"sample_id": sample_id}),
            "patient": row["sample__patient"],  # Phenotype graphs count distinct patients, not samples
            "vcf": row["sample__vcf__name"],
            "zygosity": row["zygosity"],
            "allele_frequency": row["allele_frequency"],
            "allele_frequency_unit": allele_frequency_unit,  # Graphs need a number, not the display value
            "allele_depth": row["allele_depth"],
            "read_depth": row["read_depth"],
            "phred_likelihood": row["phred_likelihood"],
            "filters": row["filters"],
            "sample_filters": row["sample_filters"],
            "enrichment_kit": row["sample__" + SAMPLE_ENRICHMENT_KIT_PATH],
            "patient_hpo": row["patient_hpo"],
            "patient_omim": row["patient_omim"],
            "patient_mondo": row["patient_mondo"],
        }

    def _get_locus_counts(self) -> list[dict]:
        """ Zygosity counts for every variant at this locus, this variant first """
        variant_by_id = {v.pk: v for v in Variant.objects.filter(pk__in=self.locus_counter)}

        sorted_rows = []
        for variant_id, zygosity_counts in self.locus_counter.items():
            v = variant_by_id[variant_id]
            row = {
                "variant_id": variant_id,
                "variant": str(v),
                "url": reverse("view_variant", kwargs={"variant_id": variant_id}),
                "description": "",
                "total": sum(zygosity_counts.values()),
            }
            for zygosity, label in LOCUS_COUNT_ZYGOSITY_LABELS:
                row[label] = zygosity_counts[zygosity]

            if variant_id == self.variant.pk:
                row["description"] = "This variant"
                sort_order = 1
            elif v.is_reference:
                sort_order = 0
            else:
                sort_order = sum(map(ord, v.alt.seq))
            sorted_rows.append((sort_order, row))

        sorted_rows.sort(key=lambda sort_order_row: sort_order_row[0])
        return [row for _, row in sorted_rows]

    def _get_multiallelic(self) -> list[dict]:
        """ Variants that were split off this one's VCF record onto another locus """
        by_multiallelic = ModifiedImportedVariant.get_other_loci_variants_by_multiallelic(self.variant)
        multiallelic = []
        for old_multiallelic, variants in by_multiallelic.items():
            variant_list = [{"id": v.pk, "label": str(v), "url": v.get_absolute_url()}
                            for v in sorted(variants, key=lambda v: v.pk)]
            multiallelic.append({"multiallelic": old_multiallelic, "variants": variant_list})
        return multiallelic
