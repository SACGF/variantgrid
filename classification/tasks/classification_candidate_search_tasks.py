import itertools
import operator
from collections import defaultdict
from functools import reduce
from typing import Iterable

from django.db.models import Q

from analysis.models import Candidate
from analysis.tasks.abstract_candidate_search_task import AbstractCandidateSearchTask
from classification.enums import AlleleOriginBucket
from classification.models import Classification
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.views.classification_datatables import ClassificationColumns
from patients.models_enums import Zygosity
from snpdb.models import Sample, GenomeBuild
from snpdb.sample_filters import get_sample_ontology_q, get_sample_qc_gene_list_gene_symbol_q
from variantgrid.celery import app




class CrossSampleClassificationCandidateSearchTask(AbstractCandidateSearchTask):
    @staticmethod
    def _filter_classifications_by_sample_and_patient(sample: Sample, classifications: Iterable[Classification]) -> Iterable[Classification]:
        for c in classifications:
            if sample and c.sample:
                if sample == c.sample:
                    continue
                if sample.patient and sample.patient == c.sample.patient:
                    continue
            yield c

    @staticmethod
    def limit_sample_and_results(sample_records, max_samples, max_results):
        num_results = 0
        for sample, results in itertools.islice(sample_records, max_samples):
            max_remaining = max_results - num_results
            if max_remaining <= 0:
                break
            limited = results[:max_remaining]
            yield sample, limited
            num_results += len(limited)
            if num_results >= max_results:
                break

    @staticmethod
    def _get_sample_qs(user, config: dict):
        # Sample filters
        sample_filters = []
        if ontology_terms := config.get("sample_ontology_term_id"):
            if q := get_sample_ontology_q(ontology_terms):
                sample_filters.append(q)

        if gene_symbol_str := config.get("sample_gene_symbol"):
            if q := get_sample_qc_gene_list_gene_symbol_q(gene_symbol_str):
                sample_filters.append(q)

        sample_qs = Sample.filter_for_user(user)
        if sample_filters:
            sample_qs = sample_qs.filter(*sample_filters)
        return sample_qs

    @staticmethod
    def _get_classification_qs(user, config: dict):
        cs_data = {}
        for cs in ClassificationColumns.CLINICAL_SIGNIFICANCE_FILTERS:
            if config.get(cs):
                cs_data[cs] = True

        classification_filters = []
        if q := ClassificationColumns.get_clinical_significance_q(cs_data):
            classification_filters.append(q)

        if allele_origin := config.get("allele_origin"):
            if allele_origin != "A":
                classification_filters.append(Q(allele_origin_bucket__in=[allele_origin, AlleleOriginBucket.UNKNOWN]))

        if gene_symbol_str := config.get("gene_symbol"):
            if q := classification_gene_symbol_filter(gene_symbol_str):
                classification_filters.append(q)

        # request.GET.get("classification_id_filter")
        if lab_id := config.get("lab"):
            lab_list = lab_id.split(",")
            classification_filters.append(Q(lab__pk__in=lab_list))

        if ontology_terms := config.get("ontology_term_id"):
            terms = []
            for term_id in ontology_terms.split(","):
                terms.append(Q(condition_resolution__resolved_terms__contains=[{"term_id": term_id}]))
            if terms:
                q = reduce(operator.or_, terms)
                classification_filters.append(q)

        if config.get("internal_requires_sample"):
            classification_filters.append(Q(lab__external=True) | Q(sample__isnull=False))

        # TODO
        config.get("user")

        # Needs a variant to look in samples
        classification_qs = Classification.filter_for_user(user).filter(variant__isnull=False)
        if classification_filters:
            classification_qs = classification_qs.filter(*classification_filters)
        return classification_qs

    @staticmethod
    def _sample_classification_overlaps(samples_qs, classification_qs, zygosities):
        classifications_by_allele = defaultdict(set)
        for c in classification_qs:
            classifications_by_allele[c.variant.allele].add(c)

        for genome_build in GenomeBuild.builds_with_annotation():
            variant_q = Classification.get_variant_q_from_classification_qs(classification_qs, genome_build)
            # variant_qs = Variant.objects.filter(variant_q)
            # variants_qs = Variant.objects.filter(variantallele__allele__in=classifications_by_allele,
            #                                     variantallele__genome_build=genome_build)

            # Easier when allele is on Classification (like master)
            print(f"{genome_build=}")
            for sample in samples_qs.filter(vcf__genome_build=genome_build):
                sample_variant_zyg_and_classifications = []
                filter_kwargs = {}
                if zygosities:
                    filter_kwargs[f"{sample.zygosity_alias}__in"] = zygosities

                for v in sample.get_variant_qs().filter(variant_q, **filter_kwargs).distinct():
                    sample_zygosity = getattr(v, sample.zygosity_alias)

                    # existing_class = ret_classifications_path.filter(sample=s, variant__variantallele__allele__variantallele__variant=v=v)
                    # print(f"{s} has {v}")

                    classifications = classifications_by_allele[v.allele]
                    classifications = list(CrossSampleClassificationCandidateSearchTask._filter_classifications_by_sample_and_patient(sample, classifications))
                    if classifications:
                        sample_variant_zyg_and_classifications.append(
                            (v, sample_zygosity, classifications))

                if sample_variant_zyg_and_classifications:
                    yield sample, sample_variant_zyg_and_classifications

    def get_candidate_records(self, candidate_search_run):
        config = candidate_search_run.config_snapshot
        sample_qs = self._get_sample_qs(candidate_search_run.user, config)
        classification_qs = self._get_classification_qs(candidate_search_run.user, config)

        # Search
        search_max_results = int(config.get("max_results"))
        search_max_samples = int(config.get("max_samples"))
        ZYG_NAMES = {
            "hom_ref": Zygosity.HOM_REF,
            "het": Zygosity.HET,
            "hom_alt": Zygosity.HOM_ALT,
        }

        zygosities = []
        for zyg, code in ZYG_NAMES.items():
            if config.get(zyg):
                zygosities.append(code)

        records = []
        sample_records = self._sample_classification_overlaps(sample_qs, classification_qs, zygosities)
        for sample, data in self.limit_sample_and_results(sample_records, search_max_results, search_max_samples):
            for (v, sample_zygosity, classifications) in data:
                for classification in classifications:
                    notes = None
                    evidence = {}
                    records.append(Candidate(
                        search_run=candidate_search_run,
                        sample=sample,
                        variant=v,
                        classification=classification,
                        notes=notes,
                        evidence=evidence,
                        zygosity=sample_zygosity,
                    ))
        return records


class ClassificationEvidenceUpdateCandidateSearchTask(AbstractCandidateSearchTask):
    def get_candidate_records(self, candidate_search_run):
        # TODO
        return []

CrossSampleClassificationCandidateSearchTask = app.register_task(CrossSampleClassificationCandidateSearchTask())
ClassificationEvidenceUpdateCandidateSearchTask = app.register_task(ClassificationEvidenceUpdateCandidateSearchTask())

