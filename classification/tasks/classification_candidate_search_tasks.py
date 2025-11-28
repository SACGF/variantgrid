import operator
from collections import defaultdict
from functools import reduce
from typing import Iterable

from django.db.models import Q, QuerySet

from analysis.models import Candidate
from analysis.tasks.abstract_candidate_search_task import AbstractCandidateSearchTask
from annotation.models import AnnotationVersion, VariantAnnotation, ClinVar, ClinVarReviewStatus
from classification.enums import AlleleOriginBucket, ClinicalSignificance
from classification.models import Classification, ClassificationModification, EvidenceKey
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.views.classification_datatables import ClassificationColumns
from snpdb.models import Sample, GenomeBuild
from snpdb.sample_filters import get_sample_ontology_q, get_sample_qc_gene_list_gene_symbol_q
from variantgrid.celery import app


class ClassificationCandidateSearchMixin:
    @staticmethod
    def _get_classification_modifications_qs(user, config: dict) -> QuerySet[ClassificationModification]:
        cs_data = {}
        for cs in ClassificationColumns.CLINICAL_SIGNIFICANCE_FILTERS:
            if config.get(cs):
                cs_data[cs] = True

        cm_filters = []
        if q := ClassificationColumns.get_clinical_significance_q(cs_data):
            cm_filters.append(q)

        if allele_origin := config.get("allele_origin"):
            if allele_origin != "A":
                cm_filters.append(Q(classification__allele_origin_bucket__in=[allele_origin, AlleleOriginBucket.UNKNOWN]))

        if gene_symbol_str := config.get("gene_symbol"):
            if q := classification_gene_symbol_filter(gene_symbol_str):
                cm_filters.append(q)

        # request.GET.get("classification_id_filter")
        if lab_id := config.get("lab"):
            lab_list = lab_id.split(",")
            cm_filters.append(Q(classification__lab__pk__in=lab_list))

        if ontology_terms := config.get("ontology_term_id"):
            terms = []
            for term_id in ontology_terms.split(","):
                terms.append(Q(classification__condition_resolution__resolved_terms__contains=[{"term_id": term_id}]))
            if terms:
                q = reduce(operator.or_, terms)
                cm_filters.append(q)

        if config.get("internal_requires_sample"):
            cm_filters.append(Q(classification__lab__external=True) | Q(classification__sample__isnull=False))

        # TODO
        config.get("user")

        # Classifications must be local (not external) and matched to a variant
        cm_qs = ClassificationModification.latest_for_user(
            user=user,
            published=True).filter(classification__lab__external=False, classification__variant__isnull=False)
        if cm_filters:
            cm_qs = cm_qs.filter(*cm_filters)
        return cm_qs



class CrossSampleClassificationCandidateSearchTask(ClassificationCandidateSearchMixin, AbstractCandidateSearchTask):
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
    def _get_sample_qs(user, config: dict):
        # Sample filters - TODO: This is basically a copy of SampleColumns filtering - extract to common code?
        sample_filters = []
        ontology_filters = []
        for ontology_service in ["hpo", "omim", "mondo"]:
            if ontology_terms := config.get(f"sample-{ontology_service}"):
                if q := get_sample_ontology_q(ontology_terms):
                    ontology_filters.append(q)
        if ontology_filters:
            # Do an or (meaning any of the phenos)
            q = reduce(operator.or_, ontology_filters)
            sample_filters.append(q)

        if gene_symbol_str := config.get("sample_gene_symbol"):
            if q := get_sample_qc_gene_list_gene_symbol_q(gene_symbol_str):
                sample_filters.append(q)

        if project := config.get("project"):
            sample_filters.append(Q(vcf__project=project))

        if vcf := config.get("vcf"):
            sample_filters.append(Q(vcf=vcf))

        sample_qs = Sample.filter_for_user(user)
        if sample_filters:
            sample_qs = sample_qs.filter(*sample_filters)
        return sample_qs

    @staticmethod
    def _sample_classification_overlaps(samples_qs, cm_qs, zygosities):
        classifications_by_allele = defaultdict(set)
        for cm in cm_qs:
            classifications_by_allele[cm.classification.variant.allele].add(cm.classification)

        classification_qs = Classification.objects.filter(classificationmodification__in=cm_qs)

        for genome_build in GenomeBuild.builds_with_annotation():
            variant_q = Classification.get_variant_q_from_classification_qs(classification_qs, genome_build)
            # variant_qs = Variant.objects.filter(variant_q)
            # variants_qs = Variant.objects.filter(variantallele__allele__in=classifications_by_allele,
            #                                     variantallele__genome_build=genome_build)

            # Easier when allele is on Classification (like master)
            print(f"{genome_build=}")
            for sample in samples_qs.filter(vcf__genome_build=genome_build):
                print(f"Searching {sample}")
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
        cm_qs = self._get_classification_modifications_qs(candidate_search_run.user, config)

        # Search
        search_max_results = int(config.get("max_results"))
        search_max_samples = int(config.get("max_samples"))
        zygosities = candidate_search_run.get_zygosities_from_config()

        records = []
        sample_records = self._sample_classification_overlaps(sample_qs, cm_qs, zygosities)
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


class ClassificationEvidenceUpdateCandidateSearchTask(ClassificationCandidateSearchMixin, AbstractCandidateSearchTask):
    @staticmethod
    def maybe_splicing_related(evidence):
        # TODO: we are also trying to do similar "look for splicing" in rna4rd so perhaps split this into common library code
        if molecular_consequence := EvidenceKey.get_value(evidence.get("molecular_consequence")):
            for mc in molecular_consequence:
                mc = mc.lower()
                if "splice" in mc or "intron" in mc:
                    return True
        if EvidenceKey.get_value(evidence.get("intron")):
            return True
        return False

    def get_candidate_records(self, candidate_search_run):
        config = candidate_search_run.config_snapshot
        records = []

        check_population = config.get("population")
        check_clinvar = config.get("clinvar")
        check_computational = config.get("computational")
        check_gene_disease = config.get("gene_disease")

        def none_or_cast(v, op):
            if v is not None:
                return op(v)

        pop_no_ba1_min_af = none_or_cast(config.get("pop_no_ba1_min_af"), float)
        pop_no_bs1_min_af = none_or_cast(config.get("pop_no_bs1_min_af"), float)
        pop_recessive_no_bs2_min_homozygotes = none_or_cast(config.get("pop_recessive_no_bs2_min_homozygotes"), int)
        pop_pm2_max_af = none_or_cast(config.get("pop_pm2_max_af"), float)
        clinvar_min_conflict_distance = none_or_cast(config.get("clinvar_min_conflict_distance"), int)
        clinvar_min_stars = none_or_cast(config.get("clinvar_min_stars"), int)
        computational_vus_spliceai_min = none_or_cast(config.get("computational_vus_spliceai_min"), float)

        # We need to use genome builds because we're going to pull in the annotations
        cm_qs = self._get_classification_modifications_qs(candidate_search_run.user, config)

        for genome_build in GenomeBuild.builds_with_annotation():
            av = AnnotationVersion.latest(genome_build)

            cm_build_qs = cm_qs.filter(classification__allele_info__imported_genome_build_patch_version__genome_build=genome_build)

            # Read in annotation / ClinVar per variant
            variant_ids = list(cm_build_qs.values_list("classification__variant_id", flat=True).distinct())
            va_by_variant_id = {}
            if check_population or check_computational:
                va_qs = VariantAnnotation.objects.filter(version=av.variant_annotation_version, variant__in=variant_ids)
                va_by_variant_id = {va.variant_id: va for va in va_qs}

            clinvar_by_variant_id = {}
            if check_clinvar:
                clinvar_qs = ClinVar.objects.filter(version=av.clinvar_version, variant__in=variant_ids)
                if clinvar_min_stars is not None:
                    review_statuses = ClinVarReviewStatus.statuses_gte_stars(clinvar_min_stars)
                    clinvar_qs = clinvar_qs.filter(review_status__in=review_statuses)

                clinvar_by_variant_id = {cv.variant_id: cv for cv in clinvar_qs}

            for cm in cm_build_qs:
                classification = cm.classification
                variant_id = classification.variant_id
                notes = None
                candidate_evidence = {}
                # clinvar = None
                va = va_by_variant_id.get(variant_id)
                cv = clinvar_by_variant_id.get(variant_id)

                def evidence(key):
                    return EvidenceKey.get_value(cm.published_evidence.get(key))

                if check_population and va:
                    if cm.clinical_significance in (ClinicalSignificance.VUS, ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.PATHOGENIC):
                        popmax_af = va.gnomad_popmax_af
                        if popmax_af is not None:
                            if pop_no_ba1_min_af is not None:
                                if not evidence("acmg:ba1") and popmax_af >= pop_no_ba1_min_af:
                                    candidate_evidence["acmg:ba1 missing"] = f"{popmax_af=} >= {pop_no_ba1_min_af}"

                            if pop_no_bs1_min_af is not None:
                                if not evidence("acmg:bs1") and popmax_af >= pop_no_bs1_min_af:
                                    candidate_evidence["acmg:bs1 missing"] = f"{popmax_af=} >= {pop_no_bs1_min_af}"

                            if pop_pm2_max_af is not None:
                                if evidence("acmg:pm2") and popmax_af >= pop_pm2_max_af:
                                    candidate_evidence["acmg:pm2"] = f"PM2 set with {popmax_af=} >= {pop_pm2_max_af}"

                        if pop_recessive_no_bs2_min_homozygotes is not None:
                            if not evidence("acmg:bs2"):
                                if evidence("mode_of_inheritance") == "autosomal_recessive":
                                    if gnomad_num_homozygotes := va.gnomad_hom_alt:
                                        if gnomad_num_homozygotes >= pop_recessive_no_bs2_min_homozygotes:
                                            candidate_evidence["acmg:bs2"] = f"{gnomad_num_homozygotes=} >= {pop_recessive_no_bs2_min_homozygotes}"


                if check_clinvar and cv:
                    if clinvar_min_conflict_distance is not None:
                        distance = ClinicalSignificance.distance(cm.clinical_significance, cv.highest_pathogenicity)
                        if abs(distance) >= clinvar_min_conflict_distance:
                            json_summary = cv.json_summary()
                            json_summary["distance"] = distance
                            candidate_evidence["clinvar"] = json_summary

                if check_computational and va:

                    if computational_vus_spliceai_min is not None:
                        if cm.clinical_significance in (ClinicalSignificance.VUS,
                                                        ClinicalSignificance.LIKELY_PATHOGENIC):

                            if not evidence("spliceai"):  # Missing
                                splice_flag_reasons = []
                                # If pp3 is flagged nothing to do
                                pp3 = evidence("acmg:pp3")
                                if pp3 is None:
                                    splice_flag_reasons.append("PP3 not applied")
                                    splicing_assertion = evidence("splicing_assertion")
                                    bp4 = evidence("acmg:bp4")
                                    if bp4 is not None:
                                        splice_flag_reasons.append("BP4 applied")
                                    if splicing_assertion is not None:
                                        if splicing_assertion in ("no_effect", "No effect"):
                                            splice_flag_reasons.append(f"{splicing_assertion=}")
                                    else:
                                        splice_flag_reasons.append("No splicing assertion applied")

                                if splice_flag_reasons and self.maybe_splicing_related(cm.published_evidence):
                                    if highest_spliceai := va.highest_spliceai():
                                        if highest_spliceai >= computational_vus_spliceai_min:
                                            reason = ", ".join(splice_flag_reasons)
                                            spliceai_evidence = f"{reason}, {highest_spliceai=} >= {computational_vus_spliceai_min}"
                                            candidate_evidence["splicing"] = spliceai_evidence

                if check_gene_disease:
                    # EKey: gene_disease_validity
                    # see MOINode.get_gene_disease_relations
                    # moderate : Moderate
                    # strong : Strong
                    # definitive : Definitive
                    # TODO: Need to read in Gene disease associations etc
                    pass

                if notes or candidate_evidence:
                    records.append(Candidate(
                        search_run=candidate_search_run,
                        variant=classification.variant,
                        classification=classification,
                        annotation_version=av,
                        notes=notes,
                        evidence=candidate_evidence,
                        # clinvar=clinvar,
                    ))
        return records


CrossSampleClassificationCandidateSearchTask = app.register_task(CrossSampleClassificationCandidateSearchTask())
ClassificationEvidenceUpdateCandidateSearchTask = app.register_task(ClassificationEvidenceUpdateCandidateSearchTask())

