import itertools
import operator
from collections import defaultdict
from functools import reduce
from typing import Iterable

from crispy_forms.layout import Layout, Field

from django.conf import settings
from django.db.models.query_utils import Q
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.shortcuts import render

from classification.enums import AlleleOriginBucket
from classification.forms import ClassificationAlleleOriginForm, SampleClassificationForm, ClinicalSignificanceForm
from classification.models import Classification
from classification.models.classification_utils import classification_gene_symbol_filter
from classification.views.classification_datatables import ClassificationColumns
from genes.forms import GeneSymbolForm
from genes.models import SampleGeneList
from ontology.forms import PhenotypeMultipleSelectForm
from patients.models_enums import Zygosity
from snpdb.forms import UserSelectForm, LabSelectForm, LabMultiSelectForm
from snpdb.models import Lab, Sample, GenomeBuild
from snpdb.sample_filters import get_sample_ontology_q, get_sample_qc_gene_list_gene_symbol_q
from snpdb.user_settings_manager import UserSettingsManager


def sample_classification_search(request) -> HttpResponse:
    # A lot of this code is copy/pasted from classifications - classification listing page
    # TODO: We should go and refactor this into a tag later

    # is cached on the request
    user_settings = UserSettingsManager.get_user_settings()

    if settings.CLASSIFICATION_GRID_MULTI_LAB_FILTER:
        lab_form = LabMultiSelectForm()
    else:
        lab_form = LabSelectForm()

    sample_gene_form = GeneSymbolForm(prefix="sample")
    sample_gene_form.fields['gene_symbol'].label = "Sample Gene List Gene Symbol"

    layout_fields = []
    for field_name in ["omim", "hpo", "mondo"]:
        layout_fields.append(Field(field_name, wrapper_class=field_name))
    layout = Layout(*layout_fields)

    sample_phenotype_form = PhenotypeMultipleSelectForm(prefix="sample")
    sample_phenotype_form.helper.layout = layout

    classification_phenotype_form = PhenotypeMultipleSelectForm(prefix="classification")
    classification_phenotype_form.helper.layout = layout

    cs_form = ClinicalSignificanceForm(initial={
        "likely_pathogenic": True,
        "pathogenic": True,
    })

    context = {
        "sample_phenotype_form": sample_phenotype_form,
        "gene_form": GeneSymbolForm(prefix="classification"),
        "user_form": UserSelectForm(),
        "lab_form": lab_form,
        "allele_origin_form": ClassificationAlleleOriginForm(),
        "labs": Lab.valid_labs_qs(request.user),
        "clinical_significance_form": cs_form,
        "classification_phenotype_form": classification_phenotype_form,
        "sample_classification_form": SampleClassificationForm(),
        "user_settings": user_settings,
    }

    samples_qs = Sample.filter_for_user(request.user)
    context["num_visible_samples"] = samples_qs.count()
    num_sample_gene_lists = SampleGeneList.objects.filter(sample__in=samples_qs).count()
    context["num_sample_gene_lists"] = num_sample_gene_lists
    if num_sample_gene_lists:
        context["sample_gene_form"] = sample_gene_form
    return render(request, 'classification/sample_classification_search.html', context)


def _filter_classifications_by_sample_and_patient(sample: Sample, classifications: Iterable[Classification]) -> Iterable[Classification]:
    for c in classifications:
        if sample and c.sample:
            if sample == c.sample:
                continue
            if sample.patient and sample.patient == c.sample.patient:
                continue
        yield c


def _sample_classification_overlaps(samples_qs, classification_qs, zygosities):
    classifications_by_allele = defaultdict(set)
    for c in classification_qs:
        classifications_by_allele[c.variant.allele].add(c)

    for genome_build in GenomeBuild.builds_with_annotation():
        variant_q = Classification.get_variant_q_from_classification_qs(classification_qs, genome_build)
        # variant_qs = Variant.objects.filter(variant_q)
        #variants_qs = Variant.objects.filter(variantallele__allele__in=classifications_by_allele,
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
                classifications = list(_filter_classifications_by_sample_and_patient(sample, classifications))
                if classifications:
                    sample_variant_zyg_and_classifications.append((v, Zygosity.display(sample_zygosity), classifications))

                # for c in classifications_by_allele[v.allele]:
                #     sample_msg = ""
                #     if c.sample == s:
                #         sample_msg = "This sample"
                #     else:
                #         if c.sample:
                #             sample_msg = "Class has other sample {c.sample}"
                #             if c.sample.patient:
                #                 if c.sample.patient == s.patient:
                #                     sample_msg = "Same patient, diff sample"
                #                 else:
                #                     print(f"diff patient: {c.sample.patient} vs {s.patient}")
                #         else:
                #             sample_id = c.evidence.get("sample_id", {}).get("value")
                #             sample_msg = f"Class has no sample - but {sample_id=}"
                #
                #     sample_alleles[v.allele] = (c, s, sample_msg)

            if sample_variant_zyg_and_classifications:
                yield sample, sample_variant_zyg_and_classifications


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




def sample_classification_search_results(request: HttpRequest) -> HttpResponse:
    # Sample filters
    sample_filters = []
    if ontology_terms := request.GET.get("sample_ontology_term_id"):
        if q := get_sample_ontology_q(ontology_terms):
            sample_filters.append(q)

    if gene_symbol_str := request.GET.get("sample_gene_symbol"):
        if q := get_sample_qc_gene_list_gene_symbol_q(gene_symbol_str):
            sample_filters.append(q)

    sample_qs = Sample.filter_for_user(request.user)
    num_unfiltered_samples = sample_qs.count()
    if sample_filters:
        sample_qs = sample_qs.filter(*sample_filters)

    # Classification filters
    classification_filters = []

    cs_data = {}
    for cs in ClassificationColumns.CLINICAL_SIGNIFICANCE_FILTERS:
        if request.GET.get(f"classification_{cs}"):
            cs_data[cs] = True

    if q := ClassificationColumns.get_clinical_significance_q(cs_data):
        classification_filters.append(q)

    if allele_origin := request.GET.get("classification_allele_origin"):
        if allele_origin != "A":
            classification_filters.append(Q(allele_origin_bucket__in=[allele_origin, AlleleOriginBucket.UNKNOWN]))

    if gene_symbol_str := request.GET.get("classification_gene_symbol"):
        if q := classification_gene_symbol_filter(gene_symbol_str):
            classification_filters.append(q)

    # request.GET.get("classification_id_filter")
    if lab_id := request.GET.get("classification_lab"):
        lab_list = lab_id.split(",")
        classification_filters.append(Q(lab__pk__in=lab_list))

    if ontology_terms := request.GET.get("classification_ontology_term_id"):
        terms = []
        for term_id in ontology_terms.split(","):
            terms.append(Q(condition_resolution__resolved_terms__contains=[{"term_id": term_id}]))
        if terms:
            q = reduce(operator.or_, terms)
            classification_filters.append(q)

    if request.GET.get("classification_internal_requires_sample"):
        classification_filters.append(Q(lab__external=True) | Q(sample__isnull=False))

    # Needs a variant to look in samples
    classification_qs = Classification.filter_for_user(request.user).filter(variant__isnull=False)
    num_unfiltered_classifications = classification_qs.count()
    if classification_filters:
        classification_qs = classification_qs.filter(*classification_filters)

    request.GET.get("classification_user")
    # Search
    search_max_results = int(request.GET.get("search_max_results"))
    search_max_samples = int(request.GET.get("search_max_samples"))
    ZYG_NAMES = {
        "hom_ref": Zygosity.HOM_REF,
        "het": Zygosity.HET,
        "hom_alt": Zygosity.HOM_ALT,
    }

    zygosities = []
    for zyg, code in ZYG_NAMES.items():
        if request.GET.get(f"search_{zyg}"):
            zygosities.append(code)

    num_samples = sample_qs.count()
    print(f"filtered - kept {num_samples=} of {num_unfiltered_samples} samples")
    num_classifications = classification_qs.count()
    print(f"filtered - kept {num_classifications=} of {num_unfiltered_classifications} classifications")

    sample_records = _sample_classification_overlaps(sample_qs, classification_qs, zygosities)
    context = {
        "search_max_samples": search_max_samples,
        "search_max_results": search_max_results,
        "sample_records": limit_sample_and_results(sample_records, search_max_samples, search_max_results),
    }
    return render(request, 'classification/sample_classification_search_results.html', context)
