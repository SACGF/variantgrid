import json
import logging

from django.core.exceptions import PermissionDenied
from django.http.response import (
    HttpResponse,
    HttpResponseRedirect,
)
from django.shortcuts import redirect, render
from django.urls.base import reverse

from annotation.forms import GeneCountTypeChoiceForm
from annotation.models import (
    AnnotationVersion,
)
from annotation.models.models import VariantAnnotationVersion
from annotation.models.models_gene_counts import (
    CohortGeneCounts,
    GeneCountType,
)
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.forms import CustomGeneListForm, GeneAndTranscriptForm, UserGeneListForm
from genes.models import CustomTextGeneList, GeneList, GeneListCategory
from library.django_utils import (
    add_save_message,
)
from snpdb import forms
from snpdb.archive import (
    DataArchivedError,
)
from snpdb.forms import (
    SampleChoiceForm,
)
from snpdb.models import (
    Cohort,
    CohortGenotypeCollection,
    CohortSample,
    Quad,
    Trio,
    UserGridConfig,
    UserSettings,
)
from snpdb.models.models_enums import (
    ProcessingStatus,
)
from snpdb.views.views_sample_gene_matrix import sample_gene_matrix


def cohorts(request):
    user_settings = UserSettings.get_for_user(request.user)
    initial = {'genome_build': user_settings.default_genome_build}
    form = forms.CreateCohortForm(request.POST or None, initial=initial)
    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            cohort = form.save(commit=False)
            cohort.user = request.user
            cohort.save()
            return HttpResponseRedirect(reverse('view_cohort', kwargs={'cohort_id': cohort.pk}))
        else:
            add_save_message(request, valid, "Cohort", created=True)

    show_group_data = UserGridConfig.get(request.user, 'Cohorts').show_group_data
    context = {"form": form, "show_group_data": show_group_data}
    return render(request, 'snpdb/patients/cohorts.html', context)


def view_cohort_details_tab(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    has_write_permission = cohort.can_write(request.user) and not cohort.data_archived
    context = {"cohort": cohort,
               "has_write_permission": has_write_permission}
    return render(request, 'snpdb/patients/view_cohort_details_tab.html', context)


def view_cohort(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.vcf:
        return redirect('view_vcf', vcf_id=cohort.vcf.pk)

    try:
        cohort_genotype_collection = cohort.cohort_genotype_collection
    except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        cohort_genotype_collection = None

    has_write_permission = cohort.can_write(request.user) and not cohort.data_archived

    form = forms.CohortForm(request.POST or None, instance=cohort)
    if request.method == "POST":
        if not has_write_permission:
            raise PermissionDenied()
        if valid := form.is_valid():
            cohort = form.save()
        add_save_message(request, valid, "Cohort")

    sample_form = SampleChoiceForm(genome_build=cohort.genome_build)
    sample_form.fields['sample'].required = False

    context = {"form": form,
               "sample_form": sample_form,
               "cohort": cohort,
               "cohort_genotype_collection": cohort_genotype_collection,
               "has_write_permission": has_write_permission}
    return render(request, 'snpdb/patients/view_cohort.html', context)


def cohort_sample_edit(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived; cohort is read-only.")

    if request.method == "POST":
        cohort_op = request.POST['cohort_op']
        sample_ids_str = request.POST['sample_ids']
        sample_ids = json.loads(sample_ids_str)
        if cohort_op == 'add':
            for sample_id in sample_ids:
                cohort.add_sample(sample_id)
        elif cohort_op == 'remove':
            for sample_id in sample_ids:
                try:
                    cohort_sample = CohortSample.objects.get(cohort=cohort, sample_id=sample_id)
                    cohort_sample.delete()
                    logging.info("Removed: %s", sample_id)
                except CohortSample.DoesNotExist:
                    pass
        else:
            raise ValueError(f"Unknown cohort_op '{cohort_op}'")

    return HttpResponse()


def cohort_hotspot(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived; cohort is read-only.")
    vav = VariantAnnotationVersion.latest(cohort.genome_build)
    form = GeneAndTranscriptForm(gene_annotation_release=vav.gene_annotation_release,
                                 has_protein_domains=True)

    try:
        cohort_genotype_collection = cohort.cohort_genotype_collection
    except Exception as e:
        cohort_genotype_collection = None
        logging.error(e)

    context = {
        "cohort": cohort,
        "cohort_genotype_collection": cohort_genotype_collection,
        "form": form,
        "gene_annotation_release": vav.gene_annotation_release,
    }
    return render(request, 'snpdb/patients/cohort_hotspot.html', context)


def cohort_gene_counts(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived; cohort is read-only.")

    COHORT_CUSTOM_GENE_LIST = f"__QC_COVERAGE_CUSTOM_GENE_LIST__{request.user}"

    # We only want to keep 1 per user
    custom_text_gene_list, _ = CustomTextGeneList.objects.get_or_create(name=COHORT_CUSTOM_GENE_LIST)

    custom_gene_list_form = CustomGeneListForm(request.POST or None,
                                               initial={"custom_gene_list_text": custom_text_gene_list.text})
    if custom_gene_list_form.is_valid():
        custom_text_gene_list.text = custom_gene_list_form.cleaned_data['custom_gene_list_text']
        custom_text_gene_list.save()
        create_custom_text_gene_list(custom_text_gene_list, request.user, GeneListCategory.QC_COVERAGE_CUSTOM_TEXT,
                                     hidden=True)
        gene_list_id = custom_text_gene_list.gene_list.pk
    else:
        gene_list_id = None

    context = {"cohort": cohort,
               'gene_list_id': gene_list_id,
               'gene_list_form': UserGeneListForm(),
               'custom_gene_list_form': custom_gene_list_form,
               'has_gene_count_types': GeneCountType.objects.filter(enabled=True).exists(),
               'gene_count_type_choice_form': GeneCountTypeChoiceForm()}
    return render(request, 'snpdb/patients/cohort_gene_counts.html', context)


def cohort_gene_counts_matrix(request, cohort_id, gene_count_type_id, gene_list_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived; cohort is read-only.")
    gene_count_type = GeneCountType.objects.get(pk=gene_count_type_id)
    gene_list = GeneList.get_for_user(request.user, gene_list_id)
    samples = list(cohort.get_samples())

    annotation_version = AnnotationVersion.latest(cohort.genome_build)
    variant_annotation_version = annotation_version.variant_annotation_version
    cgc, created = CohortGeneCounts.objects.get_or_create(variant_annotation_version=variant_annotation_version,
                                                          gene_count_type=gene_count_type,
                                                          cohort=cohort,
                                                          cohort_version=cohort.version)

    graph_kwargs = {"cohort_id": cohort_id,
                    "gene_count_type_id": gene_count_type_id,
                    "gene_list_id": gene_list_id}
    redirect_url = reverse("cohort_gene_counts_matrix", kwargs=graph_kwargs)
    if created or (cgc.processing_status not in ProcessingStatus.FINISHED_STATES):
        celery_task = cgc.launch_task()
        wait_for_task_kwargs = {"celery_task": celery_task, "sleep_ms": 2000, "redirect_url": redirect_url}
        wait_url = reverse("wait_for_task", kwargs=wait_for_task_kwargs)
        return HttpResponseRedirect(wait_url)
    else:
        if cgc.processing_status == ProcessingStatus.SUCCESS:
            return sample_gene_matrix(request, variant_annotation_version, samples, gene_list, gene_count_type)
        else:
            raise ValueError(f"{cgc} had ProcessingStatus: {cgc.processing_status}")


def trios(request):
    show_group_data = UserGridConfig.get(request.user, 'Trios').show_group_data
    context = {"show_group_data": show_group_data}
    return render(request, 'snpdb/patients/trios.html', context)


def view_trio(request, pk):
    trio = Trio.get_for_user(request.user, pk)
    context = {"trio": trio,
               "has_write_permission": trio.cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_trio.html', context)


def quads(request):
    show_group_data = UserGridConfig.get(request.user, 'Quads').show_group_data
    context = {"show_group_data": show_group_data}
    return render(request, 'snpdb/patients/quads.html', context)


def view_quad(request, pk):
    quad = Quad.get_for_user(request.user, pk)
    context = {"quad": quad,
               "has_write_permission": quad.cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_quad.html', context)


def cohort_sort(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived; cohort is read-only.")
    if request.method == "POST":
        cohort_samples_str = request.POST.get("cohort_samples")
        cohort_samples_ids = cohort_samples_str.split(',') if cohort_samples_str else []
        cohort_samples = []
        for i, cs_id in enumerate(cohort_samples_ids):
            cohort_sample = CohortSample.objects.get(pk=cs_id, cohort=cohort)
            cohort_sample.sort_order = i
            cohort_sample.save()
            cohort_samples.append(cohort_sample)
    else:
        cohort_samples = cohort.get_cohort_samples()

    context = {'cohort': cohort,
               'cohort_samples': cohort_samples}
    return render(request, 'snpdb/patients/cohort_sort.html', context)
