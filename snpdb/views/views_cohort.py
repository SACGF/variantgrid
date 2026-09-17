import json
import logging

from django.core.exceptions import PermissionDenied
from django.http.response import (
    HttpResponseRedirect,
    JsonResponse,
)
from django.shortcuts import redirect, render
from django.urls.base import reverse
from django.views.decorators.http import require_POST

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
from snpdb.models import (
    Cohort,
    Duo,
    Quad,
    Sample,
    Trio,
    UserGridConfig,
    UserSettings,
)
from snpdb.models.models_enums import (
    ImportStatus,
    ProcessingStatus,
)
from snpdb.views.vcf_cohort_page import sample_membership_rows, vcf_cohort_page_context
from snpdb.views.views_json import cohort_genotype_json_response
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


def view_cohort(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.vcf:
        return redirect('view_vcf', vcf_id=cohort.vcf.pk)

    has_write_permission = cohort.can_write(request.user) and not cohort.data_archived
    cohort_form = forms.CohortForm(request.POST or None, instance=cohort)
    if request.method == "POST":
        if not has_write_permission:
            raise PermissionDenied()
        if valid := cohort_form.is_valid():
            cohort = cohort_form.save()
        add_save_message(request, valid, "Cohort")

    context = vcf_cohort_page_context(cohort, has_write_permission)
    context["cohort_form"] = cohort_form
    return render(request, 'snpdb/data/view_vcf_cohort.html', context)


def cohort_sample_rows(request, cohort_id):
    """ Membership editor rows for samples being staged - by 'vcf_id' (add all from a VCF) or 'sample_ids' """
    cohort = Cohort.get_for_user(request.user, cohort_id)
    sample_qs = Sample.filter_for_user(request.user).filter(vcf__genome_build=cohort.genome_build,
                                                            import_status=ImportStatus.SUCCESS)
    if vcf_id := request.GET.get("vcf_id"):
        sample_qs = sample_qs.filter(vcf_id=vcf_id)
    elif sample_ids := request.GET.get("sample_ids"):
        sample_qs = sample_qs.filter(pk__in=json.loads(sample_ids))
    else:
        sample_qs = sample_qs.none()

    samples = list(sample_qs.select_related("vcf").order_by("pk"))
    return JsonResponse({"samples": sample_membership_rows(samples)})


@require_POST
def cohort_sample_edit(request, cohort_id):
    """ Batch membership save - 'sample_ids' is the whole desired membership, in display order.
        Applies the diff in one version bump then builds the genotype data (a cohort that now lives
        inside a single VCF converts to a sub cohort instead of rebuilding) """
    cohort = Cohort.get_for_user(request.user, cohort_id, write=True)
    if cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived; cohort is read-only.")

    sample_ids = json.loads(request.POST["sample_ids"])
    ordered_sample_ids = []
    for sample_id in sample_ids:
        sample = Sample.get_for_user(request.user, sample_id)
        if sample.vcf.genome_build != cohort.genome_build:
            raise ValueError(f"{sample} is {sample.vcf.genome_build}, cohort is {cohort.genome_build}")
        ordered_sample_ids.append(sample.pk)

    cohort.set_samples(ordered_sample_ids)
    logging.info("Cohort %s membership set to %d samples", cohort.pk, len(ordered_sample_ids))
    return cohort_genotype_json_response(cohort)


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


def duos(request):
    show_group_data = UserGridConfig.get(request.user, 'Duos').show_group_data
    context = {"show_group_data": show_group_data}
    return render(request, 'snpdb/patients/duos.html', context)


def view_duo(request, pk):
    duo = Duo.get_for_user(request.user, pk)
    context = {"duo": duo,
               "has_write_permission": duo.cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_duo.html', context)
