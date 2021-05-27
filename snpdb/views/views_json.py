from celery.result import AsyncResult
from django.conf import settings
from django.http.response import JsonResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_POST
import json
import os

from library.django_utils import require_superuser
from snpdb.models import CachedGeneratedFile, Cohort, Sample, VCF, VCFAlleleSource, CustomColumnsCollection
from snpdb.tasks.cohort_genotype_tasks import create_cohort_genotype_and_launch_task
from snpdb.tasks.clingen_tasks import populate_clingen_alleles_from_allele_source
from snpdb.tasks.vcf_zygosity_count_tasks import update_variant_zygosity_count_for_vcf_task


def job_status(request, job_id):
    async_result = AsyncResult(job_id)
    data = {"status": async_result.status}
    if async_result.successful():
        data['result'] = async_result.result
    else:
        data['exception'] = str(async_result.result)

    return JsonResponse(data)


def cached_generated_file_check(request, cgf_id):
    cgf = get_object_or_404(CachedGeneratedFile, pk=cgf_id)
    cgf.check_ready()
    data = {"status": cgf.task_status,
            "cgf_id": cgf.id}

    if cgf.exception:
        data["exception"] = str(cgf.exception)
    elif cgf.task_status == "SUCCESS":
        media_root_dir = os.path.join(settings.MEDIA_ROOT, "")  # with end slash
        file_path = cgf.filename.replace(media_root_dir, "")
        data["url"] = os.path.join(settings.MEDIA_URL, file_path)

    return JsonResponse(data)


def create_cohort_genotype(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    status, celery_task = create_cohort_genotype_and_launch_task(cohort)

    data = {}
    if status:
        data["status"] = status
    elif celery_task:
        data["celery_task"] = celery_task
    return JsonResponse(data)


@require_POST
def create_sub_cohort(request, cohort_id):
    parent_cohort = Cohort.get_for_user(request.user, cohort_id)
    sample_id_list_str = request.POST["sample_id_list"]
    sample_id_list = json.loads(sample_id_list_str)
    sample_list = []
    for sample_id in sample_id_list:
        sample_list.append(Sample.get_for_user(request.user, sample_id))

    sub_cohort = parent_cohort.create_sub_cohort(request.user, sample_list)
    return JsonResponse({"cohort_id": sub_cohort.pk})


def cohort_sample_count(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    return JsonResponse(cohort.cohortsample_set.count(), safe=False)


@require_POST
@require_superuser
def vcf_populate_clingen_alleles(request, vcf_id):
    vcf = VCF.get_for_user(request.user, vcf_id)
    vcf_as, _ = VCFAlleleSource.objects.get_or_create(vcf=vcf)
    populate_clingen_alleles_from_allele_source.si(vcf_as.pk, settings.CLINGEN_ALLELE_REGISTRY_MAX_MANUAL_REQUESTS).apply_async()
    return JsonResponse({})


@require_POST
def vcf_change_zygosity_count(request, vcf_id, vzcc_id, operation):
    vcf = VCF.get_for_user(request.user, vcf_id)
    if not vcf.can_write(request.user):
        raise PermissionError(f"You do not have write permission on VCF pk={vcf_id}")

    OPERATION = {"add": "+", "del": "-"}
    operation = OPERATION[operation]
    task = update_variant_zygosity_count_for_vcf_task.si(vzcc_id, vcf_id, operation)
    task.apply_async()
    return JsonResponse({})


@require_POST
def clone_custom_columns(request, custom_columns_collection_id):
    ccc = CustomColumnsCollection.get_for_user(request.user, custom_columns_collection_id)
    cloned_ccc = ccc.clone_for_user(request.user)
    return JsonResponse({"pk": cloned_ccc.pk})
