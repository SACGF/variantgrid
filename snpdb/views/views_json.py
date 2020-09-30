from celery.result import AsyncResult
from django.conf import settings
from django.http.response import JsonResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_POST
import json
import os

from snpdb.models import CachedGeneratedFile, Cohort, Sample
from snpdb.tasks.cohort_genotype_tasks import create_cohort_genotype_and_launch_task


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
