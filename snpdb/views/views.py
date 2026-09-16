import json

from celery.result import AsyncResult
from django.http.response import (
    HttpResponse,
    HttpResponseRedirect,
    JsonResponse,
)
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse
from django.views.decorators.http import require_POST
from global_login_required import login_not_required
from termsandconditions.decorators import terms_required

from library import uptime_check
from snpdb.models import (
    CachedGeneratedFile,
    Contig,
    GenomeBuild,
    UserSettings,
    Wiki,
)
from user_messages.models import Message


@terms_required
def index(request):
    return render(request, 'index.html')


def maps(request):
    return render(request, 'maps.html')


@require_POST
def cached_generated_file_delete(request):
    cgf_id = request.POST["cgf_id"]
    cgf = get_object_or_404(CachedGeneratedFile, pk=cgf_id)
    cgf.delete()
    return HttpResponse()


@require_POST
def messages_bulk_delete(request):
    messages_str = request.POST['message_ids']
    message_ids = json.loads(messages_str)
    user_messages_qs = Message.objects.filter(recipient=request.user)
    user_messages_qs.filter(pk__in=message_ids).delete()
    return HttpResponse()


def help_static_page(request, page_name):
    """ This embeds static pages in a help template """
    context = {"page_name": page_name}
    return render(request, 'snpdb/help/help_static_page.html', context)


def staff_only(request):
    return render(request, 'snpdb/staff_only.html')


def wait_for_task(request, celery_task, sleep_ms, redirect_url):
    async_result = AsyncResult(celery_task)
    if async_result.successful():
        return HttpResponseRedirect(redirect_url)

    kwargs = {"celery_task": celery_task, "sleep_ms": sleep_ms, "redirect_url": redirect_url}
    url = reverse("wait_for_task", kwargs=kwargs)

    context = {"url": url,
               "sleep_ms": sleep_ms,
               "async_result": async_result}
    return render(request, 'snpdb/wait_for_task.html', context)


@require_POST
def wiki_save(request, class_name, unique_keyword, unique_value):
    # get_or_create checks write permission (throws 403) before creating the row
    wiki = Wiki.get_or_create(class_name, unique_keyword, unique_value, user=request.user)
    markdown = request.POST["markdown"]

    wiki.markdown = markdown
    wiki.last_edited_by = request.user
    wiki.save()
    return JsonResponse({})


def view_genome_build(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    context = {
        "genome_build": genome_build,
        "standard_contigs": genome_build.standard_contigs,
    }
    return render(request, "snpdb/genomics/view_genome_build.html", context)


def view_contig(request, contig_accession):
    q = Contig.get_q(contig_accession)
    contig = get_object_or_404(Contig, q)
    # Prefer builds with annotation
    builds_with_annotation = list(contig.get_genome_builds(require_annotation=True))
    if builds_with_annotation:
        has_annotation = True
        builds = builds_with_annotation
    else:
        has_annotation = False
        builds = list(contig.get_genome_builds(require_annotation=False))

    if builds:
        # If multiple builds (eg MT), we'll use user default if they have it, doesn't really matter
        genome_build = builds[0]
        if len(builds) > 1:
            try:
                user_settings = UserSettings.get_for_user(request.user)
                if user_settings.default_genome_build in builds:
                    genome_build = user_settings.default_genome_build
            except:
                pass
    else:
        genome_build = None

    context = {
        "contig": contig,
        "genome_build": genome_build,
        "has_annotation": has_annotation,
    }
    return render(request, "snpdb/genomics/view_contig.html", context)


@login_not_required
def view_uptime(request):
    uptime_response = uptime_check.retrieve_uptime_response()
    return render(request, "snpdb/uptime_check.html", {"uptime_response": uptime_response})
