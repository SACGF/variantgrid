"""
Miscellaneous site pages that belong to no other snpdb views module: the index, help/static pages, `wait_for_task`
(polls a celery task, then redirects to a same-site URL), wiki saves, genome build and contig pages, deleting a
CachedGeneratedFile, bulk message delete and the public `view_uptime` check.
"""
import json

from celery.result import AsyncResult
from django.contrib.auth.decorators import login_not_required
from django.core.exceptions import BadRequest
from django.http.response import (
    HttpResponse,
    HttpResponseRedirect,
    JsonResponse,
)
from django.shortcuts import get_object_or_404, render
from django.urls.base import reverse
from django.utils.http import url_has_allowed_host_and_scheme
from django.views.decorators.http import require_POST
from termsandconditions.decorators import terms_required

from library import uptime_check
from snpdb.cytobands import parse_cytoband
from snpdb.models import (
    CachedGeneratedFile,
    Contig,
    GenomeBuild,
    UserSettings,
    Wiki,
)
from snpdb.models.models_genome import Cytoband
from snpdb.models.models_user_settings import get_igv_data
from snpdb.serializers import ContigSerializer, CytobandSerializer
from snpdb.variant_filters import get_default_all_variants_filters, get_variant_type_choices
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
    if not url_has_allowed_host_and_scheme(redirect_url, allowed_hosts={request.get_host()},
                                           require_https=request.is_secure()):
        raise BadRequest("redirect_url must be on this site")
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
    data = {}
    if (citation_ids := wiki.citation_ids()) is not None:
        data["citation_ids"] = citation_ids
    return JsonResponse(data)


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
            user_settings = UserSettings.get_for_user(request.user)
            if user_settings.default_genome_build in builds:
                genome_build = user_settings.default_genome_build
    else:
        genome_build = None

    # ?band=p13 (from cytoband search) or a range, ?band=p13-p12
    band = request.GET.get("band")
    cytoband_interval = None
    if band and genome_build:
        if cytoband_range := parse_cytoband(f"{contig.name}{band}"):
            cytoband_interval = cytoband_range.get_interval(genome_build)

    # The file's bands - the levels above them are generated, with no stain
    cytobands = list(Cytoband.objects.filter(contig=contig).exclude(stain="").select_related("contig").order_by("start"))
    if cytoband_interval:
        for cytoband in cytobands:
            cytoband.in_selected_band = cytoband_interval.start <= cytoband.start and cytoband.end <= cytoband_interval.end

    selected_interval = None
    if cytoband_interval:
        selected_interval = {"start": cytoband_interval.start, "end": cytoband_interval.end}

    context = {
        "contig": contig,
        "genome_build": genome_build,
        "has_annotation": has_annotation,
        "band": band,
        "cytoband_interval": cytoband_interval,
        "cytobands": cytobands,
        "contig_page_data": {
            "contig": ContigSerializer(contig).data,
            "cytobands": CytobandSerializer(cytobands, many=True).data,
            "band": band if selected_interval else None,
            "selectedInterval": selected_interval,
        },
    }
    if has_annotation:
        # The All Variants grid, locked to this contig - not saved as the user's All Variants page filters
        position_range = [selected_interval["start"], selected_interval["end"]] if selected_interval else None
        locked_filters = {"contig_ids": [contig.pk], "position_range": position_range}
        variant_types = get_variant_type_choices()
        context["contig_page_data"]["numVariantTypes"] = len(variant_types)
        context.update({
            "variant_types": variant_types,
            "all_variants_filters": get_default_all_variants_filters(genome_build) | locked_filters,
            "all_variants_locked_filters": locked_filters,
            "igv_data": get_igv_data(request.user, genome_build=genome_build),
        })
    return render(request, "snpdb/genomics/view_contig.html", context)


@login_not_required
def view_uptime(request):
    uptime_response = uptime_check.retrieve_uptime_response()
    return render(request, "snpdb/uptime_check.html", {"uptime_response": uptime_response})
