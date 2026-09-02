import logging

from django.contrib import messages
from django.http.response import (
    HttpResponse,
)
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST

from library.constants import WEEK_SECS
from library.django_utils import (
    add_save_message,
    require_superuser,
)
from snpdb import forms
from snpdb.forms import (
    TagForm,
)
from snpdb.models import (
    TAG_ALLELE_ORIGIN_CHOICES,
    Tag,
    TagColorsCollection,
    UserSettings,
)
from snpdb.tag_operations import (
    TagOperation,
    describe_tag_operation,
    get_merge_suggestions,
    get_tag_operations,
    get_tag_usage,
    get_tag_usage_by_tag_id,
    merge_tag,
    reinstate_tag,
    retire_tag,
    set_tag_allele_origin,
)
from snpdb.utils import get_tag_styles_and_colors


TAG_OPERATION_HISTORY_LIMIT = 50


def _get_tag_operation_history() -> list[dict]:
    """ Recent operations on the tag vocabulary, for the tag settings page """
    return [{
        "tag_id": log_entry.object_pk,
        "operation": TagOperation(log_entry.additional_data["operation"]).label,
        "detail": describe_tag_operation(log_entry),
        "user": log_entry.actor,
        "timestamp": log_entry.timestamp,
    } for log_entry in get_tag_operations()[:TAG_OPERATION_HISTORY_LIMIT]]


def _get_tag_summaries() -> list[dict]:
    """ Per-tag usage counts + possible typo twins, for the tag settings table """
    usage_by_tag_id = get_tag_usage_by_tag_id()
    merge_suggestions = get_merge_suggestions()
    tags_by_id = Tag.live_qs().in_bulk()
    return [{
        "tag": tags_by_id[tag_id],
        "tag_id": tag_id,
        "usage": usage,
        "suggestions": merge_suggestions.get(tag_id, []),
    } for tag_id, usage in sorted(usage_by_tag_id.items())]


def tag_settings(request):
    form = forms.CreateTagForm(request.POST or None)
    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            tag = form.save()
            name = f"Tag {tag}"
        else:
            name = "Tag"
        add_save_message(request, valid, name, created=True)

    user_settings = UserSettings.get_for_user(request.user)
    user_settings_tag_colors = user_settings.tag_colors
    user_tag_styles, user_tag_colors = get_tag_styles_and_colors(request.user, user_settings_tag_colors)
    context_dict = {
        'form': form,
        "user_settings_tag_colors": user_settings_tag_colors,
        'user_tag_styles': user_tag_styles,
        'user_tag_colors': user_tag_colors,
        'tag_summaries': _get_tag_summaries(),
        'allele_origin_choices': TAG_ALLELE_ORIGIN_CHOICES,
        'retired_tags': Tag.objects.filter(retired__isnull=False).order_by("-retired"),
        'tag_operations': _get_tag_operation_history(),
    }
    return render(request, 'snpdb/settings/tag_settings.html', context_dict)


@require_superuser
def tag_merge(request, tag_id):
    """ Merges the tag in the URL into another - always stated from the dying tag's side """
    dying_tag = get_object_or_404(Tag, pk=tag_id, retired__isnull=True)
    surviving_tag = None
    if surviving_tag_id := request.GET.get("into") or request.POST.get("surviving_tag"):
        surviving_tag = Tag.live_qs().filter(pk=surviving_tag_id).exclude(pk=dying_tag.pk).first()

    if request.method == "POST":
        confirm = request.POST.get("confirm_tag_name", "")
        if surviving_tag is None:
            messages.add_message(request, messages.ERROR, "Choose an existing tag to merge into")
        elif confirm != dying_tag.pk:
            messages.add_message(request, messages.ERROR,
                                 f"Type '{dying_tag}' exactly to confirm the merge")
        else:
            result = merge_tag(dying_tag, surviving_tag, request.user)
            messages.add_message(request, messages.INFO, result.description())
            return redirect('tag_settings')

    context = {
        "dying_tag": dying_tag,
        "surviving_tag": surviving_tag,
        "usage": get_tag_usage(dying_tag),
        "other_tags": Tag.live_qs().exclude(pk=dying_tag.pk).order_by("pk"),
        "suggestions": get_merge_suggestions().get(dying_tag.pk, []),
    }
    return render(request, 'snpdb/settings/tag_merge.html', context)


@require_superuser
@require_POST
def tag_retire(request, tag_id):
    """ Stop offering a tag - what already uses it is left alone """
    tag = get_object_or_404(Tag, pk=tag_id)
    try:
        retire_tag(tag, request.user, reason=request.POST.get("reason") or None)
        messages.add_message(request, messages.INFO,
                             f"Retired '{tag}' - it's no longer offered, existing tagging is unchanged")
    except ValueError as ve:
        messages.add_message(request, messages.ERROR, str(ve))
    return redirect('tag_settings')


@require_superuser
@require_POST
def tag_set_allele_origin(request, tag_id):
    """ Which side of the house a tag belongs to - the tag stats page filters on it """
    tag = get_object_or_404(Tag, pk=tag_id)
    try:
        set_tag_allele_origin(tag, request.POST.get("allele_origin_bucket"), request.user)
        messages.add_message(request, messages.INFO,
                             f"'{tag}' allele origin is now {tag.get_allele_origin_bucket_display()}")
    except ValueError as ve:
        messages.add_message(request, messages.ERROR, str(ve))
    return redirect('tag_settings')


@require_superuser
@require_POST
def tag_reinstate(request, tag_id):
    """ Put a retired tag back into circulation """
    tag = get_object_or_404(Tag, pk=tag_id)
    try:
        reinstate_tag(tag, request.user)
        messages.add_message(request, messages.INFO, f"Reinstated '{tag}'")
    except ValueError as ve:
        messages.add_message(request, messages.ERROR, str(ve))
    return redirect('tag_settings')


def view_tag_colors_collection(request, tag_colors_collection_id):
    tag_colors_collection = TagColorsCollection.get_for_user(request.user, pk=tag_colors_collection_id)
    has_write_permission = tag_colors_collection.can_write(request.user)
    if not has_write_permission:
        msg = "You do not have permission to edit these columns. " \
              "If you wish to customise them, click 'clone' and modify the copy"
        messages.add_message(request, messages.WARNING, msg)

    if request.method == "POST":
        tag_colors_collection.check_can_write(request.user)
        if name := request.POST.get("name"):
            tag_colors_collection.name = name
            tag_colors_collection.save()
        elif tag_order_str := request.POST.get("tag_order"):
            valid_tag_ids = set(Tag.objects.values_list("pk", flat=True))
            for i, tag_id in enumerate(tag_order_str.split(",")):
                if tag_id in valid_tag_ids:
                    tag_colors_collection.tagcolor_set.update_or_create(tag_id=tag_id,
                                                                        defaults={"sort_order": i})
        return HttpResponse()  # Nobody ever looks at this

    user_tag_styles, user_tag_colors = get_tag_styles_and_colors(request.user, tag_colors_collection)
    context = {
        "tag_colors_collection": tag_colors_collection,
        "has_write_permission": has_write_permission,
        'user_tag_styles': user_tag_styles,
        'user_tag_colors': user_tag_colors,
    }
    return render(request, 'snpdb/settings/view_tag_colors_collection.html', context)


@require_POST
def set_tag_color(request, tag_colors_collection_id):
    tcc = TagColorsCollection.get_for_user(request.user, pk=tag_colors_collection_id)
    tag = request.POST['tag']
    rgb = request.POST['rgb']
    tc = tcc.tagcolor_set.get_or_create(tag_id=tag)[0]
    tc.rgb = rgb
    tc.save()
    logging.info("set_tag_color: saved %s", tc)
    return HttpResponse()


@cache_page(WEEK_SECS)
def tag_autocomplete_form(request):
    """ This is an absolutely minimal HTML to create a Tag autocomplete form (used for load()) """

    context = {"tag_form": TagForm()}
    return render(request, 'snpdb/tag_autocomplete_form.html', context)
