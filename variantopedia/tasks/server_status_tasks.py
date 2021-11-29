from random import randint

import celery
from django.conf import settings
from django.urls import reverse
from django.utils.timezone import now

from classification.models import Classification
from library.django_utils import get_url_from_view_path
from library.guardian_utils import admin_bot
from library.log_utils import NotificationBuilder
from library.utils import pretty_label, count
from ontology.models import OntologyImport
from variantgrid.perm_path import get_visible_url_names
from variantgrid.tasks.server_monitoring_tasks import get_disk_messages
from variantopedia.server_status import get_dashboard_notices


@celery.shared_task
def notify_server_status():
    if not settings.HEALTH_CHECK_ENABLED:
        return
    notify_server_status_now()


def notify_server_status_now(detailed: bool = True):
    dashboard_notices = get_dashboard_notices(admin_bot(), days_ago=1)
    url = get_url_from_view_path(reverse('server_status')) + '?activeTab=server_status_activity_detail_1'

    emoji = ":male-doctor:" if randint(0, 1) else ":female-doctor:"
    nb = NotificationBuilder(message="Health Check", emoji=emoji)

    keys = set(dashboard_notices.keys())
    keys.discard('events')
    keys.discard('notice_header')
    visible_urls = get_visible_url_names()
    if not visible_urls.get('analyses'):
        for exclude_key in ['vcfs', 'analyses_created', 'analyses_modified']:
            keys.discard(exclude_key)
    sorted_keys = sorted(list(keys))

    lines = list()
    zeros = list()

    for key in sorted_keys:
        values = dashboard_notices.get(key)
        count_display = count(values)

        display_individuals = False
        emoji = ":blue_book:"
        if 'analyses' in key:
            emoji = ":orange_book:"
        elif 'vcf' in key:
            emoji = ":green_book:"
        elif 'active_users' in key:
            if count_display == 0:
                emoji = ":ghost:"
            else:
                emojis = [":nerd_face:", ":thinking_face:", ":face_with_monocle:", ":face_with_cowboy_hat:"]
                emoji = emojis[randint(0, len(emojis) - 1)]
                display_individuals = True

        if count_display:
            count_display = f"*{count_display}*"
        elif not 'active_users' in key:
            zeros.append(pretty_label(key))
            continue

        line = f"{emoji} {count_display} : {pretty_label(key)}"
        if display_individuals:
            line = f"{line} : {', '.join(str(value) for value in values)}"

        lines.append(line)
    if zeros:
        lines.append(f":open_file_folder: 0 : {', '.join(zeros)}")

    nb.add_header("Health Check")
    nb.add_markdown(f"*In the <{url}|last 24 hours>*")
    nb.add_markdown("\n".join(lines), indented=True)

    if detailed:
        overall_lines = list()
        for _, message in get_disk_messages(info_messages=True):
            overall_lines.append(f":floppy_disk: {message}")
        if not overall_lines:
            overall_lines.append(f":floppy_disk: _Disk Usage Unknown_")

        total_shared = Classification.dashboard_total_shared_classifications()
        total_unshared = Classification.dashboard_total_unshared_classifications()
        total = total_unshared + total_shared
        if total:
            percent_shared = 100.0 * float(total_shared) / float(total)
            overall_lines.append(
                f":blue_book: {total:,} : Classifications - {int(percent_shared)}% shared"
            )

        def emoji_for_age(days: int):
            if days <= 60:
                return ":smile:"
            if days <= 120:
                return ":neutral_face:"
            if days <= 180:
                return ":cry:"
            return ":rage:"

        annotation_ages = list()
        # others we might want to check date of
        right_now = now()
        for context, label in [("mondo_file", "MONDO"), ("omim_file", "OMIM")]:
            if last_import := OntologyImport.objects.filter(context=context).order_by('-created').first():
                time_delta = right_now - last_import.created
                annotation_ages.append(f"{emoji_for_age(time_delta.days)} {time_delta.days} days old : {label}")

        overall_lines.extend(annotation_ages)

    # TODO add other important annotation ages like ClinVar

        nb.add_markdown("*Overall*")
        nb.add_markdown("\n".join(overall_lines), indented=True)
    nb.send()
