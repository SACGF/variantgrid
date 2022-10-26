
from random import randint

import celery
from django.conf import settings
from django.urls import reverse

from library.django_utils import get_url_from_view_path
from library.health_check import populate_health_check
from library.log_utils import NotificationBuilder


# all the dependencies required by injecting everything into this, no longer needed using signals
# from django.utils.timezone import now
# from typing import Optional
# from annotation.models import ClinVarVersion
# from classification.models import Classification
# from flags.models.flag_health_check import flag_chanced_since
# from collections import defaultdict
# from datetime import timedelta
# from library.guardian_utils import admin_bot
# from library.utils import pretty_label, count
# from ontology.models import OntologyImport
# from snpdb.models import GenomeBuild
# from variantgrid.perm_path import get_visible_url_names
# from variantgrid.tasks.server_monitoring_tasks import get_disk_messages
# from variantopedia.server_status import get_dashboard_notices


@celery.shared_task
def notify_server_status():
    if not settings.HEALTH_CHECK_ENABLED:
        return
    notify_server_status_now()


def notify_server_status_now():
    url = get_url_from_view_path(reverse('server_status')) + '?activeTab=server_status_activity_detail_1'
    nb = NotificationBuilder(message="Health Check")
    heading_emoji = ":male-doctor:" if randint(0, 1) else ":female-doctor:"
    nb.add_header(f"{heading_emoji} Health Check")
    nb.add_markdown(f"*In the <{url}|last 24 hours>*")
    populate_health_check(nb)
    nb.send()


# def notify_server_status_now(detailed: bool = True):
#     dashboard_notices = get_dashboard_notices(admin_bot(), days_ago=1)
#     url = get_url_from_view_path(reverse('server_status')) + '?activeTab=server_status_activity_detail_1'
#     right_now = now()
#     day_ago = right_now - timedelta(days=1)
#
#     heading_emoji = ":male-doctor:" if randint(0, 1) else ":female-doctor:"
#     nb = NotificationBuilder(message="Health Check")
#
#     keys = set(dashboard_notices.keys())
#     keys.discard('events')
#     keys.discard('notice_header')
#     visible_urls = get_visible_url_names()
#     if not visible_urls.get('analyses'):
#         for exclude_key in ['vcfs', 'analyses_created', 'analyses_modified']:
#             keys.discard(exclude_key)
#     sorted_keys = sorted(list(keys))
#
#     lines = []]
#     zeros = []
#
#     for key in sorted_keys:
#         values = dashboard_notices.get(key)
#         count_display = count(values)
#
#         display_individuals = False
#         emoji = ":blue_book:"
#         if 'analyses' in key:
#             emoji = ":orange_book:"
#         elif 'vcf' in key:
#             emoji = ":green_book:"
#         elif 'active_users' in key:
#             if count_display == 0:
#                 emoji = ":ghost:"
#             else:
#                 emojis = [":nerd_face:", ":thinking_face:", ":face_with_monocle:", ":face_with_cowboy_hat:"]
#                 emoji = emojis[randint(0, len(emojis) - 1)]
#                 display_individuals = True
#
#         if count_display:
#             count_display = f"*{count_display}*"
#         elif 'active_users' not in key:
#             zeros.append(pretty_label(key))
#             continue
#
#         line = f"{emoji} {count_display} : {pretty_label(key)}"
#         if display_individuals:
#             line = f"{line} : {', '.join(str(value) for value in values)}"
#
#         lines.append(line)
#     if zeros:
#         lines.append(f":open_file_folder: 0 : {', '.join(zeros)}")
#
#     # only put this in health check and not regular
#     # server status since
#     flag_deltas = flag_chanced_since(day_ago)
#     for flag_delta in flag_deltas:
#         parts = []
#         if flag_delta.added:
#             parts.append(f"*+{flag_delta.added}*")
#         if flag_delta.resolved:
#             parts.append(f"*-{flag_delta.resolved}*")
#         joined_parts = "/".join(parts)
#         lines.append(f":flags: {joined_parts} : Flag {flag_delta.flag_type}")
#
#     nb.add_header(f"{heading_emoji} Health Check")
#     nb.add_markdown(f"*In the <{url}|last 24 hours>*")
#
#     # just need to sort out the order and grouping of health check
#     # and it can replace all of these hardcoded references
#     populate_health_check(nb)
#
#     nb.add_markdown("\n".join(lines), indented=True)
#
#     if detailed:
#         overall_lines = []
#         for _, message in get_disk_messages(info_messages=True):
#             overall_lines.append(f":floppy_disk: {message}")
#         if not overall_lines:
#             overall_lines.append(":floppy_disk: _Disk Usage Unknown_")
#
#         def emoji_for_age(days: int) -> str:
#             if days <= 60:
#                 return ":smile:"
#             if days <= 120:
#                 return ":neutral_face:"
#             if days <= 180:
#                 return ":cry:"
#             if days <= 365:
#                 return ":rage:"
#             return ":exploding_head:"
#
#         def message_for_annotation(annotation_name: str, days: int) -> Optional[str]:
#             # this is the window where we've seen that the annotations have been updated
#             # and not old enough for us to worry about again yet
#             if 2 <= days <= 59:
#                 return None
#             return f"{emoji_for_age(days)} {days} day{'s' if days != 1 else ''} old : {annotation_name}"
#
#         annotation_ages = []
#         # others we might want to check date of
#
#         # when we have an array, only report on the first instance that's been imported
#         # e.g. if we import OMIM directly from OMIM, we don't care how old the biomart omim file is
#
#         age_days_to_annotations = defaultdict(list)
#         for contexts in [("mondo_file", "MONDO"), [("omim_file", "OMIM"), ("biomart_omim_aliases", "OMIM - via biomart")], ("hpo_file", "HPO"), ("gencc_file", "GenCC")]:
#             if not isinstance(contexts, list):
#                 contexts = [contexts]
#             for context, label in contexts:
#                 if last_import := OntologyImport.objects.filter(context=context).order_by('-processed_date').first():
#                     time_delta = right_now - last_import.processed_date
#                     age_days_to_annotations[time_delta.days].append(label)
#                     break
#
#         for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
#             if latest_clinvar := ClinVarVersion.objects.filter(genome_build=genome_build).order_by('-annotation_date').first():
#                 time_delta = right_now - latest_clinvar.annotation_date
#                 age_days_to_annotations[time_delta.days].append(f"ClinVar {genome_build}")
#
#         for days in sorted(age_days_to_annotations.keys()):
#             annotations = age_days_to_annotations[days]
#             if message := message_for_annotation(", ".join(sorted(annotations)), days):
#                 annotation_ages.append(message)
#
#         overall_lines.extend(annotation_ages)
#
#         nb.add_markdown("*Overall*")
#         nb.add_markdown("\n".join(overall_lines), indented=True)
#     nb.send()
