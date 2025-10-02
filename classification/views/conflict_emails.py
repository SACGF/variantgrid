### EMAILS
from collections import Counter
from typing import Optional, Iterable

from more_itertools.more import last
from pysam.libcvcf import defaultdict
from rest_framework.reverse import reverse

from classification.enums import AlleleOriginBucket
from classification.models import ConflictNotificationStatus, ConflictNotificationRun, ConflictNotification, \
    ConflictHistory
from classification.services.conflict_services import ConflictCompareType, ConflictCompare, ConflictDataRow
from classification.views.conflict_view import ConflictFeed
from library.django_utils import get_url_from_view_path
from snpdb.models import Lab
from snpdb.utils import LabNotificationBuilder


class ConflictCompareFeed(ConflictFeed):

    def __init__(self, conflict_compare: ConflictCompare):
        self.conflict_compare = conflict_compare
        super().__init__(conflict=conflict_compare.conflict, user=None)

    def history_generator(self) -> Iterable[ConflictHistory]:
        if previous := self.conflict_compare.previous_state:
            return [previous, self.conflict_compare.current_state]
        else:
            return [self.conflict_compare.current_state]


def send_emails_for_conflict_notification_run(notification_run: ConflictNotificationRun):
    conflict_notifications = sorted(ConflictNotification.objects.filter(notification_run=notification_run))

    lab_to_conflict_compares: dict[int, list[ConflictCompare]] = defaultdict(list)
    for cn in conflict_notifications:
        cc = ConflictCompare.from_conflict_notification(cn)
        if cc.is_notifiable_difference:
            for lab_id in cc.involved_lab_ids:
                lab_to_conflict_compares[lab_id].append(cc)

    # now we have the notifications to send...
    for lab_id, compares in lab_to_conflict_compares.items():
        change_types: dict[ConflictCompareType, int] = Counter()
        gene_symbols: set[str] = set()
        allele_origins: set[AlleleOriginBucket] = set()
        for conflict_compare in compares:
            notification_change_type = conflict_compare.notification_change_type
            if notification_change_type != ConflictCompareType.NO_CHANGE:
                allele_origins.add(AlleleOriginBucket(conflict_compare.current_state.conflict.allele_origin_bucket))
                change_types[conflict_compare.notification_change_type] += 1
                found_gene_symbol: Optional[str] = None
                if c_hgvses := conflict_compare.conflict.c_hgvses():
                    for c_hgvs in c_hgvses:
                        if gene_symbol := c_hgvs.gene_symbol:
                            found_gene_symbol = gene_symbol
                            break
                if found_gene_symbol:
                    gene_symbols.add(found_gene_symbol)
                else:
                    gene_symbols.add("Unknown Gene Symbol")

        lab = Lab.objects.get(id=lab_id)
        if len(compares) > 0:

            gene_symbol_str = ", ".join(sorted(gene_symbols))
            allele_origin_str = " and ".join(ao.label for ao in sorted(allele_origins))
            title: str
            if len(compares) == 1:
                title = f"Discordance update for {allele_origin_str}. Gene symbol ({gene_symbol_str}) - 1 : {compares[0].notification_change_type.label}"
            else:
                update_parts = []
                for change in ConflictCompareType:
                    if change_count := change_types.get(change):
                        update_parts.append(f"{change_count} : {change.label}")
                update_parts_str = ", ".join(update_parts)
                title = f"Discordance updates for {allele_origin_str}. Gene symbols ({gene_symbol_str}) - {update_parts_str}"

            lab_notification = LabNotificationBuilder(lab, title)

            for compare in compares:
                lab_notification.add_divider()
                # TODO move this into code that can be re-used for telling admins
                lab_notification.add_markdown(f"Discordance Report (CR_{compare.conflict.pk}) Update")

                conflict_link = get_url_from_view_path(
                    reverse('conflict', kwargs={'conflict_id': compare.conflict.pk})
                )
                parts = []
                parts.append(f"<{conflict_link}|CR_{compare.conflict.pk}>")
                parts.append(compare.conflict.context_summary)
                parts.append(compare.notification_change_type.label)

                lab_notification.add_field("Discordance Detail", " - ".join(parts))

                feed_item = last(ConflictCompareFeed(compare).history_iterator())
                if explanations := feed_item.extra_comments:
                    for explanation in explanations:
                        lab_notification.add_markdown(explanation)

                row: ConflictDataRow
                for row in feed_item.data_rows:
                    if row.exclude:
                        continue
                    lab = row.lab
                    # TODO add lab date here
                    value = row.value_label
                    lab_notification.add_markdown(f"**{lab}**: {value}")
                    if extra_comment := row.message:
                        lab_notification.add_markdown(extra_comment)
                    lab_notification.add_markdown(f"\n")

            lab_notification.send()

    notification_run.status = ConflictNotificationStatus.COMPLETE
    notification_run.save()