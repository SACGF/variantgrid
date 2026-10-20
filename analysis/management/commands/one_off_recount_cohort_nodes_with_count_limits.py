"""
Reloads every CohortNode whose zygosity count limits filter something. Their counts came from the cohort's
unfiltered precomputed stats rather than their own query (#2370), and a node keeps its counts until its
version changes - saving it queryset-dirty bumps it and its descendants, which reload next time they're viewed.
"""
from django.core.management import BaseCommand
from django.db.models import Q

from analysis.models import CohortNode


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        any_limit_q = Q()
        for _, _, min_field, max_field in CohortNode.COUNT_COLUMNS:
            any_limit_q |= Q(**{f"{min_field}__gt": 0}) | Q(**{f"{max_field}__isnull": False})

        qs = CohortNode.objects.filter(any_limit_q, accordion_panel=CohortNode.COUNT, cohort__isnull=False)
        num_reloaded = 0
        for node in qs.select_related("cohort").order_by("pk"):
            if not node.has_effective_count_limits():
                continue
            self.stdout.write(f"{'Would reload' if dry_run else 'Reloading'} {node.analysis_id}/{node.pk} "
                              f"({node.count} counted): {node._get_zygosity_count_description()}")
            if not dry_run:
                node.queryset_dirty = True
                node.save()
            num_reloaded += 1
        self.stdout.write(f"{'Would reload' if dry_run else 'Reloaded'} {num_reloaded} cohort nodes")
