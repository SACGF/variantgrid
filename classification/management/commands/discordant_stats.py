from collections import defaultdict
from datetime import date
from enum import StrEnum

from django.core.management.base import BaseCommand

from classification.enums import OverlapType, OverlapStatus, OverlapOverrideStatus, ClassificationResultValue
from classification.models import Overlap, OverlapContribution
from classification.services.overlap_calculator import OverlapCalculatorOncPath


class OverlapCategories(StrEnum):
    COMPLEX = "COMPLEX"
    CONTINUED = "CONTINUED"
    CLINVAR_INCLUDED = "CLINVAR_INCLUDED"
    CLINVAR_NEWER_DISCORDANCE = "CLINVAR_NEWER_DISCORDANCE"
    CLINVAR_OLDER_DISCORDANCE = "CLINVAR_OLDER_DISCORDANCE"
    REGULAR = "REGULAR"


class Command(BaseCommand):

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        overlap_counts = defaultdict(set)
        clinvar_date_diffs = []
        for overlap in Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                valid=True,
                overlap_status__gte=OverlapStatus.MAJOR_DIFFERENCES,
                value_type=ClassificationResultValue.ONC_PATH
        ).iterator():
            if override_status := overlap.overlap_override_status:
                if override_status == OverlapOverrideStatus.COMPLEX:
                    overlap_counts[OverlapCategories.COMPLEX].add(overlap.pk)
                elif override_status == OverlapOverrideStatus.CONTINUED_DISCORDANCE:
                    overlap_counts[OverlapCategories.CONTINUED].add(overlap.pk)
            else:
                clinvar_contributions: list[OverlapContribution] = []
                non_clinvar_contributions: list[OverlapContribution] = []
                for contribution in overlap.contributions_list:
                    if contribution.scv:
                        clinvar_contributions.append(contribution)
                    else:
                        non_clinvar_contributions.append(contribution)
                if len(clinvar_contributions) == 0:
                    overlap_counts[OverlapCategories.REGULAR].add(overlap.pk)
                else:
                    # clinvar_values = {cc.value for cc in clinvar_contributions}
                    non_clinvar_values = {cc.value for cc in non_clinvar_contributions}

                    #combined_status = OverlapCalculatorOncPath.calculate_status_for_multiple_entries(clinvar_values.union(non_clinvar_values))
                    sans_clinvar_status = OverlapStatus.SINGLE_SUBMITTER
                    if len(non_clinvar_contributions) > 1:
                        if len(non_clinvar_values) == 1:
                            sans_clinvar_status = OverlapStatus.EXACT_AGREEMENT
                        else:
                            sans_clinvar_status = OverlapCalculatorOncPath.calculate_status_for_multiple_entries(non_clinvar_values)

                    if sans_clinvar_status.is_discordant:
                        # the overlap has ClinVar Expert Panel in it, but would still be discordant without it
                        overlap_counts[OverlapCategories.CLINVAR_INCLUDED].add(overlap.pk)
                    else:
                        latest_lab_curated = max(cont.effective_date_obj for cont in non_clinvar_contributions)
                        latest_clinvar_curated = max(cont.effective_date_obj for cont in clinvar_contributions)

                        date_diff = date.fromisoformat(latest_lab_curated.date) - date.fromisoformat(latest_clinvar_curated.date)
                        clinvar_date_diffs.append(date_diff)

                        if latest_clinvar_curated < latest_lab_curated:
                            overlap_counts[OverlapCategories.CLINVAR_OLDER_DISCORDANCE].add(overlap.pk)
                        else:
                            overlap_counts[OverlapCategories.CLINVAR_NEWER_DISCORDANCE].add(overlap.pk)

        for dd in sorted(clinvar_date_diffs):
            print(dd)

        for category, overlaps in overlap_counts.items():
            print(f"{category} : {len(overlaps)} - e.g. {list(overlaps)[0:5]}")
