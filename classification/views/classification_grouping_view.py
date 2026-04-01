from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Iterable, Optional, Self
from xml.sax.handler import property_dom_node

from django.http import HttpRequest
from django.shortcuts import render
from more_itertools.more import first
from requests import Response

from classification.enums import OverlapStatus, TestingContextBucket
from classification.models import ClassificationGrouping, Overlap, ClassificationResultValue, OverlapContributionStatus, \
    OverlapType, OverlapContribution, OverlapContributionSkew
from classification.services.overlaps_services import OverlapEntryCompare, OverlapGrouping


@dataclass(frozen=True)
class OverlapEntry:
    overlap_contribution: OverlapContribution
    context: OverlapType
    is_primary: bool = False

    @property
    def is_primary_context(self) -> bool:
        return self.context == OverlapType.SINGLE_CONTEXT

    @property
    def sort_index(self):
        return (self.context.priority_order,
                self.overlap_contribution.value_sort_index,
                self.overlap_contribution.lab.name if self.overlap_contribution.lab else "")

    def __lt__(self, other):
        return self.sort_index < other.sort_index


@dataclass(frozen=True)
class OverlapSummary:
    skews: list[OverlapContributionSkew]
    value_type: ClassificationResultValue
    overlaps: list[Overlap]

    @property
    def value_type_label(self):
        return first(self.overlaps).value_type_label

    @cached_property
    def overlap_entries(self) -> list[OverlapEntry]:
        used_contributions = set([first(self.skews).contribution.pk])
        overlap_entries = [OverlapEntry(overlap_contribution=first(self.skews).contribution, context=OverlapType.SINGLE_CONTEXT, is_primary=True)]

        overlaps = sorted(self.overlaps, key=lambda x: OverlapType(x.overlap_type))
        for overlap in overlaps:
            for contribution in overlap.contributions:
                if contribution.pk not in used_contributions:
                    if contribution.value_type == self.value_type:
                        used_contributions.add(contribution.pk)
                        overlap_entries.append(
                            OverlapEntry(overlap_contribution=contribution, context=OverlapType(overlap.overlap_type))
                        )
        return sorted(overlap_entries)

    @staticmethod
    def overlap_summaries_for(cg: ClassificationGrouping) -> list[Self]:
        overlap_contributions = OverlapContribution.objects.filter(classification_grouping=cg,
                                                                   contribution_status=OverlapContributionStatus.CONTRIBUTING)
        skews = OverlapContributionSkew.objects.filter(contribution__in=overlap_contributions)
        skew_by_value_type = defaultdict(list)
        overlap_summary_builders: dict[ClassificationResultValue, list] = defaultdict(list)
        summaries = []
        for skew in skews:
            skew_by_value_type[skew.contribution.value_type].append(skew)

            overlap = skew.overlap
            if overlap.valid:
                overlap_summary_builders[overlap.value_type].append(overlap)

        for value_type, overlaps in overlap_summary_builders.items():
            if skews := skew_by_value_type.get(value_type):
                summaries.append(
                    OverlapSummary(skews=skews, value_type=value_type, overlaps=overlaps)
                )
            else:
                print(f"No skews for {value_type}")
        return summaries

def view_overlaps_for_classification_grouping(request: HttpRequest, classification_grouping_id: int) -> Response:
    cg = ClassificationGrouping.objects.filter(pk=classification_grouping_id).get()
    # overlap_contributions = OverlapContribution.objects.filter(classification_grouping=cg, contribution_status=OverlapContributionStatus.CONTRIBUTING)
    # skews = OverlapContributionSkew.objects.filter(contribution__in=overlap_contributions)
    summaries = OverlapSummary.overlap_summaries_for(cg)

    return render(request, 'classification/snippets/overlaps_for_classification_grouping_detail.html', {
        "allele": cg.allele,
        "summaries": summaries
    })
