from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from typing import Optional

from pysam.libcvcf import defaultdict

from classification.enums import AlleleOriginBucket, OverlapStatus
from classification.models import ClassificationGrouping, AlleleOriginGrouping, Overlap, OverlapType, \
    ClassificationResultValue, ClassificationSummaryCacheDict, ClassificationGroupingOverlapContribution, \
    OverlapContributionStatus, EvidenceKeyMap
from classification.services.overlap_calculator_clinsig import OverlapCalculatorClinSig
from classification.services.overlap_calculator_oncpath import OverlapCalculatorOncPath


class OverlapServices:

    """
    overlap_type = models.TextField(choices=OverlapType.choices)
    value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
    overlap_status = models.TextField(max_length=1, choices=OverlapStatus.choices, default=OverlapStatus.NO_CONTRIBUTIONS.value)

    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
    testing_context = models.TextField(max_length=1, choices=TestingContextBucket.choices, null=True, blank=True)
    condition = models.TextField(null=True, blank=True)  # condition isn't always relevant
    lab = models.ForeignKey(Lab, on_delete=models.CASCADE, null=True, blank=True)  # only use for
    """

    def calculate_and_apply_overlaps_for_ao(self, allele_origin_grouping: AlleleOriginGrouping):
        # get the default context for the classification overlaps
        context_overlaps = []
        onc_path_context_overlap, _ = Overlap.objects.get_or_create(
            overlap_type=OverlapType.SINGLE_CONTEXT,
            value_type=ClassificationResultValue.ONC_PATH,
            allele=allele_origin_grouping.allele_grouping.allele,
            testing_context=allele_origin_grouping.testing_context_bucket,
            tumor_type_category=allele_origin_grouping.tumor_type_category
        )
        context_overlaps.append(onc_path_context_overlap)
        if allele_origin_grouping.allele_origin_bucket == AlleleOriginBucket.SOMATIC:
            clinsig_context_overlap, _ = Overlap.objects.get_or_create(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE,
                allele=allele_origin_grouping.allele_grouping.allele,
                testing_context=allele_origin_grouping.testing_context_bucket,
                tumor_type_category=allele_origin_grouping.tumor_type_category
            )
            context_overlaps.append(clinsig_context_overlap)

        # TODO optimise as a batch
        for overlap in context_overlaps:
            for class_group in allele_origin_grouping.classificationgrouping_set.prefetch_related('classificationgroupingoverlapcontribution_set').all():
                # need to make sure it is link to all Overlaps that it should be linked to
                # Overlap could be new, or the Overlap could have already existed but the classification grouping is new
                class_group.classificationgroupingoverlapcontribution_set.get_or_create(classification_grouping=class_group, overlap=overlap)

        # now that everything is linked re-calculate the overlap
        for overlap in context_overlaps:
            self.calculate_overlap(overlap)

    def calculate_overlap(self, overlap: Overlap):
        for contribution in overlap.classificationgroupingoverlapcontribution_set.select_related('classification_grouping').all():
            # ensure triage of the required value type is created
            contribution.classification_grouping.classificationgroupingvaluetriage_set.get_or_create(
                classification_grouping=contribution,
                result_value_type=overlap.value_type
            )

        if overlap.value_type == ClassificationResultValue.ONC_PATH:
            OverlapCalculatorOncPath().calculate_and_apply_overlaps_for_ao(overlap)
        elif overlap.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            OverlapCalculatorClinSig().calculate_and_apply_overlaps_for_ao(overlap)

    def calculate_and_apply_overlaps_for(self, classification_grouping: ClassificationGrouping):
        pass
