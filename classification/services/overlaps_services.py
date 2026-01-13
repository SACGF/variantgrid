from classification.enums import TestingContextBucket, OverlapStatus
from classification.models import ClassificationGrouping, ClassificationResultValue, OverlapContributionStatus, \
    OverlapContribution, OverlapEntrySourceTextChoices, Overlap, OverlapType
from classification.services.overlap_calculator import calculator_for_value_type
from snpdb.models import Allele


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
    ONC_PATH_MULTI_CONTEXT = [
        {TestingContextBucket.GERMLINE, TestingContextBucket.NON_CANCER},
        {TestingContextBucket.GERMLINE, TestingContextBucket.HAEMATOLOGY, TestingContextBucket.SOLID_TUMOR}
    ]


    def calculate_and_apply_overlaps_for_allele(self, allele: Allele):
        print("FIXME: NOT RE-SUPPORTED YET")

    @staticmethod
    def update_classification_grouping_overlap_contribution(classification_grouping: ClassificationGrouping):
        value_types: list[ClassificationResultValue] = [ClassificationResultValue.ONC_PATH]
        if classification_grouping.testing_context != TestingContextBucket.GERMLINE:
            value_types.append(ClassificationResultValue.CLINICAL_SIGNIFICANCE)

        for value_type in value_types:
            calc = calculator_for_value_type(value_type)
            value = calc.value_from_summary(classification_grouping.latest_cached_summary)
            is_comparable = calc.is_comparable_value(value)
            is_shared = classification_grouping.share_level_obj.is_discordant_level

            contribution: OverlapContributionStatus
            if value is None:
                contribution = OverlapContributionStatus.NO_VALUE
            elif not is_shared:
                contribution = OverlapContributionStatus.NOT_SHARED
            elif not is_comparable:
                contribution = OverlapContributionStatus.NON_COMPARABLE_VALUE
            else:
                contribution = OverlapContributionStatus.CONTRIBUTING

            overlap_contribution, created = OverlapContribution.objects.update_or_create(
                source=OverlapEntrySourceTextChoices.CLASSIFICATION,
                allele=classification_grouping.allele,
                classification_grouping=classification_grouping,
                testing_context_bucket=classification_grouping.testing_context,
                tumor_type_category=classification_grouping.tumor_type,
                value_type=value_type,
                defaults={
                    "value": value,
                    "contribution": contribution,
                    "effective_date": classification_grouping.latest_classification_modification.curated_date
                }
            )
            if created:
                OverlapServices.link_overlap_contribution(overlap_contribution)
                # make sure this is added to or creates the relevant overlaps
                overlap_contribution.refresh_from_db()

            # now update status of any created overlaps or existing linked overlaps
            for overlap in overlap_contribution.overlap_set.all():
                OverlapServices.recalc_overlap(overlap)

    @staticmethod
    def link_overlap_contribution(overlap_contribution: OverlapContribution):
        # get single context overlap
        """
        overlap_type = models.TextField(choices=OverlapType.choices)
        value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
        allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
        testing_contexts = ArrayField(models.TextField(max_length=1, choices=TestingContextBucket.choices), null=True, blank=True)
        tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
        overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_CONTRIBUTIONS.value)
        valid = models.BooleanField(default=False)  # if it's cross context but only has contributions from 1 context, or if it's NO_SUBMITTERS it shouldn't be valid

        # have to cache the values
        cached_values = JSONField(null=True, blank=True)
        contributions = models.ManyToManyField(OverlapContribution)
        """

        single_context_overlap, created = Overlap.objects.get_or_create(
            overlap_type=OverlapType.SINGLE_CONTEXT,
            value_type=overlap_contribution.value_type,
            allele=overlap_contribution.allele,
            testing_contexts=[overlap_contribution.testing_context_bucket],
            tumor_type_category=overlap_contribution.tumor_type_category,
            defaults={
                "overlap_status": OverlapStatus.NO_CONTRIBUTIONS,
                "valid": False
            }
        )
        single_context_overlap.contributions.add(overlap_contribution)

        if overlap_contribution.value_type == ClassificationResultValue.ONC_PATH:
            for cross_context in OverlapServices.ONC_PATH_MULTI_CONTEXT:
                if overlap_contribution.testing_context_bucket in cross_context:
                    cross_context_overlap, created = Overlap.objects.get_or_create(
                        overlap_type=OverlapType.CROSS_CONTEXT,
                        value_type=ClassificationResultValue.ONC_PATH,
                        allele=overlap_contribution.allele,
                        testing_contexts=list(sorted(cross_context)),
                        tumor_type_category=None,
                        defaults={
                            "overlap_status": OverlapStatus.NO_CONTRIBUTIONS,
                            "valid": False
                        }
                    )
                    cross_context_overlap.contributions.add(overlap_contribution)


    @staticmethod
    def recalc_overlap(overlap: Overlap):
        calculator = calculator_for_value_type(overlap.value_type)
        overlap_status = calculator.calculate_entries(
            list(overlap.contributions.all()))
        overlap.overlap_status = overlap_status

        if overlap.overlap_type == OverlapType.SINGLE_CONTEXT:
            overlap.valid = True
        elif overlap.overlap_type == OverlapType.CROSS_CONTEXT:
            # cross contexts need at least 2 different contexts to be considered valid
            valid = False
            testing_contexts = set()
            entry: OverlapContribution
            for entry in overlap.contributions.all():
                if entry.contribution == OverlapContributionStatus.CONTRIBUTING:
                    testing_contexts.add(entry.testing_context_bucket)
                    if len(testing_contexts) >= 2:
                        valid = True
                        break
            overlap.valid = valid
        overlap.save()

        """
        class OverlapContribution(TimeStampedModel):
        source = models.TextField(choices=OverlapEntrySourceTextChoices.choices)
        scv = models.TextField(null=True, blank=True) # could SCV change?
        # lab_id = models.ForeignKey(Lab, on_delete=CASCADE)
        classification_grouping = models.ForeignKey(ClassificationGrouping, null=True, blank=True, on_delete=CASCADE)
        value_type = models.TextField(choices=ClassificationResultValue.choices)
        value = models.TextField(null=True, blank=True)
        # annoying thing about contribution is it takes a little bit of context knowledge to work out
        contribution = models.TextField(choices=OverlapContributionStatus.choices)
        testing_context_bucket = models.TextField(choices=TestingContextBucket.choices)
        tumor_type_category = models.TextField(null=True, blank=True)
        # TODO do we want to keep date type somewhere?
        effective_date = models.DateField(null=True, blank=True)

        :param sender:n  
        :param instance:
        :param kwargs:
        :return:
        """
