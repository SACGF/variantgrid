import operator
from collections import Counter
from dataclasses import dataclass, field
from functools import cached_property, reduce
from typing import Optional, Self, Tuple

import django
from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import PermissionDenied
from django.db.models import CASCADE, TextChoices, SET_NULL, IntegerChoices, Q, QuerySet
from django.urls import reverse
from django_extensions.db.models import TimeStampedModel
from frozendict import frozendict
from more_itertools import last

from classification.enums import AlleleOriginBucket, ShareLevel, SpecialEKeys, TestingContextBucket
from django.db import models, transaction
from classification.models import Classification, ImportedAlleleInfo, EvidenceKeyMap, ClassificationModification, \
    ConditionResolved, ConditionReference
from classification.models.evidence_mixin_summary_cache import ClassificationSummaryCacheDict, \
    ClassificationSummaryCacheDictPathogenicity, ClassificationSummaryCacheDictSomatic
from genes.models import GeneSymbol
from library.utils import strip_json
from snpdb.models import Allele, Lab


classification_grouping_search_term_signal = django.dispatch.Signal()  # args: "grouping", expects iterable of ClassificationGroupingSearchTermStub


# TODO this needs to be moved to Classification
# class ClassificationQualityLevel(TextChoices):
#     STANDARD = "S", "Standard"
#     LEGACY = "L", "Legacy"
#     INCOMPLETE = "I", "Incomplete"


class ClassificationClassificationBucket(TextChoices):
    BENIGN = "B", "Benign"
    VUS = "V", "VUS"
    PATHOGENIC = "P", "Pathogenic"
    ONCOGENIC = "O", "Oncogenic"  # should oncogenic just be an alias for pathogenic?
    OTHER = "X", "Other"
    CONFLICTING = "C", "Conflicting"
    NO_DATA = "N", "No Data"

    @staticmethod
    def bucket_for_classification(value: Optional[str]):
        if value is None:
            return ClassificationClassificationBucket.NO_DATA
        if value.startswith("VUS"):
            return ClassificationClassificationBucket.VUS
        if value in {"B", "LB"}:
            return ClassificationClassificationBucket.BENIGN
        if value in {"P", "LP"}:
            return ClassificationClassificationBucket.PATHOGENIC
        if value in {"O", "LO"}:
            return ClassificationClassificationBucket.ONCOGENIC
        return ClassificationClassificationBucket.OTHER

        # this goes a little against the buckets that we store directly into


def classification_sort_order(clin_sig: str) -> int:
    return EvidenceKeyMap.instance().get(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_indexes.get(clin_sig, 0)


class OverlapStatus(IntegerChoices):
    NO_SHARED_RECORDS = 0, "No Shared Records"
    SINGLE_SUBMITTER = 10, "Single Shared Submitter"
    NOT_COMPARABLE_OVERLAP = 20, "Multiple Submitters"  # e.g., no method to work out discordance
    AGREEMENT = 30, "Agreement"
    CONFIDENCE = 40, "Confidence"
    DISCORDANCE = 50, "Discordance"
    DISCORDANCE_MEDICALLY_SIGNIFICANT = 60, "Discordance"


class AlleleGrouping(TimeStampedModel):
    allele = models.OneToOneField(Allele, on_delete=models.CASCADE)
    # TODO probably some more summary fields ew could have here?
    # otherwise what is this serving that Allele isn't?

    def get_absolute_url(self) -> str:
        return reverse('allele_grouping_detail', kwargs={"allele_grouping_id": self.pk})

    @cached_property
    def allele_origin_dict(self) -> dict[AlleleOriginBucket, 'AlleleOriginGrouping']:
        by_bucket = {}
        for ao in AlleleOriginGrouping.objects.filter(allele_grouping=self).prefetch_related("classificationgrouping_set"):
            by_bucket[ao.allele_origin_bucket] = ao
        return by_bucket

    def allele_origin_grouping(self, allele_origin_bucket: AlleleOriginBucket) -> Optional['AlleleOriginGrouping']:
        return self.allele_origin_dict.get(allele_origin_bucket)


class AlleleOriginGrouping(TimeStampedModel):
    allele_grouping = models.ForeignKey(AlleleGrouping, on_delete=models.CASCADE)
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices, default=AlleleOriginBucket.UNKNOWN)
    testing_context_bucket = models.CharField(max_length=1, choices=TestingContextBucket.choices, default=TestingContextBucket.UNKNOWN)
    tumor_type_category = models.TextField(null=True, blank=True)

    # class Meta:
    #     unique_together = ("lab", "allele_origin_grouping", "allele_origin_bucket", "testing_context")
    #     indexes = [
    #         models.Index(fields=["allele_origin_grouping"]),
    #         models.Index(fields=["testing_context_bucket", "tumor_type_category"])
    #     ]

    @property
    def allele_origin_bucket_obj(self):
        return AlleleOriginBucket(self.allele_origin_bucket)

    class Meta:
        unique_together = ("allele_grouping", "allele_origin_bucket", "testing_context_bucket", "tumor_type_category")

    overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_SHARED_RECORDS)
    dirty = models.BooleanField(default=True)
    # classification_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)
    # somatic_clinical_significance_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)

    def __lt__(self, other: Self):
        if id_diff := self.allele_grouping.pk - other.allele_grouping.pk:
            return id_diff
        return self.allele_origin_bucket < other.allele_origin_bucket

    def get_absolute_url(self) -> str:
        return reverse('allele_grouping_detail', kwargs={"allele_grouping_id": self.allele_grouping_id})

    def update(self):
        # FIX-ME redo this once we work out how we handle discordance within a single lab ClassificationGrouping

        # # not sure if this is the best solution, or if an AlleleOriginGrouping should just refuse to update if attached allele groupings are dirty
        # all_classification_values = set()
        # all_somatic_clinical_significance_values = set()
        # classification_groupings: list[ClassificationGrouping] = list(self.classificationgrouping_set.all())
        # for cg in classification_groupings:
        #     if cg.dirty:
        #         cg.update()
        # shared_groupings = [cg for cg in classification_groupings if cg.share_level_obj.is_discordant_level]
        # for cg in shared_groupings:
        #     if classification_values := cg.classification_values:
        #         all_classification_values |= set(classification_values)
        #     if somatic_values := cg.somatic_clinical_significance_values:
        #         all_somatic_clinical_significance_values |= set(somatic_values)
        #
        # self.classification_values = list(EvidenceKeyMap.instance().get(SpecialEKeys.CLINICAL_SIGNIFICANCE).sort_values(all_classification_values))
        # self.somatic_clinical_significance_values = [sg.as_str for sg in sorted(SomaticClinicalSignificanceValue.from_str(sg) for sg in all_somatic_clinical_significance_values)]
        #
        # overlap_status: OverlapStatus
        # if len(shared_groupings) == 0:
        #     overlap_status = OverlapStatus.NO_SHARED_RECORDS
        # elif len(shared_groupings) == 1:
        #     overlap_status = OverlapStatus.SINGLE_SUBMITTER
        # elif self.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
        #     overlap_status = OverlapStatus.NOT_COMPARABLE_OVERLAP
        # else:
        #     bucket_mapping = EvidenceKeyMap.instance().get(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("bucket")
        #     buckets = {bucket_mapping.get(class_value) for class_value in self.classification_values}
        #     if None in buckets:
        #         buckets.remove(None)
        #
        #     if len(buckets) > 1:
        #         # discordant
        #         if "P" in all_classification_values or "LP" in all_classification_values:
        #             overlap_status = OverlapStatus.DISCORDANCE_MEDICALLY_SIGNIFICANT
        #         else:
        #             overlap_status = OverlapStatus.DISCORDANCE
        #     else:
        #         if len(all_classification_values) > 1:
        #             overlap_status = OverlapStatus.CONFIDENCE
        #         else:
        #             # complete agreement
        #             overlap_status = OverlapStatus.AGREEMENT
        #
        # self.overlap_status = overlap_status
        self.dirty = False
        self.save()

        # TODO manage pending classifications


# @dataclass
# class ClassificationSubGrouping:
#     latest_modification: ClassificationModification
#     count: int
#
#     @staticmethod
#     def from_modifications(modifications: list[ClassificationModification]) -> list['ClassificationSubGrouping']:
#         # subgroup by classification value for germline and
#         by_status: dict[str, ClassificationSubGrouping] = {}
#         for mod in modifications:
#             key = str(mod.somatic_clinical_significance_value) + str(mod.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
#             if existing := by_status.get(key):
#                 existing.count += 1
#                 existing.latest_modification = max(mod, existing.latest_modification, key=lambda m: (m.curated_date_check, m.pk))
#             else:
#                 by_status[key] = ClassificationSubGrouping(latest_modification=mod, count=1)
#         # FIXME sort
#         return list(by_status.values())


class ClassificationGroupingPathogenicDifference(IntegerChoices):
    NO_DIFF = 0, "No Differences"
    SMALL_DIFF = 1, "Differences"
    CLIN_SIG_DIFFS = 2, "Clinical Significance Difference"


class ClassificationGroupingSomaticDifference(IntegerChoices):
    NO_DIFF = 0, "No Differences"
    AMP_DIFF = 1, "Same Tier, different AMP Level"
    TIER_DIFF = 2, "Different Tier"


class ClassificationGrouping(TimeStampedModel):
    # key
    allele_origin_grouping = models.ForeignKey(AlleleOriginGrouping, on_delete=models.CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    # FIXME remove allele_origin_bucket from ClassificationGrouping, redundant to allele_origin_grouping
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices)
    share_level = models.CharField(max_length=16, choices=ShareLevel.choices())

    @property
    def share_level_obj(self):
        return ShareLevel(self.share_level)

    def can_view(self, user: User) -> bool:
        # inefficient but doesn't duplicate code
        return ClassificationGrouping.filter_for_user(user, ClassificationGrouping.objects.filter(pk=self.pk)).exists()

    def check_can_view(self, user: User):
        if not self.can_view(user):
            raise PermissionDenied("You do not have permission to view this classification grouping.")

    @staticmethod
    def filter_for_user(user: User, qs: QuerySet['ClassificationGrouping']) -> QuerySet['ClassificationGrouping']:
        # TODO, consider making the groups GuardianPermission rather than this manual security check
        if not user.is_superuser:
            permission_q: list[Q] = []
            labs = Lab.valid_labs_qs(user, admin_check=True).select_related("organization")
            orgs = {lab.organization for lab in labs}
            permission_q.append(Q(share_level=ShareLevel.LAB) & Q(lab__in=labs))
            permission_q.append(Q(share_level=ShareLevel.INSTITUTION) & Q(lab__organization__in=orgs))
            permission_q.append(Q(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS))
            return qs.filter(reduce(operator.or_, permission_q))
        else:
            return qs

    classification_count = models.IntegerField(default=0)
    pathogenic_difference = models.IntegerField(choices=ClassificationGroupingPathogenicDifference.choices, default=ClassificationGroupingPathogenicDifference.NO_DIFF)
    somatic_difference = models.IntegerField(choices=ClassificationGroupingSomaticDifference.choices, default=ClassificationGroupingSomaticDifference.NO_DIFF)

    dirty = models.BooleanField(default=True)

    def dirty_up(self):
        self.dirty = True
        self.save(update_fields=["dirty"])
        self.allele_origin_grouping.dirty = True
        self.allele_origin_grouping.save(update_fields=["dirty"])

    # def sub_groupings(self) -> list[ClassificationSubGrouping]:
    #     return ClassificationSubGrouping.from_modifications(self.classification_modifications)

    conditions = models.JSONField(null=True, blank=True)

    zygosity_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)
    # # for discordances
    # classification_bucket
    #
    # # search details
    # gene_symbols
    # conditions
    # curated
    #
    # # cached details
    # summary_curated
    # summary_classification
    # summary_somatic_clinical_significance
    # summary_c_hgvses
    # summary_criteria
    latest_classification_modification = models.ForeignKey(ClassificationModification, on_delete=SET_NULL, null=True, blank=True)
    latest_allele_info = models.ForeignKey(ImportedAlleleInfo, on_delete=SET_NULL, null=True, blank=True)

    def __lt__(self, other):
        def _sort_value(obj: ClassificationGrouping):
            return obj.share_level_obj, obj.lab
        return _sort_value(self) < _sort_value(other)

    def get_absolute_url(self):
        return reverse('classification_grouping_detail', kwargs={"classification_grouping_id": self.pk})

    @staticmethod
    def _desired_grouping_for_classification(classification: Classification) -> Tuple[Optional['ClassificationGrouping'], bool]:
        # withdrawn classifications are removed from groupings
        if classification.withdrawn:
            return None, False

        allele = classification.allele_object
        lab = classification.lab
        share_level = classification.share_level
        if allele:
            allele_grouping, _ = AlleleGrouping.objects.get_or_create(allele=allele)
            allele_origin_grouping, _ = AlleleOriginGrouping.objects.get_or_create(
                allele_grouping=allele_grouping,
                allele_origin_bucket=classification.allele_origin_bucket,
                testing_context_bucket=classification.testing_context_bucket,
                tumor_type_category=classification.tumor_type_category
            )
            allele_origin_grouping.dirty = True
            allele_origin_grouping.save(update_fields=["dirty"])

            grouping, is_new = ClassificationGrouping.objects.get_or_create(
                allele_origin_grouping=allele_origin_grouping,
                lab=lab,
                allele_origin_bucket=classification.allele_origin_bucket,
                share_level=share_level
            )
            return grouping, is_new
        return None, False

    @staticmethod
    def assign_grouping_for_classification(classification: Classification, force_dirty_up=True) -> bool:
        """
        :param classification: The classification that needs to go into a grouping
        :param force_dirty_up: If grouping for the classification needs to be marked as dirty even if the classification is not changing groupings
        :return: A boolean indicating if the classification changed groupings
        """
        desired_grouping, is_new_grouping = ClassificationGrouping._desired_grouping_for_classification(classification)
        if desired_grouping:
            entry, is_new_entry = ClassificationGroupingEntry.objects.get_or_create(
                classification=classification,
                defaults={"grouping": desired_grouping}
            )
            if is_new_grouping and is_new_entry:
                # if we've got the first record in the grouping, process it right now, so we can see it during the import process
                desired_grouping.update()
                return True
            elif is_new_entry:
                entry.dirty_up()
                return True
            else:
                # entry previous existed, make sure it's pointing to the correct grouping
                old_grouping = entry.grouping
                if entry.grouping != desired_grouping:
                    entry.grouping = desired_grouping
                    entry.save()

                    old_grouping.dirty_up()
                    desired_grouping.dirty_up()
                    return True
                elif force_dirty_up:
                    entry.dirty_up()
                return False
        else:
            # if we don't even have an allele, make sure we are removed from any grouping
            if existing := ClassificationGroupingEntry.objects.filter(classification=classification).first():
                grouping = existing.grouping
                grouping.dirty_up()
                existing.delete()
                return True
            return False

    @cached_property
    def classification_modifications(self) -> list[ClassificationModification]:
        # show in date order
        all_classifications = self.classificationgroupingentry_set.values_list("classification", flat=True)
        all_modifications = ClassificationModification.objects.filter(classification_id__in=all_classifications, is_last_published=True)
        all_modifications = all_modifications.select_related("classification")
        return list(sorted(all_modifications, key=lambda mod: mod.curated_date_check))

    @cached_property
    def allele(self) -> Allele:
        return self.allele_origin_grouping.allele_grouping.allele

    # def to_json(self):
    #     scs = self.latest_classification_modification.somatic_clinical_significance_value
    #     return {
    #         "lab": str(self.lab),
    #         "classification": self.latest_classification_modification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
    #         "somatic_clinical_significance": scs.as_json() if scs else None,
    #         "curated_date": str(self.latest_curation_date)  # fixme, need a quick way for yyyy-mm-dd
    #     }

    @transaction.atomic
    def update(self):

        self.allele_origin_grouping.dirty = True
        self.allele_origin_grouping.save(update_fields=["dirty"])

        if all_modifications := self.classification_modifications:
            # FIND THE MOST RECENT CLASSIFICATION
            best_classification = last(all_modifications)

            self.latest_classification_modification = best_classification
            self.latest_allele_info = best_classification.classification.allele_info

            # TODO check for dirty values
            all_terms = Counter()
            all_free_text_conditions = Counter()

            all_terms.keys()

            # UPDATE CLASSIFICATION / CLINICAL SIGNIFICANCE
            # TODO - CLINICAL SIGNIFICANCE
            all_zygosities: set[str] = set()

            all_buckets = set()
            all_pathogenic_values = set()
            all_tiers = set()
            all_levels = set()
            condition_references: list[ConditionReference] = []

            for modification in all_modifications:
                if condition := modification.classification.condition_resolution_obj:
                    for term in condition.terms:
                        all_terms[term] += 1
                    # all_terms |= set(condition.terms)
                    # if plain_text := condition.plain_text:
                    #     all_free_text_conditions.add(plain_text)
                elif condition_text := modification.get(SpecialEKeys.CONDITION):
                    all_free_text_conditions[condition_text] += 1
                    # all_free_text_conditions.add(condition_text)

                all_zygosities |= set(modification.get_value_list(SpecialEKeys.ZYGOSITY))

                # only store valid terms as quick links to the classification
                # all_terms = {term for term in all_terms if term.is_valid_for_condition}
                # self._update_conditions(all_terms)
                summary_dict: ClassificationSummaryCacheDict = modification.classification.summary

                pathogenicity_dict: ClassificationSummaryCacheDictPathogenicity = summary_dict.get("pathogenicity", {})
                somatic_dict: ClassificationSummaryCacheDictSomatic = summary_dict.get("somatic", {})

                if bucket := pathogenicity_dict.get("bucket"):
                    all_buckets.add(bucket)
                if path_val := pathogenicity_dict.get("classification"):
                    all_pathogenic_values.add(path_val)
                if tier := somatic_dict.get("clinical_significance"):
                    all_tiers.add(tier)
                if level := somatic_dict.get("amp_level"):
                    all_levels.add(level)

            evidence_map = EvidenceKeyMap.instance()

            # the below shouldn't happen but has in development environemnts
            if None in all_zygosities:
                all_zygosities.remove(None)

            all_zygosities = list(sorted(evidence_map[SpecialEKeys.ZYGOSITY].sort_values(all_zygosities)))
            self.zygosity_values = all_zygosities

            # self.conditions = strip_json(ConditionResolved.from_uncounted_terms(
            #     terms=list(sorted(all_terms)),
            #     plain_text_terms=list(sorted(all_free_text_conditions))
            # ).to_json(include_join=False))
            for term, count in all_terms.items():
                condition_references.append(ConditionReference(term=term, count=count))
            for text, count in all_free_text_conditions.items():
                condition_references.append(ConditionReference(name=text, count=count))
            condition_references.sort()
            self.conditions = strip_json(ConditionResolved(references=condition_references).to_json())

            self.classification_count = len(all_modifications)

            pathogenic_difference = ClassificationGroupingPathogenicDifference.NO_DIFF
            if len(all_buckets) > 1:
                pathogenic_difference = ClassificationGroupingPathogenicDifference.CLIN_SIG_DIFFS
            elif len(all_pathogenic_values) > 1:
                pathogenic_difference = ClassificationGroupingPathogenicDifference.SMALL_DIFF

            somatic_difference = ClassificationGroupingSomaticDifference.NO_DIFF
            if len(all_tiers) > 1:
                somatic_difference = ClassificationGroupingSomaticDifference.TIER_DIFF
            elif len(all_levels) > 1:
                somatic_difference = ClassificationGroupingSomaticDifference.AMP_DIFF

            self.pathogenic_difference = pathogenic_difference
            self.somatic_difference = somatic_difference

            all_term_stubs: list[ClassificationGroupingSearchTermStub] = []
            for _, term_stubs in classification_grouping_search_term_signal.send(sender=ClassificationGrouping, grouping=self):
                if term_stubs:
                    all_term_stubs += [ts.to_search_term(grouping=self) for ts in term_stubs]
            ClassificationGroupingSearchTerm.objects.filter(grouping=self).delete()
            ClassificationGroupingSearchTerm.objects.bulk_create(all_term_stubs)

            self.dirty = False
            self.save()
        else:
            # there are no classifications, time to die
            self.delete()

    def gene_symbols(self):
        terms = set(self.classificationgroupingsearchterm_set.filter(term_type=ClassificationGroupingSearchTermType.GENE_SYMBOL).values_list("term", flat=True))
        return GeneSymbol.objects.filter(symbol__in=terms)

    @staticmethod
    def update_all_dirty():
        # maybe move this out since it does AlleleOriginGroupings too
        for dirty in ClassificationGrouping.objects.filter(dirty=True).iterator():
            dirty.update()
        for dirty in AlleleOriginGrouping.objects.filter(dirty=True).iterator():
            dirty.update()


class ClassificationGroupingSearchTermType(TextChoices):
    # CLASSIFICATION_ID = "CR_ID", "Classification Record ID"
    # CLASSIFICATION_LAB_ID = "CR_LAB_ID", "Classification Lab Record ID"
    # SOURCE_ID = "SOURCE_ID", "Classification Source ID"
    CONDITION_ID = "CON_ID", "Condition ID"
    CLINVAR_SCV = "SCV", "Clinvar SCV"
    GENE_SYMBOL = "GENE_SYMBOL", "Gene Symbol"
    DISCORDANCE_REPORT = "DR", "Discordance Report"

    @property
    def is_partial_text(self):
        return False


class ClassificationGroupingSearchTerm(TimeStampedModel):
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)
    term = models.TextField()
    term_type = models.CharField(max_length=20, choices=ClassificationGroupingSearchTermType.choices)
    extra = models.JSONField(null=True, blank=True)

    class Meta:
        unique_together = ("grouping", "term", "term_type")
        indexes = [models.Index(fields=["term_type", "term"])]

    @staticmethod
    def filter_q(term_type: ClassificationGroupingSearchTermType, terms: str | list[str]) -> Q:
        # TODO add support for partial search based on type
        # TODO add support for different privacy levels based on type
        if isinstance(terms, str):
            terms = [terms]
        if term_type.is_partial_text:
            matching_values = ClassificationGroupingSearchTerm.objects.filter(
                Q(term_type=term_type) & reduce(operator.or_, [Q(term__icontains=t.upper()) for t in terms])
            ).values_list('grouping_id', flat=True)
        else:
            matching_values = ClassificationGroupingSearchTerm.objects.filter(term_type=term_type, term__in=[t.upper() for t in terms]).values_list('grouping_id', flat=True)
        return Q(pk__in=matching_values)


@dataclass(frozen=True)
class ClassificationGroupingSearchTermStub:
    term_type: ClassificationGroupingSearchTermType
    term: str
    extra: Optional[frozendict] = None

    def to_search_term(self, grouping: ClassificationGrouping) -> 'ClassificationGroupingSearchTerm':
        return ClassificationGroupingSearchTerm(
            grouping=grouping,
            term=self.term,
            term_type=self.term_type,
            extra=self.extra
        )


@dataclass
class ClassificationGroupingSearchTermBuilder:
    term: str
    term_type: ClassificationGroupingSearchTermType
    extra: dict = field(default_factory=dict)

    def as_stub(self) -> ClassificationGroupingSearchTermStub:
        return ClassificationGroupingSearchTermStub(
            term=self.term,
            term_type=ClassificationGroupingSearchTermType(self.term_type),
            extra=frozendict(self.extra) if self.extra else None
        )


class ClassificationGroupingEntry(TimeStampedModel):
    # this is just here so this model can stay completely separate from Classification
    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)

    def dirty_up(self):
        self.grouping.dirty_up()

    def __str__(self):
        return f"Grouping {self.grouping} : Classification {self.classification}"

    class Meta:
        unique_together = ("grouping", "classification")
