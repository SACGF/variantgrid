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
classification_grouping_update_signal = django.dispatch.Signal()  # TODO is there a point for this over ClassificationGrouping
#classification_grouping_onc_path_signal = django.dispatch.Signal()  # args: "instance", expects Classification
#classification_grouping_clin_sig_signal = django.dispatch.Signal()  # args: "instance", expects Classification

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

#
# class OverlapStatus(IntegerChoices):
#     NO_SHARED_RECORDS = 0, "No Shared Records"
#     SINGLE_SUBMITTER = 10, "Single Shared Submitter"
#     NOT_COMPARABLE_OVERLAP = 20, "Multiple Submitters"  # e.g., no method to work out discordance
#     AGREEMENT = 30, "Agreement"
#     CONFIDENCE = 40, "Confidence"
#     DISCORDANCE = 50, "Discordance"
#     DISCORDANCE_MEDICALLY_SIGNIFICANT = 60, "Discordance"


class AlleleGrouping(TimeStampedModel):
    allele = models.OneToOneField(Allele, on_delete=models.CASCADE)
    # TODO probably some more summary fields ew could have here?
    # otherwise what is this serving that Allele isn't?

    def __str__(self):
        return f"({self.pk}) Allele Grouping for {self.allele}"

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

    @property
    def testing_context_bucket_obj(self):
        return TestingContextBucket(self.testing_context_bucket)

    def __str__(self):
        return f"{self.allele_grouping.allele} {self.get_allele_origin_bucket_display()} Testing Context: {self.get_testing_context_bucket_display()} Sub-Type: {self.tumor_type_category}"

    def labels(self, include_allele_origin=True) -> list[Optional[str]]:
        parts = []
        if include_allele_origin:
            parts.append(AlleleOriginBucket(self.allele_origin_bucket).label)
        if self.allele_origin_bucket == AlleleOriginBucket.SOMATIC:
            if self.testing_context_bucket in {TestingContextBucket.UNKNOWN.value, TestingContextBucket.OTHER.value}:
                parts.append("Testing Context")
            parts.append(TestingContextBucket(self.testing_context_bucket).label)
            if TestingContextBucket(self.testing_context_bucket).should_have_subdivide:
                parts.append(self.tumor_type_category)
        return parts

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

    # overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_SHARED_RECORDS)
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
        # FIXME as overlaps now belong in
        # from classification.services.overlaps_services import OverlapServices
        # OverlapServices().calculate_and_apply_overlaps_for_ao(self)
        # self.dirty = False
        # self.save()
        print("FIXME: AlleleOriginGrouping.Update currently doesn't do anything")


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
    share_level = models.CharField(max_length=16, choices=ShareLevel.choices())
    classification_count = models.IntegerField(default=0)
    # these differences are internal within a lab, even if it's medical significant the overall
    pathogenic_difference = models.IntegerField(choices=ClassificationGroupingPathogenicDifference.choices, default=ClassificationGroupingPathogenicDifference.NO_DIFF)
    somatic_difference = models.IntegerField(choices=ClassificationGroupingSomaticDifference.choices, default=ClassificationGroupingSomaticDifference.NO_DIFF)
    dirty = models.BooleanField(default=True)
    conditions = models.JSONField(null=True, blank=True)
    zygosity_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)
    latest_classification_modification = models.ForeignKey(ClassificationModification, on_delete=SET_NULL, null=True, blank=True)
    latest_cached_summary = models.JSONField(null=False, blank=True, default=dict)
    latest_allele_info = models.ForeignKey(ImportedAlleleInfo, on_delete=SET_NULL, null=True, blank=True)

    def __str__(self):
        parts = [
            f"({self.pk})",
            "Classification-Grouping",
            f"{self.allele_origin_grouping.allele_grouping.allele:CA}",
            self.allele_origin_grouping.get_testing_context_bucket_display(),
            str(self.lab),
            "Shared" if self.share_level_obj.is_discordant_level else "Not-shared",
        ]
        if classification := self.latest_cached_summary.get("pathogenicity").get("classification"):
            parts.append(f"Class({classification})")
        if clin_sig := self.latest_cached_summary.get("somatic").get("clin_sig"):
            parts.append(f"ClinSig({clin_sig})")
        return " ".join(parts)

    @property
    def category_text_compact(self) -> str:
        parts = self.allele_origin_grouping.labels(include_allele_origin=True)
        if parts[0] == "Somatic":
            parts = parts[1:]
        return " - ".join(p for p in parts if p)

    @property
    def testing_context(self) -> TestingContextBucket:
        return TestingContextBucket(self.allele_origin_grouping.testing_context_bucket)

    @property
    def tumor_type(self) -> str:
        return self.allele_origin_grouping.tumor_type_category

    @property
    def allele_origin_bucket(self):
        return self.allele_origin_grouping.allele_origin_bucket

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

    def dirty_up(self):
        self.dirty = True
        self.save(update_fields=["dirty"])
        self.allele_origin_grouping.dirty = True
        self.allele_origin_grouping.save(update_fields=["dirty"])

    # def sub_groupings(self) -> list[ClassificationSubGrouping]:
    #     return ClassificationSubGrouping.from_modifications(self.classification_modifications)

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

        # FIXME there's a bunch of overlap calculation here that doesn't need to happen
        # as that is now managed by Conflicts

        self.allele_origin_grouping.dirty = True
        self.allele_origin_grouping.save(update_fields=["dirty"])

        if all_modifications := self.classification_modifications:
            # FIND THE MOST RECENT CLASSIFICATION
            best_classification = last(all_modifications)

            primary_onc_path_change = False
            primary_clin_sig_change = False
            previous_summary: ClassificationSummaryCacheDict
            if previous_summary := self.latest_cached_summary:
                new_summary: ClassificationSummaryCacheDict = best_classification.classification.summary
                if previous_summary.get("pathogenicity", {}).get("classification") != new_summary.get("pathogenicity", {}).get("classification"):
                    primary_onc_path_change = True
                if previous_summary.get("somatic", {}).get("clinical_significance", {}) != new_summary.get("somatic", {}).get("clinical_significance"):
                    primary_clin_sig_change = True

            self.latest_classification_modification = best_classification
            # TODO now see if there's a change in overall classification / clinical significance - for the grouping
            self.latest_cached_summary = best_classification.classification.summary
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
            classification_grouping_update_signal.send(sender=ClassificationGrouping, instance=self)
            # if primary_clin_sig_change:
            #     classification_grouping_clin_sig_signal.send(sender=ClassificationGrouping, instance=self)
            # if primary_onc_path_change:
            #     classification_grouping_onc_path_signal.send(sender=ClassificationGrouping, instance=self)
        else:
            # there are no classifications, time to die
            # FIXME add soft delete
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
    # should only ever be a single grouping per classification (should classification be made unique??)
    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)

    def dirty_up(self):
        self.grouping.dirty_up()

    def __str__(self):
        return f"Grouping {self.grouping} : Classification {self.classification}"

    class Meta:
        unique_together = ("grouping", "classification")

    @staticmethod
    def grouping_for(classification: Classification) -> Optional[ClassificationGrouping]:
        if entry := ClassificationGroupingEntry.objects.filter(classification=classification).select_related("grouping").first():
            return entry.grouping
        return None


# @dataclass(frozen=True)
# class ConflictKey:
#     conflict_type: ConflictType
#     allele_origin_bucket: Optional[AlleleOriginBucket]
#     testing_context_bucket: Optional[TestingContextBucket] = None
#     tumor_type_category: Optional[str] = None
#
#     severity: Optional[OverlapStatus] = None  # just used when we're determining if we can merge 2 conflicts
#
#
# class ConflictQuerySet(QuerySet['Conflict']):
#     def for_lab(self, lab: Lab):
#         return self.filter(
#             Exists(ConflictLab.objects.filter(conflict=OuterRef('pk'), lab=lab))
#         )
#
#     def for_labs(self, labs: Union[Iterable[int], Iterable[Lab]]) -> Self:
#         lab_ids = labs
#
#         # when_resolved = When(
#         #     Q(Exists(Subquery(ConflictHistory.objects.filter(
#         #         conflict=OuterRef('pk'),
#         #         severity__gte=ConflictSeverity.MEDIUM,
#         #     )))) &
#         #     Q(severity__lt=ConflictSeverity.MEDIUM),
#         #     then=Value("X")
#         # )
#
#         when_resolved = When(
#             Q(Exists(Subquery(ConflictLab.objects.filter(
#                 conflict=OuterRef('pk'),
#                 lab__in=lab_ids
#             )))) &
#             Q(severity__lt=ConflictSeverity.MEDIUM) &
#             Q(Exists(Subquery(ConflictHistory.objects.filter(
#                 conflict=OuterRef('pk'),
#                 severity__gte=ConflictSeverity.MEDIUM
#             )))),
#             then=Value(DiscordanceReportNextStep.RESOLVED)
#         )
#
#         when_not_active = When(
#             ~Q(Exists(Subquery(ConflictLab.objects.filter(
#                 conflict=OuterRef('pk'),
#                 active=True,
#                 lab__in=lab_ids
#             )))),
#             then=Value(DiscordanceReportNextStep.NOT_INVOLVED)
#         )
#
#         when_no_conflict = When(
#             severity__lt=ConflictSeverity.MEDIUM,
#             then=Value(DiscordanceReportNextStep.NO_CONFLICT)
#         )
#
#         when_waiting_on_your_triage = When(
#             Exists(Subquery(ConflictLab.objects.filter(
#                 conflict=OuterRef('pk'),
#                 status=DiscordanceReportTriageStatus.PENDING,
#                 active=True,
#                 lab__in=lab_ids
#             ))),
#             then=Value(DiscordanceReportNextStep.AWAITING_YOUR_TRIAGE)
#         )
#
#         when_waiting_on_your_amend = When(
#             Exists(Subquery(ConflictLab.objects.filter(
#                 conflict=OuterRef('pk'),
#                 status=DiscordanceReportTriageStatus.REVIEWED_WILL_FIX,
#                 active=True,
#                 lab__in=lab_ids
#             ))),
#             then=Value(DiscordanceReportNextStep.AWAITING_YOUR_AMEND)
#         )
#
#         when_waiting_on_other_lab = When(
#             Exists(Subquery(ConflictLab.objects.filter(
#                 conflict=OuterRef('pk'),
#                 status__in={DiscordanceReportTriageStatus.REVIEWED_WILL_FIX, DiscordanceReportTriageStatus.PENDING},
#                 active=True
#             ).exclude(lab__in=lab_ids))),
#             then=Value(DiscordanceReportNextStep.AWAITING_OTHER_LAB)
#         )
#
#         when_complex = When(
#             ~Q(
#                 Exists(Subquery(ConflictLab.objects.filter(
#                     conflict=OuterRef('pk'),
#                     active=True
#                 ).exclude(status=DiscordanceReportTriageStatus.COMPLEX)))
#             ), then=Value(DiscordanceReportNextStep.UNANIMOUSLY_COMPLEX)
#         )
#
#         qs = self.annotate(overall_status=Case(
#             when_not_active,
#             when_resolved,
#             when_no_conflict,
#             when_waiting_on_your_triage,
#             when_waiting_on_your_amend,
#             when_waiting_on_other_lab,
#             when_complex,
#             default=Value(DiscordanceReportNextStep.TO_DISCUSS)))
#
#         return qs.filter(
#             Exists(ConflictLab.objects.filter(conflict=OuterRef('pk'), lab__in=labs))
#         )
#
#
# class ConflictObjectManager(Manager):
#
#     def get_queryset(self):
#         qs = ConflictQuerySet(self.model, using=self._db)
#         qs = qs.annotate(severity=Subquery(ConflictHistory.objects.filter(conflict_id=OuterRef("pk"), is_latest=True).values('severity')))
#         qs = qs.annotate(change_date=Subquery(ConflictHistory.objects.filter(conflict_id=OuterRef("pk"), is_latest=True).values('created')))
#         qs = qs.annotate(data=Subquery(ConflictHistory.objects.filter(conflict_id=OuterRef("pk"), is_latest=True).values('data')))
#         qs = qs.annotate(past_issue=Exists(Subquery(ConflictHistory.objects.filter(conflict_id=OuterRef("pk"), severity__gte=ConflictSeverity.MEDIUM))))
#         return qs
#
#     # def with_severity(self):
#     #     return self.annotate(severity=Subquery(ConflictHistory.objects.filter(conflict_id=OuterRef("pk"), is_latest=True).values('severity')))
#
#
# class Conflict(ReviewableModelMixin, PreviewModelMixin, TimeStampedModel):
#     objects = ConflictObjectManager()
#     allele = models.ForeignKey(Allele, on_delete=CASCADE)
#     conflict_type = models.CharField(max_length=1, choices=ConflictType.choices)
#     # FIXME while is allele origin bucket and testing context bucket nullable?
#     allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices, null=True, blank=True)
#     testing_context_bucket = models.CharField(max_length=1, choices=TestingContextBucket.choices, null=True, blank=True)
#     tumor_type_category = models.TextField(null=True, blank=True)
#     meta_data = models.JSONField(null=False, blank=False, default=dict)
#
#     def get_absolute_url(self) -> str:
#         return reverse('conflict', kwargs={"conflict_id": self.pk})
#
#     @classmethod
#     def preview_category(cls) -> str:
#         return "Overlap"
#
#     @classmethod
#     def preview_icon(cls) -> str:
#         return "fa-solid fa-arrow-down-up-across-line"
#
#     @property
#     def preview(self) -> 'PreviewData':
#
#         c_hgvs_key_values = []
#         for c_hgvs in self.c_hgvses():
#             c_hgvs_key_values.append(
#                 PreviewKeyValue(key=f"{c_hgvs.genome_build} c.HGVS", value=str(c_hgvs), dedicated_row=True)
#             )
#
#         status_text=ConflictSeverity(self.latest.severity).label
#
#         # note there's also preview_extra_signal that provides the lab data
#         return self.preview_with(
#             identifier=f"CR_{self.pk}",
#             summary_extra=
#                 [PreviewKeyValue(key="Status", value=status_text, dedicated_row=True)] +
#                 [PreviewKeyValue(key="Allele", value=f"{self.allele:CA}", dedicated_row=True)] +
#                 c_hgvs_key_values
#         )
#
#     @property
#     def show_triage(self) -> bool:
#         return self.latest.severity >= ConflictSeverity.MEDIUM
#
#     def get_conflict_type_display(self):
#         return ConflictType(self.conflict_type).label_for_context(AlleleOriginBucket(self.allele_origin_bucket))
#
#     @property
#     def _sort_key(self):
#         # put testing conflict_type bucket last so we group conflicts that are for the same data but different contexts e.g.
#         return (
#             self.allele,
#             self.allele_origin_bucket or AlleleOriginBucket.UNKNOWN,
#             self.testing_context_bucket or TestingContextBucket.UNKNOWN,
#             self.tumor_type_category or "",
#             self.conflict_type,
#         )
#
#     def __lt__(self, other):
#         return self._sort_key < other._sort_key
#
#     def c_hgvses(self) -> list[CHGVS]:
#         if c_hgvs_values := self.meta_data.get("c_hgvs"):
#             return [CHGVS.from_json_short(c_hgvs_value) for c_hgvs_value in c_hgvs_values]
#         return []
#
#     class Meta:
#         unique_together = ("allele", "conflict_type", "allele_origin_bucket", "testing_context_bucket", "tumor_type_category")
#
#     @property
#     def context_summary_short(self) -> str:
#         parts = [self.allele_origin_bucket]
#         if self.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
#             parts.append(self.get_testing_context_bucket_display())
#         if self.tumor_type_category:
#             parts.append(self.tumor_type_category)
#         parts.append(self.get_conflict_type_display())
#         return " ".join(parts)
#
#     @property
#     def context_summary(self) -> str:
#         parts = [self.get_allele_origin_bucket_display()]
#         if self.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
#             parts.append(self.get_testing_context_bucket_display())
#         if self.tumor_type_category:
#             parts.append(self.tumor_type_category)
#         parts.append(self.get_conflict_type_display())
#         return " ".join(parts)
#
#     # @property
#     # def concise_str(self) -> str:
#     #     parts = []
#     #     parts.append(self.allele_origin_bucket) # just a single letter for Allele Origin
#     #     if self.tumor_type_category:
#     #         parts.append(self.tumor_type_category)
#
#     def __str__(self) -> str:
#         parts = [f"{self.allele:CA}", self.get_allele_origin_bucket_display()]
#         if self.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
#             parts.append(self.get_testing_context_bucket_display())
#         if self.tumor_type_category:
#             parts.append(self.tumor_type_category)
#         parts.append(self.get_conflict_type_display())
#         parts.append("-")
#         try:
#             parts.append(self.latest.get_severity_display())
#         except ConflictHistory.DoesNotExist:
#             pass
#
#         return " ".join(parts)
#
#     @cached_property
#     def conflict_labs(self) -> list['ConflictLab']:
#         return list(self.conflictlab_set.order_by('lab__organization__name', 'lab__name'))
#
#     @cached_property
#     def comments(self) -> list['ConflictComment']:
#         return list(self.conflictcomment_set.order_by('-created'))
#
#     @cached_property
#     def latest(self) -> 'ConflictHistory':
#         if the_latest := ConflictHistory.objects.filter(conflict=self, is_latest=True).first():
#             return the_latest
#         raise ConflictHistory.DoesNotExist(f"Conflict {self.pk} has no ConflictHistory marked as is_latest")
#
#     @property
#     def current_severity_as_of(self) -> datetime:
#         severity: Optional[ConflictSeverity] = None
#         severity_date: Optional[datetime] = None
#         for history in self.conflicthistory_set.order_by('-created'):
#             if severity is None:
#                 severity = history.severity
#                 severity_date = history.created
#             else:
#                 if severity == history.severity:
#                     severity_date = history.created
#                 else:
#                     break
#
#         return severity_date
#
#     def history(self, newest_to_oldest: bool = True) -> Iterable['ConflictHistory']:
#         qs = self.conflicthistory_set.all()
#         if newest_to_oldest:
#             qs = qs.order_by('-created')
#         else:
#             qs = qs.order_by('created')
#         return qs
#
#     """
#     allele_grouping = models.ForeignKey(AlleleGrouping, on_delete=models.CASCADE)
#     allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices, default=AlleleOriginBucket.UNKNOWN)
#     testing_context_bucket = models.CharField(max_length=1, choices=TestingContextBucket.choices, default=TestingContextBucket.UNKNOWN)
#     tumor_type_category = models.TextField(null=True, blank=True)
#     """
#
#     def allele_origin_grouping(self) -> Optional[AlleleOriginGrouping]:
#         if allele_grouping := AlleleGrouping.objects.filter(allele=self.allele).get():
#             if allele_origin_grouping := AlleleOriginGrouping.objects.filter(
#                 allele_grouping=allele_grouping,
#                 allele_origin_bucket=self.allele_origin_bucket,
#                 testing_context_bucket=self.testing_context_bucket,
#                 tumor_type_category=self.tumor_type_category
#             ).get():
#                 return allele_origin_grouping
#         return None
#
#     # Methods for review
#
#     @property
#     def reviewing_labs(self) -> list[Lab]:
#         return [cl.lab for cl in ConflictLab.objects.filter(conflict=self, active=True).select_related("lab")]
#
#     def post_review_url(self, review: Review) -> str:
#         return reverse('conflict_review_complete', kwargs={'review_id': review.pk})
#
#     # def grouped_data(self) -> 'ConflictLabGrouped':
#     #     lab_comments: dict[Lab, list[ConflictLabComment]] = defaultdict(list)
#     #     for comment in self.conflictlab_set.select_related("lab").all():
#     #         lab_comments[comment.lab].append(comment)
#     #
#     #     excluded: list[ConflictDataRow] = []
#     #     lab_data: dict[Lab, list[ConflictDataRow]] = defaultdict(list)
#     #     if latest_history := self.conflicthistory_set.select_related("lab").filter(is_latest=True).first():
#     #         for data_row in latest_history.data_rows():
#     #             if data_row.exclude:
#     #                 excluded.append(data_row)
#     #             else:
#     #                 lab_data[data_row.lab].append(data_row)
#     #
#     #     lab_groupings: list[ConflictLabGrouping] = []
#     #     for lab, data_rows in lab_data.items():
#     #         lab_groupings.append(ConflictLabGrouping(
#     #             lab=lab,
#     #             data=data_rows,
#     #             comments=lab_comments[lab]
#     #         ))
#     #     return list(sorted(lab_groupings)), excluded
#
#
# class ConflictHistory(TimeStampedModel):
#     conflict = models.ForeignKey(Conflict, on_delete=CASCADE)
#     data = models.JSONField(null=False, blank=False)
#     severity = models.IntegerField(choices=ConflictSeverity.choices)
#     is_latest = models.BooleanField(default=False)
#
#     class Meta:
#         indexes = [models.Index(fields=["conflict", "is_latest"])]
#
#     @cached_property
#     def involved_lab_ids(self) -> set[int]:
#         lab_ids: set[int] = set()
#         for row in self.data_rows():
#             if not row.exclude:
#                 lab_ids.add(row.lab_id)
#         return lab_ids
#
#     @cached_property
#     def involved_labs(self) -> list[Lab]:
#         return list(sorted(Lab.objects.filter(pk__in=self.involved_lab_ids)))
#
#     @property
#     def date_detected_str(self) -> str:
#         date_str = f"{self.created:%Y-%m-%d}"
#         if (timezone.now() - self.created) <= timedelta(days=1):
#             date_str = f"{date_str} (NEW)"
#         return date_str
#
#     # TODO, caching this and letting other methods annotate it is a bit messy in some places
#     # but a lot cleaner in others
#     def data_rows(self) -> list['ConflictDataRow']:
#         rows = self.data.get("rows")
#         from classification.services.conflict_services import ConflictDataRow
#         # TODO sort
#         return [ConflictDataRow.from_json(row) for row in rows]
#
#     def data_rows_for_user(self, user: User):
#         data_rows = self.data_rows()
#         return [dr for dr in data_rows if dr.can_view(user)]
#
#
# class ConflictLab(TimeStampedModel):
#     conflict = models.ForeignKey(Conflict, on_delete=CASCADE)
#     classification_grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE, null=True)
#     lab = models.ForeignKey(Lab, on_delete=CASCADE)
#     active = models.BooleanField(default=True)  # set to False if lab has withdrawn
#     status = models.TextField(choices=DiscordanceReportTriageStatus.choices, default=DiscordanceReportTriageStatus.PENDING)
#
#     class Meta:
#         unique_together = ("conflict", "lab")
#
#
# class ConflictComment(TimeStampedModel):
#     conflict = models.ForeignKey(Conflict, on_delete=CASCADE)
#     lab = models.ForeignKey(Lab, on_delete=CASCADE, null=True, blank=True)
#     user = models.ForeignKey(User, on_delete=PROTECT)
#     comment = models.TextField(null=False, blank=False)
#     meta_data = models.JSONField(null=False, blank=False, default=dict)
#
#     def __str__(self):
#         return f"{self.user}: {self.comment}"
#
#     @property
#     def meta_data_html(self) -> Optional[str]:
#         if not self.meta_data:
#             return None
#         else:
#             parts = []
#             for lab_id, status in self.meta_data.items():
#                 lab = Lab.objects.get(pk=lab_id)
#                 status_obj = DiscordanceReportTriageStatus(status)
#                 parts.append(f"{html.escape(str(lab))} -> {status_obj.label}")
#             return SafeString("<br/>".join(parts))
#
#
# class ConflictNotificationStatus(TextChoices):
#     QUEUED = "Q", "Queued"
#     PROCESSING = "P", "Processing"
#     COMPLETE = "C", "Complete"
#
#
# class ConflictNotificationRun(TimeStampedModel):
#     status = models.TextField(choices=ConflictNotificationStatus.choices, default=ConflictNotificationStatus.QUEUED)
#     # TODO maybe add some overall stats
#
#
# class ConflictNotification(TimeStampedModel):
#     conflict = models.ForeignKey(Conflict, on_delete=CASCADE)
#     current_state = models.ForeignKey(ConflictHistory, on_delete=CASCADE, related_name='+')
#     previous_state = models.ForeignKey(ConflictHistory, on_delete=CASCADE, related_name='+', null=True, blank=True)
#     notification_run = models.ForeignKey(ConflictNotificationRun, on_delete=SET_NULL, null=True, blank=True)
#
#     def __lt__(self, other):
#         return self.conflict < other.conflict