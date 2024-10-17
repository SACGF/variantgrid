import operator
from dataclasses import dataclass, field
from functools import cached_property, reduce
from typing import Optional, Set

import django
from django.contrib.postgres.fields import ArrayField
from django.db.models import CASCADE, TextChoices, SET_NULL, IntegerChoices, Q
from django.urls import reverse
from django_extensions.db.models import TimeStampedModel
from frozendict import frozendict
from classification.enums import AlleleOriginBucket, ShareLevel, SpecialEKeys
from django.db import models, transaction
from classification.models import Classification, ImportedAlleleInfo, EvidenceKeyMap, ClassificationModification, \
    ConditionResolved
from genes.models import GeneSymbol
from library.utils import first
from ontology.models import OntologyTerm
from snpdb.models import Allele, Lab


classification_grouping_search_term_signal = django.dispatch.Signal()  # args: "grouping", expects iterable of ClassificationGroupingSearchTermStub


# TODO this needs to be moved to Classification
class ClassificationQualityLevel(TextChoices):
    STANDARD = "S", "Standard"
    LEGACY = "L", "Legacy"
    INCOMPLETE = "I", "Incomplete"


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
    NOT_COMPARABLE_OVERLAP = 20, "Multiple Submitters"  # e.g. no method to work out discordance
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

    @property
    def allele_origin_bucket_obj(self):
        return AlleleOriginBucket(self.allele_origin_bucket)

    class Meta:
        unique_together = ("allele_grouping", "allele_origin_bucket")

    overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_SHARED_RECORDS)
    dirty = models.BooleanField(default=True)
    # classification_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)
    # somatic_clinical_significance_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)

    def get_absolute_url(self) -> str:
        return reverse('allele_grouping_detail', kwargs={"allele_grouping_id": self.allele_grouping_id})

    def update(self):
        # FIX-ME re-do this once we work out how we handle discordance within a single lab ClassificationGrouping

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


@dataclass
class ClassificationSubGrouping:
    latest_modification: ClassificationModification
    count: int

    @staticmethod
    def from_modifications(modifications: list[ClassificationModification]) -> list['ClassificationSubGrouping']:
        # subgroup by classification value for germline and
        by_status: dict[str, ClassificationSubGrouping] = {}
        for mod in modifications:
            key = str(mod.somatic_clinical_significance_value) + str(mod.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
            if existing := by_status.get(key):
                existing.count += 1
                existing.latest_modification = max(mod, existing.latest_modification, key=lambda m: (m.curated_date_check, m.pk))
            else:
                by_status[key] = ClassificationSubGrouping(latest_modification=mod, count=1)
        # FIXME sort
        return list(by_status.values())


class ClassificationGrouping(TimeStampedModel):
    # key
    allele_origin_grouping = models.ForeignKey(AlleleOriginGrouping, on_delete=models.CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices)
    share_level = models.CharField(max_length=16, choices=ShareLevel.choices())

    @property
    def share_level_obj(self):
        return ShareLevel(self.share_level)

    quality_level = models.CharField(max_length=1, choices=ClassificationQualityLevel.choices, default=ClassificationQualityLevel.STANDARD)
    classification_bucket = models.CharField(max_length=1, choices=ClassificationClassificationBucket.choices, default=ClassificationClassificationBucket.NO_DATA)
    classification_count = models.IntegerField(default=0)

    dirty = models.BooleanField(default=True)

    def dirty_up(self):
        self.dirty = True
        self.save(update_fields=["dirty"])
        self.allele_origin_grouping.dirty = True
        self.allele_origin_grouping.save(update_fields=["dirty"])

    def sub_groupings(self) -> list[ClassificationSubGrouping]:
        return ClassificationSubGrouping.from_modifications(self.classification_modifications)

    conditions = models.JSONField(null=True, blank=True)
    conflicting_ratings = models.BooleanField(default=False)

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
    latest_criteria = ArrayField(models.CharField(max_length=50), null=True, blank=True)
    latest_curation_date = models.DateField(null=True, blank=True)

    def __lt__(self, other):
        def _sort_value(obj: ClassificationGrouping):
            return obj.share_level_obj, obj.lab
        return _sort_value(self) < _sort_value(other)

    def get_absolute_url(self):
        return reverse('classification_grouping_detail', kwargs={"classification_grouping_id": self.pk})

    @staticmethod
    def desired_grouping_for_classification(classification: Classification) -> Optional['ClassificationGrouping']:
        # withdrawn classifications are removed from groupings
        if classification.withdrawn:
            return None

        allele = classification.allele_object
        lab = classification.lab
        share_level = classification.share_level
        if allele:
            allele_grouping, _ = AlleleGrouping.objects.get_or_create(allele=allele)
            allele_origin_grouping, _ = AlleleOriginGrouping.objects.get_or_create(allele_grouping=allele_grouping, allele_origin_bucket=classification.allele_origin_bucket)
            allele_origin_grouping.dirty = True
            allele_origin_grouping.save(update_fields=["dirty"])

            grouping, _ = ClassificationGrouping.objects.get_or_create(
                allele_origin_grouping=allele_origin_grouping,
                lab=lab,
                allele_origin_bucket=classification.allele_origin_bucket,
                share_level=share_level
            )
            return grouping
        return None

    @staticmethod
    def assign_grouping_for_classification(classification: Classification, force_dirty_up=True):
        if desired_grouping := ClassificationGrouping.desired_grouping_for_classification(classification):
            entry, is_new = ClassificationGroupingEntry.objects.get_or_create(
                classification=classification,
                defaults={"grouping": desired_grouping}
            )
            if not is_new:
                old_grouping = entry.grouping
                if entry.grouping != desired_grouping:
                    entry.grouping = desired_grouping
                    entry.save()

                    old_grouping.dirty_up()
                    desired_grouping.dirty_up()
                elif force_dirty_up:
                    entry.dirty_up()
        else:
            # if we don't even have an allele, make sure we are removed from any grouping
            if existing := ClassificationGroupingEntry.objects.filter(classification=classification).first():
                grouping = existing.grouping
                grouping.dirty_up()
                existing.delete()

    @cached_property
    def classification_modifications(self) -> list[ClassificationModification]:
        # show in date order
        all_classifications = self.classificationgroupingentry_set.values_list("classification", flat=True)
        return list(sorted(ClassificationModification.objects.filter(classification_id__in=all_classifications, is_last_published=True), key=lambda mod: mod.curated_date_check))

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
            best_classification = first(all_modifications)

            self.latest_classification_modification = best_classification
            self.latest_curation_date = best_classification.curated_date
            self.latest_allele_info = best_classification.classification.allele_info

            # TODO check for dirty values
            all_terms: Set[OntologyTerm] = set()
            all_free_text_conditions: Set[str] = set()

            # UPDATE CLASSIFICATION / CLINICAL SIGNIFICANCE
            # TODO - CLINICAL SIGNIFICANCE
            all_classification_buckets: set[ClassificationClassificationBucket] = set()
            all_zygosities: set[str] = set()

            for modification in all_modifications:
                if condition := modification.classification.condition_resolution_obj:
                    all_terms |= set(condition.terms)
                    if plain_text := condition.plain_text:
                        all_free_text_conditions.add(plain_text)
                elif condition_text := modification.get(SpecialEKeys.CONDITION):
                    all_free_text_conditions.add(condition_text)

                all_zygosities |= set(modification.get_value_list(SpecialEKeys.ZYGOSITY))

                # only store valid terms as quick links to the classification
                all_terms = {term for term in all_terms if term.is_valid_for_condition and not term.is_stub}
                # self._update_conditions(all_terms)

            evidence_map = EvidenceKeyMap.instance()

            # the below shouldn't happen but has in development environemnts
            if None in all_zygosities:
                all_zygosities.remove(None)

            all_zygosities = list(sorted(evidence_map[SpecialEKeys.ZYGOSITY].sort_values(all_zygosities)))
            self.zygosity_values = all_zygosities

            self.conditions = ConditionResolved(
                terms=list(sorted(all_terms)),
                plain_text=list(sorted(all_free_text_conditions))
            ).to_json(include_join=False)

            bucket_count = len(all_classification_buckets)
            if bucket_count > 1 and ClassificationClassificationBucket.NO_DATA in all_classification_buckets:
                all_classification_buckets.remove(ClassificationClassificationBucket.NO_DATA)
                bucket_count = len(all_classification_buckets)

            use_bucket: ClassificationClassificationBucket
            match bucket_count:
                case 0:
                    use_bucket = ClassificationClassificationBucket.OTHER
                case 1:
                    use_bucket = first(all_classification_buckets)
                case _:
                    use_bucket = ClassificationClassificationBucket.CONFLICTING

            self.classification_bucket = use_bucket
            self.classification_count = len(all_modifications)

            # update search terms

            # TODO do a better syncing of existing values

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
