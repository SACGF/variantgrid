from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date, datetime
from itertools import groupby
from typing import Optional

from django.db import models, transaction
from django.db.models import CASCADE, SET_NULL, F, QuerySet, Window
from django.db.models.fields.json import KeyTextTransform
from django.db.models.functions import Coalesce, RowNumber
from django.utils.timezone import localdate, now

from classification.enums import AlleleOriginBucket, ClinicalSignificance, SpecialEKeys
from classification.models import Classification, ClassificationModification, EvidenceKeyMap
from genes.models import GeneSymbol
from snpdb.models import Lab


class ReclassificationEvent(models.Model):
    """
    A timeline of clinical significance states, one row per state a classification has held. The first
    row for a classification (from_clinical_significance = None) is its initial classification, every
    subsequent row is a reclassification. Records curated on the somatic axis are left out - see
    ReclassificationEventBuilder.tracked_classifications_qs.

    Derived from published modification history by ReclassificationEventBuilder.bring_up_to_date,
    which the analytics page and a nightly task both call.
    @see /classification/reclassification_analytics
    """

    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    """ Denormalised from classification, so the page can filter without a join """

    gene_symbol = models.ForeignKey(GeneSymbol, null=True, blank=True, on_delete=SET_NULL)
    """ Denormalised from allele_info at write time, drives the VUS burden by gene chart """

    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices,
                                            default=AlleleOriginBucket.GERMLINE)

    from_clinical_significance = models.CharField(max_length=1, choices=ClinicalSignificance.CHOICES,
                                                  null=True, blank=True)
    """ None marks the initial classification """

    to_clinical_significance = models.CharField(max_length=1, choices=ClinicalSignificance.CHOICES)

    from_modification = models.ForeignKey(ClassificationModification, null=True, blank=True,
                                          on_delete=CASCADE, related_name='+')
    to_modification = models.ForeignKey(ClassificationModification, on_delete=CASCADE, related_name='+')

    reclassified_date = models.DateField(db_index=True)
    """ Date of to_modification - for external labs this is the sync date, not the curation date """

    step = models.PositiveSmallIntegerField()
    """ 1-based position in this classification's timeline """

    significance_delta = models.SmallIntegerField(null=True, blank=True)
    """ Signed ClinicalSignificance.distance, positive is the benign direction, None for OTHER or initial """

    class Meta:
        indexes = [
            models.Index(fields=["lab", "reclassified_date"]),
            models.Index(fields=["from_clinical_significance", "to_clinical_significance"]),
        ]
        constraints = [
            models.UniqueConstraint(fields=["classification", "step"], name="unique_reclassification_step")
        ]

    def __str__(self):
        return f"{self.classification_id} step {self.step} : {self.from_label} -> {self.to_label}"

    @property
    def from_label(self) -> str:
        return ClinicalSignificance.SHORT_LABELS.get(self.from_clinical_significance)

    @property
    def to_label(self) -> str:
        return ClinicalSignificance.SHORT_LABELS.get(self.to_clinical_significance)

    @property
    def is_reclassification(self) -> bool:
        return self.from_clinical_significance is not None

    @staticmethod
    def reclassifications_qs() -> QuerySet['ReclassificationEvent']:
        """ Excludes initial classifications, i.e. only rows where the significance actually moved """
        return ReclassificationEvent.objects.filter(from_clinical_significance__isnull=False)

    @staticmethod
    def latest_state_qs(as_of: Optional[date] = None,
                        base_qs: Optional[QuerySet['ReclassificationEvent']] = None) -> QuerySet['ReclassificationEvent']:
        """
        The state every classification was in on a given date (now, if as_of is None) - the last row per
        classification that had happened by then. Drives the denominators for the VUS reclassification rate
        and the current-significance counts for the VUS burden by gene.
        """
        qs = base_qs if base_qs is not None else ReclassificationEvent.objects.all()
        if as_of:
            qs = qs.filter(reclassified_date__lte=as_of)
        latest_pks = qs.alias(
            _step_from_latest=Window(
                expression=RowNumber(),
                partition_by=[F('classification_id')],
                order_by=F('step').desc()
            )
        ).filter(_step_from_latest=1).values('pk')
        # pinned by pk so callers can keep filtering without the extra criteria folding into the window
        return ReclassificationEvent.objects.filter(pk__in=latest_pks)


class ReclassificationEventBuildState(models.Model):
    """ How far through classification history the timelines have been built. @see ReclassificationEvent """

    built_to = models.DateTimeField(null=True, blank=True)
    """ Timelines are current for every classification modified at or before this """

    last_run = models.DateTimeField(null=True, blank=True)

    def __str__(self):
        return f"Reclassification timelines built to {self.built_to}"

    @staticmethod
    def instance() -> 'ReclassificationEventBuildState':
        state, _ = ReclassificationEventBuildState.objects.get_or_create(pk=1)
        return state


@dataclass(frozen=True)
class ReclassificationBuildResult:
    """ What one call to ReclassificationEventBuilder.bring_up_to_date got through """
    classifications: int
    """ Timelines rebuilt this run """
    events: int
    """ Event rows written """
    outstanding: int
    """ Classifications still waiting, 0 once caught up """
    built_to: Optional[datetime]
    last_run: Optional[datetime]


@dataclass(frozen=True)
class ReclassificationStep:
    """ One state in a classification's significance timeline, before it's saved """
    from_clinical_significance: Optional[str]
    to_clinical_significance: str
    from_modification_id: Optional[int]
    to_modification_id: int
    reclassified_date: date
    step: int

    @property
    def significance_delta(self) -> Optional[int]:
        if self.from_clinical_significance:
            return ClinicalSignificance.distance(self.from_clinical_significance, self.to_clinical_significance)
        return None


class ReclassificationTimeline:
    """
    Walks a classification's published modifications and works out where its clinical significance moved.
    The significance column is set during patch rather than publish, so where it hasn't caught up we read
    the value straight out of published_evidence.
    """

    @staticmethod
    def published_modifications_qs(classification_qs: QuerySet[Classification]) -> QuerySet:
        """
        Just enough of each published modification to build a timeline - the evidence significance is pulled
        out in SQL so we never drag the whole published_evidence blob back.
        """
        evidence_clinical_significance = Coalesce(
            KeyTextTransform('value', KeyTextTransform(SpecialEKeys.CLINICAL_SIGNIFICANCE, 'published_evidence')),
            KeyTextTransform(SpecialEKeys.CLINICAL_SIGNIFICANCE, 'published_evidence')
        )
        return ClassificationModification.objects \
            .filter(published=True, classification__in=classification_qs) \
            .annotate(evidence_clinical_significance=evidence_clinical_significance) \
            .values('id', 'classification_id', 'created', 'clinical_significance', 'evidence_clinical_significance') \
            .order_by('classification_id', 'created')

    @staticmethod
    def significance_of(clinical_significance: Optional[str], evidence_clinical_significance: Optional[str]) -> Optional[str]:
        if clinical_significance:
            return clinical_significance
        if evidence_clinical_significance:
            options = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE) \
                .matched_options(evidence_clinical_significance)
            return options[0].get('vg')
        return None

    @staticmethod
    def steps_for(modification_rows: Iterable[dict]) -> list[ReclassificationStep]:
        """
        :param modification_rows: published modification rows for a single classification, oldest first
        """
        steps: list[ReclassificationStep] = []
        previous_significance: Optional[str] = None
        previous_modification_id: Optional[int] = None

        for row in modification_rows:
            significance = ReclassificationTimeline.significance_of(
                row['clinical_significance'], row['evidence_clinical_significance'])
            if significance and significance != previous_significance:
                steps.append(ReclassificationStep(
                    from_clinical_significance=previous_significance,
                    to_clinical_significance=significance,
                    from_modification_id=previous_modification_id if previous_significance else None,
                    to_modification_id=row['id'],
                    reclassified_date=localdate(row['created']),
                    step=len(steps) + 1
                ))
                previous_significance = significance
            previous_modification_id = row['id']

        return steps


class ReclassificationEventBuilder:
    """ Turns timelines into ReclassificationEvent rows """

    BATCH_SIZE = 2000

    @staticmethod
    def tracked_classifications_qs() -> QuerySet[Classification]:
        """
        Everything curated on the germline clinical_significance axis. Somatic records are curated on
        somatic:clinical_significance instead, so they get their own timeline once that data matures -
        allele_origin_bucket is stored on each event so the page can tell germline from unknown.
        """
        return Classification.objects.exclude(allele_origin_bucket=AlleleOriginBucket.SOMATIC)

    @staticmethod
    def denormalised_columns(classification_ids: list[int]) -> dict[int, dict]:
        """ Only the tracked classifications come back, so callers can hand in a plain "touched" set """
        columns = {}
        for pk, lab_id, allele_origin_bucket, gene_symbol_38, gene_symbol_37 in \
                ReclassificationEventBuilder.tracked_classifications_qs() \
                .filter(pk__in=classification_ids) \
                .values_list('pk', 'lab_id', 'allele_origin_bucket',
                             'allele_info__grch38__gene_symbol_id', 'allele_info__grch37__gene_symbol_id'):
            columns[pk] = {
                "lab_id": lab_id,
                "allele_origin_bucket": allele_origin_bucket,
                "gene_symbol_id": gene_symbol_38 or gene_symbol_37
            }
        return columns

    @staticmethod
    def events_for(classification_id: int, steps: list[ReclassificationStep], columns: dict) -> list[ReclassificationEvent]:
        return [
            ReclassificationEvent(
                classification_id=classification_id,
                from_clinical_significance=step.from_clinical_significance,
                to_clinical_significance=step.to_clinical_significance,
                from_modification_id=step.from_modification_id,
                to_modification_id=step.to_modification_id,
                reclassified_date=step.reclassified_date,
                step=step.step,
                significance_delta=step.significance_delta,
                **columns
            ) for step in steps
        ]

    @staticmethod
    def rebuild(classification_qs: QuerySet[Classification]) -> int:
        """
        Deletes the timeline of every classification in the queryset and rebuilds it for those still on the
        germline axis, so a record that turns somatic loses its timeline. Modifications stream out of a
        single query, and each batch is deleted only as it is rewritten, so the queryset can be one that the
        rewrite itself would otherwise change the answer to.
        :return: the number of event rows written
        """
        rows_written = 0
        classification_ids: list[int] = []
        steps_by_classification: dict[int, list[ReclassificationStep]] = {}

        def flush() -> int:
            if not classification_ids:
                return 0
            ReclassificationEvent.objects.filter(classification_id__in=classification_ids).delete()
            columns_by_classification = ReclassificationEventBuilder.denormalised_columns(classification_ids)
            events = []
            for classification_id, steps in steps_by_classification.items():
                if columns := columns_by_classification.get(classification_id):
                    events += ReclassificationEventBuilder.events_for(classification_id, steps, columns)
            ReclassificationEvent.objects.bulk_create(events, batch_size=ReclassificationEventBuilder.BATCH_SIZE)
            classification_ids.clear()
            steps_by_classification.clear()
            return len(events)

        modification_rows = ReclassificationTimeline.published_modifications_qs(classification_qs)
        for classification_id, rows in groupby(modification_rows.iterator(chunk_size=10000),
                                               key=lambda row: row['classification_id']):
            classification_ids.append(classification_id)
            if steps := ReclassificationTimeline.steps_for(rows):
                steps_by_classification[classification_id] = steps
            if len(classification_ids) >= ReclassificationEventBuilder.BATCH_SIZE:
                rows_written += flush()

        rows_written += flush()
        return rows_written

    @staticmethod
    def bring_up_to_date(max_classifications: Optional[int] = None) -> ReclassificationBuildResult:
        """
        Rebuilds the timeline of every classification touched since the last run, then advances the
        watermark. Safe to call from any view that reads ReclassificationEvent, and from the nightly task.
        :param max_classifications: hand the batch to the nightly run rather than build more than this now
        """
        started_at = now()
        ReclassificationEventBuildState.instance()

        with transaction.atomic():
            state = ReclassificationEventBuildState.objects \
                .select_for_update(skip_locked=True).filter(pk=1).first()
            if not state:
                # another build holds the row, so report the watermark it started from
                return ReclassificationEventBuilder._result_for(
                    ReclassificationEventBuildState.instance(), outstanding=0)

            touched_qs = Classification.objects.all()
            if state.built_to:
                touched_qs = touched_qs.filter(modified__gt=state.built_to)

            classification_count = touched_qs.count()
            if max_classifications is not None and classification_count > max_classifications:
                return ReclassificationEventBuilder._result_for(state, outstanding=classification_count)

            events = ReclassificationEventBuilder.rebuild(touched_qs)
            state.built_to = started_at
            state.last_run = now()
            state.save()
            return ReclassificationEventBuilder._result_for(
                state, outstanding=0, classifications=classification_count, events=events)

    @staticmethod
    def _result_for(state: ReclassificationEventBuildState, outstanding: int,
                    classifications: int = 0, events: int = 0) -> ReclassificationBuildResult:
        return ReclassificationBuildResult(classifications=classifications, events=events,
                                           outstanding=outstanding, built_to=state.built_to,
                                           last_run=state.last_run)
