"""
The Classify & Report tab on the sample and patient pages - tags that are asking to be classified, the
classifications made for the case, and what's needed to turn one into the other.

A case is a set of samples (one for a sample page, the patient's samples for a patient page). A tagging is in
the case's queue when its tag is in the classify queue vocabulary (Tag.requires_classification) and it belongs
to one of the case's samples - either by its own sample FK, or by having been made in an analysis whose proband
is one of them.

A tagging is done once it is resolved against a classification (@see analysis.variant_tag_operations), which
happens by itself when the classification is of the tagging's own sample and by the row's button otherwise.
The queue also shows a classification the case already has for the allele, so a variant classified outside this
tab arrives with the link and the button rather than looking untouched.

@see https://github.com/SACGF/variantgrid_sapath/issues/246
"""
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional

from django.contrib.auth.models import User
from django.db.models import Q, QuerySet
from django.urls import reverse

from analysis.models import Analysis, VariantTag
from analysis.variant_tag_operations import sample_carries_variant
from classification.enums import AlleleOriginBucket, SpecialEKeys
from classification.models import Classification, ClassificationModification
from patients.models import Patient
from snpdb.models import Lab, Sample, Tag


@dataclass(frozen=True)
class PreviousClassification:
    """ An existing classification of the same allele - the scientist decides whether it applies here """
    modification: ClassificationModification
    own_lab: bool

    @property
    def can_copy(self) -> bool:
        """ An external lab's record was curated under another config, assertion method and namespaces, so it
            is context for the decision rather than something to copy evidence out of """
        return not self.modification.classification.lab.external

    @property
    def condition(self) -> str:
        return self.modification.condition_text or ""

    @property
    def interpretation_summary(self) -> str:
        return self.modification.get(SpecialEKeys.INTERPRETATION_SUMMARY) or ""


@dataclass
class ClassifyQueueRow:
    """ One tagged variant in a case, and how far it has got """
    variant_tag: VariantTag
    sample: Optional[Sample] = None
    classification: Optional[Classification] = None
    previous: list[PreviousClassification] = field(default_factory=list)

    @property
    def done(self) -> bool:
        return self.variant_tag.is_resolved

    @property
    def can_resolve(self) -> bool:
        """ The row has something to point the to-do at, but nobody has said it is the right person yet """
        return self.classification is not None and not self.done

    @property
    def needs_classification(self) -> bool:
        """ Nothing has been classified for the allele yet - a row that has one is waiting on "Clear tag",
            not on more curation, so the wizard walks past it """
        return not self.done and self.classification is None

    @property
    def full_form_url(self) -> str:
        """ The full create page for this tagging - the analysis one where the tag was made in an analysis, so
            creating there links the record to it and clears the tag """
        if analysis_id := self.variant_tag.analysis_id:
            return reverse("create_classification_for_variant_tag",
                           kwargs={"analysis_id": analysis_id, "variant_tag_id": self.variant_tag.pk})
        return reverse("create_classification_for_variant",
                       kwargs={"variant_id": self.variant_tag.variant_id,
                               "genome_build_name": self.variant_tag.genome_build.name})

    @property
    def copyable(self) -> list[PreviousClassification]:
        """ The previous classifications this row can be started from - external labs are shown, not copied """
        return [previous for previous in self.previous if previous.can_copy]

    @property
    def conditions(self) -> list[str]:
        """ Distinct conditions the allele has been classified for before - they rarely agree, which is why
            picking one is left to the scientist """
        conditions = []
        for previous in self.previous:
            if (condition := previous.condition) and condition not in conditions:
                conditions.append(condition)
        return conditions

    @cached_property
    def gene_symbol(self):
        return self.variant_tag.gene_symbol

    @property
    def allele_origin_bucket(self) -> Optional[AlleleOriginBucket]:
        """ The bucket a classification made from this tagging goes into - the tag says which work it is
            part of, and nothing germline is ever copied into a somatic record or the other way round """
        return tag_allele_origin_bucket(self.variant_tag.tag)


def tag_allele_origin_bucket(tag: Tag) -> Optional[AlleleOriginBucket]:
    """ Which bucket a tag's classifications belong in - "Both" leaves the choice to the record being copied """
    bucket = AlleleOriginBucket(tag.allele_origin_bucket)
    return bucket if bucket != AlleleOriginBucket.UNKNOWN else None


def _allele_keys(allele_id, variant_id) -> set[tuple[str, int]]:
    """ Tags and classifications match on allele where both have one, otherwise on variant - a tag's allele is
        assigned asynchronously, and a classification only has one once it's matched """
    keys = set()
    if allele_id:
        keys.add(("allele", allele_id))
    if variant_id:
        keys.add(("variant", variant_id))
    return keys


class ClassifyReportCase:
    """ The samples a Classify & Report tab covers, seen as one case """

    def __init__(self, user: User, obj, samples: list[Sample]):
        self.user = user
        self.obj = obj  # Sample or Patient the tab is on
        self.samples = samples

    @staticmethod
    def for_sample(user: User, sample: Sample) -> 'ClassifyReportCase':
        return ClassifyReportCase(user, sample, [sample])

    @staticmethod
    def for_patient(user: User, patient: Patient) -> 'ClassifyReportCase':
        samples = list(Sample.filter_for_user(user).filter(pk__in=patient.get_samples()))
        return ClassifyReportCase(user, patient, samples)

    @property
    def sample_ids(self) -> set[int]:
        return {s.pk for s in self.samples}

    @property
    def labs(self) -> list[Lab]:
        return list(Lab.valid_labs_qs(self.user))

    def _visible_variant_tags(self):
        """ Live classify-queue taggings the user can see - a tagging made in an analysis is the analysis's
            to show (@see VariantTag.can_view) """
        return VariantTag.objects.filter(tag__in=Tag.classify_queue_qs()) \
            .filter(Q(analysis__isnull=True) | Q(analysis__in=Analysis.filter_for_user(self.user))) \
            .select_related("tag", "variant", "allele", "analysis", "sample", "genome_build")

    def _analysis_samples(self, analysis_ids) -> dict[int, list[Sample]]:
        """ Which of the case's samples each analysis contains """
        samples_by_analysis = {}
        by_pk = {s.pk: s for s in self.samples}
        for analysis in Analysis.objects.filter(pk__in=analysis_ids):
            if ours := [by_pk[pk] for pk in sorted(by_pk.keys() & {s.pk for s in analysis.get_samples()})]:
                samples_by_analysis[analysis.pk] = ours
        return samples_by_analysis

    def _shown_to_carriers(self, variant_tag: VariantTag, candidates: list[Sample]) -> bool:
        """ Whether a tagging whose proband we don't know is worth showing this case. Carrying the variant
            doesn't make it a sample's to-do - a relative can be HET for the proband's variant - it just keeps
            the row off the pages it cannot be about. The sample stays unset, so the scientist picks """
        return any(sample_carries_variant(s, variant_tag) for s in candidates)

    def queue_row(self, variant_tag: VariantTag) -> ClassifyQueueRow:
        """ One tagging's row - what the classify dialog is built from """
        sample = variant_tag.sample if variant_tag.sample_id in self.sample_ids else None
        classification = None
        classifications_by_key = self._classifications_by_allele_key()
        for key in _allele_keys(variant_tag.allele_id, variant_tag.variant_id):
            if classification := classifications_by_key.get(key):
                break
        return ClassifyQueueRow(variant_tag=variant_tag, sample=sample, classification=classification,
                                previous=self._previous_by_tag([variant_tag]).get(variant_tag.pk, []))

    def variant_tags(self) -> list[tuple[VariantTag, Optional[Sample]]]:
        """ The case's taggings, each with the sample it is about where the tagging knows """
        tags_qs = self._visible_variant_tags()
        rows = [(vt, vt.sample) for vt in tags_qs.filter(sample__in=self.samples)]

        no_sample = list(tags_qs.filter(sample__isnull=True, analysis__isnull=False))
        samples_by_analysis = self._analysis_samples({vt.analysis_id for vt in no_sample})
        for variant_tag in no_sample:
            candidates = samples_by_analysis.get(variant_tag.analysis_id, [])
            if candidates and self._shown_to_carriers(variant_tag, candidates):
                rows.append((variant_tag, None))

        return rows

    def classification_modifications(self):
        """ Latest published classification of each of the case's samples """
        qs = ClassificationModification.latest_for_user(self.user, published=True,
                                                        classification__sample__in=self.samples)
        return qs.order_by("classification__pk")

    def _case_classifications(self) -> QuerySet[Classification]:
        """ What the case has classified - the published records the user can see plus their own lab's work in
            progress, because a record created from this queue is unpublished until they fill it in """
        published = self.classification_modifications().values_list("classification_id", flat=True)
        own = Classification.filter_for_user(self.user).values("pk")
        return Classification.objects.filter(sample__in=self.samples, withdrawn=False) \
            .filter(Q(pk__in=published) | Q(pk__in=own))

    def _classifications_by_allele_key(self) -> dict[tuple[str, int], Classification]:
        by_key = {}
        for classification in self._case_classifications():
            for key in _allele_keys(classification.allele_id, classification.variant_id):
                by_key[key] = classification
        return by_key

    def _previous_by_tag(self, variant_tags) -> dict[int, list[PreviousClassification]]:
        """ The classifications of each tagged allele the user can see, this user's own labs first """
        alleles = {vt.allele_id for vt in variant_tags if vt.allele_id}
        variants = {vt.variant_id for vt in variant_tags if not vt.allele_id}
        if not (alleles or variants):
            return {}

        lab_ids = {lab.pk for lab in self.labs}
        qs = ClassificationModification.latest_for_user(self.user, published=True) \
            .filter(Q(classification__allele__in=alleles) | Q(classification__variant__in=variants)) \
            .select_related("classification", "classification__lab", "classification__user")

        by_key = defaultdict(list)
        for cm in qs:
            classification = cm.classification
            previous = PreviousClassification(modification=cm,
                                              own_lab=classification.lab_id in lab_ids)
            for key in _allele_keys(classification.allele_id, classification.variant_id):
                by_key[key].append(previous)

        previous_by_tag = {}
        for variant_tag in variant_tags:
            previous = []
            bucket = tag_allele_origin_bucket(variant_tag.tag)
            for key in _allele_keys(variant_tag.allele_id, variant_tag.variant_id):
                for candidate in by_key.get(key, []):
                    if bucket and candidate.modification.classification.allele_origin_bucket != bucket:
                        continue
                    if candidate not in previous:
                        previous.append(candidate)
            # Own lab first, then most recently curated - it's the lab's own work that usually applies
            previous.sort(key=lambda p: (p.own_lab, p.modification.curated_date_check), reverse=True)
            previous_by_tag[variant_tag.pk] = previous
        return previous_by_tag

    def queue_rows(self) -> list[ClassifyQueueRow]:
        variant_tags_and_samples = self.variant_tags()
        variant_tags = [vt for vt, _ in variant_tags_and_samples]
        classifications_by_key = self._classifications_by_allele_key()
        previous_by_tag = self._previous_by_tag(variant_tags)

        rows = []
        for variant_tag, sample in variant_tags_and_samples:
            classification = None
            for key in _allele_keys(variant_tag.allele_id, variant_tag.variant_id):
                if classification := classifications_by_key.get(key):
                    break
            rows.append(ClassifyQueueRow(variant_tag=variant_tag, sample=sample,
                                         classification=classification,
                                         previous=previous_by_tag.get(variant_tag.pk, [])))
        rows.sort(key=lambda row: (str(row.gene_symbol or ""), row.variant_tag.pk))
        return rows


def outstanding_tag_count(rows: list[ClassifyQueueRow]) -> int:
    return sum(1 for row in rows if not row.done)


def tag_summary(rows: list[ClassifyQueueRow]) -> list[tuple[Tag, int]]:
    """ Outstanding work broken down by tag, for the funnel at the top of the tab """
    counts = defaultdict(int)
    for row in rows:
        if not row.done:
            counts[row.variant_tag.tag] += 1
    return sorted(counts.items(), key=lambda tag_count: tag_count[0].pk)
