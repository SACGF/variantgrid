from functools import cached_property
from typing import Optional, Set

from django.contrib.postgres.fields import ArrayField
from django.db.models import CASCADE, TextChoices, SET_NULL, IntegerChoices
from django.urls import reverse
from django_extensions.db.models import TimeStampedModel

from classification.criteria_strengths import CriteriaStrength
from classification.enums import AlleleOriginBucket, ShareLevel, SpecialEKeys, CriteriaEvaluation
from django.db import models, transaction

from classification.models import Classification, ImportedAlleleInfo, EvidenceKeyMap, ClassificationModification, \
    ConditionResolved, SomaticClinicalSignificanceValue
from genes.models import GeneSymbol
from library.utils import first
from ontology.models import OntologyTerm
from snpdb.models import Allele, Lab


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
    BENIGN_DISCORDANCE = 50, "Benign Discordance"
    PATHOGENIC_DISCORDANCE = 60, "Pathogenic Discordance"


class AlleleGrouping(TimeStampedModel):
    allele = models.OneToOneField(Allele, on_delete=models.CASCADE)
    # TODO probably some more summary fields ew could have here?
    # otherwise what is this serving that Allele isn't?

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
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices)
    @property
    def allele_origin_bucket_obj(self):
        return AlleleOriginBucket(self.allele_origin_bucket)

    class Meta:
        unique_together = ("allele_grouping", "allele_origin_bucket")

    overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_SHARED_RECORDS)
    dirty = models.BooleanField(default=True)
    classification_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)

    def update(self):
        # not sure if this is the best solution, or if an AlleleOriginGrouping should just refuse to update if attached allele groupings are dirty
        all_classification_values = set()
        classification_groupings: list[ClassificationGrouping] = list(self.classificationgrouping_set.all())
        for cg in classification_groupings:
            if cg.dirty:
                cg.update()
        shared_groupings = [cg for cg in classification_groupings if cg.share_level_obj.is_discordant_level]
        for cg in shared_groupings:
            all_classification_values |= set(cg.classification_values)

        self.classification_values = list(all_classification_values)

        overlap_status: OverlapStatus
        if len(shared_groupings) == 0:
            overlap_status = OverlapStatus.NO_SHARED_RECORDS
        elif len(shared_groupings) == 1:
            overlap_status = OverlapStatus.SINGLE_SUBMITTER
        elif self.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
            overlap_status = OverlapStatus.NOT_COMPARABLE_OVERLAP
        else:
            bucket_mapping = EvidenceKeyMap.instance().get(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("bucket")
            buckets = {bucket_mapping.get(class_value) for class_value in self.classification_values}
            if None in buckets:
                buckets.remove(None)

            if len(buckets) > 1:
                # discordant
                if "P" in all_classification_values or "LP" in all_classification_values:
                    overlap_status = OverlapStatus.PATHOGENIC_DISCORDANCE
                else:
                    overlap_status = OverlapStatus.BENIGN_DISCORDANCE
            else:
                if len(all_classification_values) > 1:
                    overlap_status = OverlapStatus.CONFIDENCE
                else:
                    # complete agreement
                    overlap_status = OverlapStatus.AGREEMENT

        self.overlap_status = overlap_status
        self.dirty = False
        self.save()

        # TODO manage pending classifications


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

    # All values before need to be nullable as they wont be populated in the time between creating the ClassificationGrouping and updating it
    classification_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)  # values from EvidenceKey.CLINICAL_SIGNIFICANCE
    somatic_clinical_significance_values = ArrayField(models.CharField(max_length=30), null=True, blank=True)
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
    latest_classification = models.ForeignKey(ClassificationModification, on_delete=SET_NULL, null=True, blank=True)
    latest_criteria = ArrayField(models.CharField(max_length=50), null=True, blank=True)
    latest_allele_info = models.ForeignKey(ImportedAlleleInfo, on_delete=SET_NULL, null=True, blank=True)
    latest_curation_date = models.DateField(null=True, blank=True)

    somatic_clinical_significance_sort = models.IntegerField(db_index=True, null=True, blank=True)
    classification_sort_value = models.TextField(null=True, blank=True)

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
    def assign_grouping_for_classification(classification: Classification):
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

                    old_grouping.dirty = True
                    old_grouping.save(update_fields=["dirty"])
                    desired_grouping.dirty = True
                    desired_grouping.save(update_fields=["dirty"])
        else:
            # if we don't even have an allele, make sure we are removed from any grouping
            if existing := ClassificationGroupingEntry.objects.filter(classification=classification).first():
                grouping = existing.grouping
                existing.delete()
                grouping.dirty = True
                grouping.save(update_fields=["dirty"])

    def _update_gene_symbols(self):
        all_gene_symbols = set()
        all_imported_allele_ids = ImportedAlleleInfo.objects.filter(
            pk__in=self.classificationgroupingentry_set.all().values_list('classification__allele_info_id',
                                                                          flat=True)
        )
        for imported_allele_info in all_imported_allele_ids:
            all_gene_symbols |= ClassificationGrouping.gene_symbols_for_imported_allele_info(imported_allele_info)

        # delete any gene symbols links that are no longer relevant
        ClassificationGroupingGeneSymbol.objects.filter(grouping=self).exclude(
            gene_symbol__in=all_gene_symbols).delete()
        # bulk create all the ones it should be in
        desired_entries = [ClassificationGroupingGeneSymbol(grouping=self, gene_symbol=gene_symbol) for gene_symbol in
                           all_gene_symbols]
        output = ClassificationGroupingGeneSymbol.objects.bulk_create(desired_entries, ignore_conflicts=True)

    @cached_property
    def gene_symbols(self) -> list[GeneSymbol]:
        return list(sorted(gss.gene_symbol for gss in self.classificationgroupinggenesymbol_set.select_related("gene_symbol")))

    def _update_conditions(self, terms: list[OntologyTerm]):
        ClassificationGroupingCondition.objects.filter(grouping=self).exclude(
            ontology_term_id__in=[term.pk for term in terms]).delete()
        desired_entries = [ClassificationGroupingCondition(grouping=self, ontology_term_id=term.pk) for term in
                           terms]
        ClassificationGroupingCondition.objects.bulk_create(desired_entries, ignore_conflicts=True)

    @transaction.atomic
    def update(self):

        self.allele_origin_grouping.dirty = True
        self.allele_origin_grouping.save(update_fields=["dirty"])

        # UPDATE GENE SYMBOLS
        self._update_gene_symbols()

        all_modifications: list[ClassificationModification] = list(ClassificationModification.objects.filter(
            is_last_published=True,
            classification__in=self.classificationgroupingentry_set.values_list("classification"),
            classification__withdrawn=False
        ).select_related('classification'))
        all_modifications: list[ClassificationModification] = list(sorted(all_modifications, key=lambda x: x.curated_date, reverse=True))

        if all_modifications:
            # FIND THE MOST RECENT CLASSIFICATION
            best_classification = first(all_modifications)

            self.latest_classification = best_classification
            self.latest_curation_date = best_classification.curated_date
            self.latest_allele_info = best_classification.classification.allele_info
            # TODO, calculate the somatic sort order here so we can remove it from classification
            self.somatic_clinical_significance_sort = best_classification.somatic_clinical_significance_sort

            best_clin_sig = best_classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            self.classification_sort_value = classification_sort_order(best_clin_sig)

            # GET ALL CRITERIA
            def criteria_converter(cm: ClassificationModification) -> set[CriteriaStrength]:
                strengths: set[CriteriaStrength] = set()
                for e_key in EvidenceKeyMap.cached().criteria():
                    strength = cm.get(e_key.key)
                    if CriteriaEvaluation.is_met(strength):
                        strengths.add(CriteriaStrength(e_key, strength))
                for amp_level, letter in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items():
                    if value := cm.get_value_list(amp_level):
                        e_key = EvidenceKeyMap.cached_key(amp_level)
                        for sub_value in value:
                            sub_value_label = e_key.pretty_value(sub_value)
                            strengths.add(CriteriaStrength(
                                ekey=EvidenceKeyMap.cached_key(amp_level),
                                custom_strength=f"{letter}_{sub_value_label}")
                            )
                return list(sorted(strengths))

            self.latest_criteria = [str(crit) for crit in criteria_converter(self.latest_classification)]

            if self.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
                all_somatic_clin_sigs: set[SomaticClinicalSignificanceValue] = set()
                for cm in all_modifications:
                    if somatic_clin_sig := cm.somatic_clinical_significance_value:
                        all_somatic_clin_sigs.add(somatic_clin_sig)
                sorted_somatic_clin_sig = sorted(all_somatic_clin_sigs)
                self.somatic_clinical_significance_values = [scs.as_str for scs in sorted_somatic_clin_sig]
            else:
                self.somatic_clinical_significance_values = None

            all_classification_values: Set[str] = set()
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

                classification_value = modification.classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                all_classification_values.add(classification_value)
                classification_bucket = ClassificationClassificationBucket.bucket_for_classification(classification_value)
                all_classification_buckets.add(classification_bucket)

                all_zygosities |= set(modification.get_value_list(SpecialEKeys.ZYGOSITY))

                # only store valid terms as quick links to the classification
                all_terms = {term for term in all_terms if term.is_valid_for_condition and not term.is_stub}
                self._update_conditions(all_terms)

            evidence_map = EvidenceKeyMap.instance()
            all_classification_values = evidence_map[SpecialEKeys.CLINICAL_SIGNIFICANCE].sort_values(all_classification_values)
            self.classification_values = all_classification_values

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

            # FIXME need the concept of an OTHER that's low priority and doesn't interfere with anything
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
            self.dirty = False
            self.save()

        else:
            # there are no classifications, time to die
            self.delete()

    @staticmethod
    def gene_symbols_for_imported_allele_info(imported_allele_info: ImportedAlleleInfo) -> Set[GeneSymbol]:
        all_gene_symbols: Set[GeneSymbol] = set()
        for rb in imported_allele_info.resolved_builds:
            if gene_symbol := rb.gene_symbol:
                all_gene_symbols.add(gene_symbol)

            # For now exclude transcript versions linking to other gene symbols, revisit
            # if transcript_version := rb.transcript_version:
            #     if gene_version := transcript_version.gene_version:
            #         if gene_symbol := gene_version.gene_symbol:
            #             all_gene_symbols.add(gene_symbol)

        if c_hgvs := imported_allele_info.imported_c_hgvs_obj:
            # TODO if it's alias look up aliases
            if included_gene_symbol := GeneSymbol.objects.filter(symbol=c_hgvs.gene_symbol).first():
                all_gene_symbols.add(included_gene_symbol)

            # # need to invert how GeneSymbolAliasesMeta works by starting from symbol and going back to
            # if annotation := VariantAnnotation.objects.filter(
            #     version=VariantAnnotationVersion.latest(genome_build=rb.genome_build),
            #     variant=rb.variant
            # ).first():
            #     annotation: VariantAnnotation
            #     gene = annotation.gene
            #     GeneSymbolAliasesMeta.objects.filter()
        return all_gene_symbols

    @cached_property
    def allele(self) -> Allele:
        return self.allele_origin_grouping.allele

    @cached_property
    def classification_modifications(self) -> list[ClassificationModification]:
        # show in date order
        all_classifications = self.classificationgroupingentry_set.values_list("classification", flat=True)
        return list(sorted(ClassificationModification.objects.filter(classification_id__in=all_classifications, is_last_published=True), key=lambda mod: mod.curated_date_check))

    def to_json(self):
        scs = self.latest_classification.somatic_clinical_significance_value
        return {
            "lab": str(self.lab),
            "classification": self.latest_classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
            "somatic_clinical_significance": scs.as_json() if scs else None,
            "curated_date": str(self.latest_curation_date)  # fixme, need a quick way for yyyy-mm-dd
        }

    @staticmethod
    def update_all_dirty():
        for dirty in ClassificationGrouping.objects.filter(dirty=True).iterator():
            dirty.update()
        for dirty in AlleleOriginGrouping.objects.filter(dirty=True).iterator():
            dirty.update()


class ClassificationGroupingGeneSymbol(TimeStampedModel):
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE, db_index=False)

    def __str__(self):
        return f"Grouping {self.grouping} : GeneSymbol {self.gene_symbol}"

    class Meta:
        unique_together = ("grouping", "gene_symbol")


class ClassificationGroupingCondition(TimeStampedModel):
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)

    def __str__(self):
        return f"Grouping {self.grouping} : OntologyTerm {self.ontology_term}"

    class Meta:
        unique_together = ("grouping", "ontology_term")


class ClassificationGroupingEntry(TimeStampedModel):
    # this is just here so this model can stay completely separate from Classification
    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)

    def dirty_up(self):
        self.grouping.dirty = True
        self.grouping.save(update_fields=["dirty"])
        self.grouping.allele_origin_grouping.dirty = True
        self.grouping.allele_origin_grouping.save(update_fields=["dirty"])

    def __str__(self):
        return f"Grouping {self.grouping} : Classification {self.classification}"

    class Meta:
        unique_together = ("grouping", "classification")
