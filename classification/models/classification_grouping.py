from typing import Optional, Set

from django.db.models import CASCADE, TextChoices
from django_extensions.db.models import TimeStampedModel

from annotation.models import VariantAnnotation, VariantAnnotationVersion
from classification.enums import AlleleOriginBucket, ShareLevel, SpecialEKeys
from django.db import models, transaction

from classification.models import Classification, ImportedAlleleInfo, EvidenceKeyMap
from genes.models import GeneSymbol, GeneSymbolAliasesMeta
from genes.signals.gene_symbol_search import gene_symbol_alias_search
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


class ClassificationGrouping(TimeStampedModel):
    # key
    allele = models.ForeignKey(Allele, on_delete=models.CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices)
    share_level = models.CharField(max_length=16, choices=ShareLevel.choices())
    quality_level = models.CharField(max_length=1, choices=ClassificationQualityLevel.choices, default=ClassificationQualityLevel.STANDARD)

    classification_bucket = models.CharField(max_length=1, choices=ClassificationClassificationBucket.choices, default=ClassificationClassificationBucket.NO_DATA)

    dirty = models.BooleanField(default=False)

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

    @staticmethod
    def desired_grouping_for_classification(classification: Classification) -> Optional['ClassificationGrouping']:
        allele = classification.allele_object
        lab = classification.lab
        share_level = classification.share_level
        if allele:
            grouping, _ = ClassificationGrouping.objects.get_or_create(
                allele=allele,
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
            pk__in=self.classificationgroupingentry_set.all().values_list('classification__imported_allele_id',
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
        ClassificationGroupingGeneSymbol.objects.bulk_create(desired_entries, ignore_conflicts=True)

    @transaction.atomic
    def update_based_on_entries(self):
        # TODO update a lot of fields

        # UPDATE GENE SYMBOLS
        # self._update_gene_symbols

        # UPDATE CLASSIFICATION
        all_classification_buckets: set[ClassificationClassificationBucket] = set()
        for classification_entry in self.classificationgroupingentry_set.select_related('classification').all():
            classification_value = classification_entry.classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            classification_bucket = ClassificationClassificationBucket.bucket_for_classification(classification_value)
            all_classification_buckets.add(classification_bucket)

        use_bucket = ClassificationClassificationBucket.OTHER
        bucket_count = len(all_classification_buckets)

        if bucket_count > 1 and ClassificationClassificationBucket.NO_DATA in bucket_count:
            all_classification_buckets.pop(ClassificationClassificationBucket.NO_DATA)
            bucket_count = len(all_classification_buckets)

        # FIXME need the concept of an OTHER that's low priority and doesn't interfere with anything
        match bucket_count:
            case 0:
                use_bucket = ClassificationClassificationBucket.OTHER
            case 1:
                use_bucket = first(all_classification_buckets)
            case _:
                use_bucket = ClassificationClassificationBucket.CONFLICTING
        self.classification_bucket = use_bucket

        self.dirty = False
        self.save()
        pass

    @staticmethod
    def gene_symbols_for_imported_allele_info(imported_allele_info: ImportedAlleleInfo) -> Set[GeneSymbol]:
        all_gene_symbols: Set[GeneSymbol] = set()
        for rb in imported_allele_info.resolved_builds:
            all_gene_symbols.add(rb.gene_symbol)
            all_gene_symbols.add(rb.transcript_version.gene_version.gene_symbol)

        if c_hgvs := imported_allele_info.imported_c_hgvs_obj:
            # TODO if it's alias look up aliases
            if included_gene_symbol := GeneSymbol.objects.filter(symbol=c_hgvs.gene_symbol):
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

    @staticmethod
    def update_all_dirty():
        for dirty in ClassificationGrouping.objects.filter(dirty=True):
            dirty.update_based_on_entries()


class ClassificationGroupingGeneSymbol(TimeStampedModel):
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)
    # gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)


class ClassificationGroupingConditions(TimeStampedModel):
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)
    condition = models.ForeignKey(OntologyTerm, on_delete=CASCADE)


class ClassificationGroupingEntry(TimeStampedModel):
    # this is just here so this model can stay completely separate from Classification
    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    grouping = models.ForeignKey(ClassificationGrouping, on_delete=CASCADE)
