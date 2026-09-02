from dataclasses import dataclass
from typing import Optional

from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.http import require_POST

from annotation.manual_variant_entry import check_can_create_variants
from annotation.models import (
    ClassificationModification,
    ClinVarRecordCollection,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from classification.enums import AlleleOriginBucket, ShareLevel, SpecialEKeys
from classification.models import (
    AlleleOriginGrouping,
    ClassificationGrouping,
    ClassificationGroupingEntry,
    DiscordanceReport,
    EvidenceKeyMap,
    OverlapStatus,
)
from classification.variant_card import AlleleCard
from classification.views.exports import ClassificationExportFormatterCSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_csv import FormatDetailsCSV
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import link_allele_to_existing_variants
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import Allele, AlleleConversionTool, AlleleOrigin, ImportSource
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.utils import get_genome_build_or_404


@dataclass
class ShareLevelRecordCounts:
    lab_count: int
    # record_count: int


@dataclass
class AlleleOriginGroupingDescription:
    allele_origin_grouping: AlleleOriginGrouping
    discordance_report: Optional[DiscordanceReport]
    overlap_status: OverlapStatus
    shared_counts: int
    unshared_counts: int

    @property
    def get_overlap_status_display(self):
        return OverlapStatus(self.overlap_status).label

    @property
    def should_show_diffs(self):
        return self.shared_counts + self.unshared_counts > 1

    @staticmethod
    def describe(allele_origin_grouping: AlleleOriginGrouping, for_user: User) -> 'AlleleOriginGroupingDescription':
        discordance_report: Optional[DiscordanceReport] = None

        allele = allele_origin_grouping.allele_grouping.allele
        if allele_origin_grouping.allele_origin_bucket == AlleleOriginBucket.GERMLINE:
            for cc in allele.clinicalcontext_set.all():
                if dr := DiscordanceReport.latest_report(cc):
                    if dr.is_active:
                        discordance_report = dr

        shared_counts = allele_origin_grouping.classificationgrouping_set.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).count()
        unshared_counts = (
            ClassificationGrouping.filter_for_user(for_user, allele_origin_grouping.classificationgrouping_set.exclude(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)).count())

        overlap_status: OverlapStatus
        if shared_counts == 0:
            overlap_status = OverlapStatus.NO_SHARED_RECORDS
        elif shared_counts == 1:
            overlap_status = OverlapStatus.SINGLE_SUBMITTER
        elif allele_origin_grouping.allele_origin_bucket != AlleleOriginBucket.GERMLINE:
            overlap_status = OverlapStatus.NOT_COMPARABLE_OVERLAP
        else:
            # TODO right now it's lookiung at every classification in the grouping
            # In future, only look at the latest
            classification_groups = allele_origin_grouping.classificationgrouping_set.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).values_list('pk', flat=True)
            classification_groups_classes = ClassificationGroupingEntry.objects.filter(grouping__in=classification_groups).values_list('classification_id', flat=True)
            classification_values = set(ClassificationModification.objects.filter(is_last_published=True, classification_id__in=classification_groups_classes).values_list(f'published_evidence__{SpecialEKeys.CLINICAL_SIGNIFICANCE}__value', flat=True).all())

            # now see if we're agreement, confidence or discordant
            bucket_mapping = EvidenceKeyMap.instance().get(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("bucket")
            buckets = {bucket_mapping.get(class_value) for class_value in classification_values}
            if None in buckets:
                buckets.remove(None)

            if len(buckets) > 1:
                # discordant
                if "P" in classification_values or "LP" in classification_values:
                    overlap_status = OverlapStatus.DISCORDANCE_MEDICALLY_SIGNIFICANT
                else:
                    overlap_status = OverlapStatus.DISCORDANCE
            else:
                if len(classification_values) > 1:
                    overlap_status = OverlapStatus.CONFIDENCE
                else:
                    # complete agreement
                    overlap_status = OverlapStatus.AGREEMENT

        return AlleleOriginGroupingDescription(
            allele_origin_grouping=allele_origin_grouping,
            discordance_report=discordance_report,
            overlap_status=overlap_status,
            shared_counts=shared_counts,
            unshared_counts=unshared_counts,
        )


def get_allele_short_label(allele: Allele, preferred_genome_build: GenomeBuild) -> str:
    """ Shortest label that still tells a human which allele this is, eg "RUNX1:p.Ala547Val" - taken from
        the representative (MANE where we have it) transcript, falling back to g.HGVS """
    variant_alleles = list(allele.variant_alleles())
    preferred_first = sorted(variant_alleles, key=lambda va: va.genome_build != preferred_genome_build)
    for variant_allele in preferred_first:
        vav = VariantAnnotationVersion.latest(variant_allele.genome_build)
        if va := variant_allele.variant.variantannotation_set.filter(version=vav).first():
            return va.get_short_label()
    if preferred_first:
        return VariantAnnotation.get_hgvs_g(preferred_first[0].variant) or str(allele)
    return str(allele)


def view_allele(request, allele_id: int):
    allele: Allele = get_object_or_404(Allele, pk=allele_id)
    user_settings = UserSettings.get_for_user(request.user)
    link_allele_to_existing_variants(allele, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
    ClinVarRecordCollection.set_allele_for_variants(allele)

    # Filter on classification grouping first, so we can find all unique AlleleGroupings
    # that the user has access to
    aog_qs = AlleleOriginGrouping.objects.filter(pk__in=\
        ClassificationGrouping.filter_for_user(
            request.user,
            ClassificationGrouping.objects.filter(allele_origin_grouping__allele_grouping__allele=allele_id)
        ).values_list("allele_origin_grouping")
    )
    aogs = [AlleleOriginGroupingDescription.describe(aog, request.user) for aog in sorted(aog_qs.all())]

    show_overall_diff = len(aogs) > 1

    context = {
        "allele_origin_groupings_desc": aogs,
        "show_overall_diff": show_overall_diff,
        "allele_card": AlleleCard(user=request.user, allele=allele),
        "allele": allele,
        "allele_short_label": get_allele_short_label(allele, user_settings.default_genome_build),
        "edit_clinical_groupings": request.GET.get('edit_clinical_groupings') == 'True'
    }
    return render(request, "variantopedia/view_allele.html", context)


def export_classifications_allele(request, allele_id: int):
    """
    CSV export of what is currently filtered into the classification grid
    """
    allele = get_object_or_404(Allele, pk=allele_id)
    return ClassificationExportFormatterCSV(
        ClassificationFilter(
            user=request.user,
            genome_build=GenomeBuildManager.get_current_genome_build(),
            allele=allele_id,
            file_prefix=f"classifications_allele_{allele:CA}"
        ),
        FormatDetailsCSV()
    ).serve()


@require_POST
def create_variant_for_allele(request, allele_id, genome_build_name):
    """ Shortcut to create manual variant, but as a POST """
    check_can_create_variants(request.user)
    allele = get_object_or_404(Allele, pk=allele_id)
    genome_build = get_genome_build_or_404(genome_build_name)
    non_liftover_origin = [AlleleOrigin.IMPORTED_TO_DATABASE, AlleleOrigin.IMPORTED_NORMALIZED]
    if variant_allele := allele.variantallele_set.filter(origin__in=non_liftover_origin).first():
        create_liftover_pipelines(admin_bot(), [allele], ImportSource.WEB, variant_allele.genome_build, [genome_build])
    return redirect(allele)
