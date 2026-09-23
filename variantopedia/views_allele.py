from dataclasses import dataclass
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.http import require_POST

from annotation.manual_variant_entry import check_can_create_variants
from annotation.models import (
    ClinVarRecordCollection,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from classification.enums import OverlapType
from classification.models import (
    OverlapStatus, Classification, Overlap,
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

    overlaps = Overlap.objects.filter(allele=allele, overlap_type=OverlapType.SINGLE_CONTEXT, valid=True, overlap_status__gte=OverlapStatus.SINGLE_SUBMITTER)
    overlaps = list(sorted(overlaps, key=lambda overlap: (overlap.testing_contexts_objs[0], overlap.value_type)))

    cross_overlaps = Overlap.objects.filter(allele=allele, overlap_type=OverlapType.CROSS_CONTEXT, valid=True, overlap_status__gte=OverlapStatus.SINGLE_SUBMITTER)
    cross_overlaps = list(sorted(cross_overlaps, key=lambda overlap: (overlap.testing_contexts_objs[0], overlap.value_type)))

    context = {
        # "allele_origin_groupings_desc": aogs,
        "overlaps": overlaps,
        "cross_overlaps": cross_overlaps,
        # "show_overall_diff": show_overall_diff,
        "allele_card": AlleleCard(user=request.user, allele=allele),
        "allele": allele,
        "allele_short_label": get_allele_short_label(allele, user_settings.default_genome_build),
        "edit_clinical_groupings": request.GET.get('edit_clinical_groupings') == 'True'
    }
    if request.user.is_superuser:
        withdrawn_count = Classification.objects.filter(allele_info__allele=allele, withdrawn=True).count()
        context["withdrawn_count"] = withdrawn_count

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
        # The user asked for this allele specifically, so retry every tool that has already failed on it
        create_liftover_pipelines(admin_bot(), [allele], ImportSource.WEB, variant_allele.genome_build, [genome_build],
                                  retry_conversion_tools=list(AlleleConversionTool))
    return redirect(allele)
