import uuid
from html import escape
from typing import Union, Optional, Iterable, Any, Collection

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Model
from django.db.models.query import QuerySet
from django.template import Library
from django.utils.safestring import mark_safe

from annotation.manual_variant_entry import check_can_create_variants, CreateManualVariantForbidden
from classification.criteria_strengths import CriteriaStrength
from classification.enums import SpecialEKeys
from classification.enums.classification_enums import ShareLevel
from classification.models import ConditionTextMatch, ConditionResolved, DiscordanceReportRowData, \
    ClassificationLabSummary
from classification.models.classification import ClassificationModification, Classification
from classification.models.classification_groups import ClassificationGroup, ClassificationGroups, \
    ClassificationGroupUtils
from classification.models.classification_ref import ClassificationRef
from classification.models.clinical_context_models import ClinicalContext
from classification.models.discordance_models import DiscordanceReport
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from classification.models.evidence_mixin import VCDbRefDict
from genes.hgvs import CHGVS
from genes.models import GeneSymbol
from library.utils import first, get_single_element
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import VariantAllele, Lab
from snpdb.models.models_genome import GenomeBuild, Contig, GenomeFasta
from snpdb.models.models_user_settings import UserSettings
from snpdb.models.models_variant import Allele, Variant, VariantAlleleSource
from snpdb.variant_links import variant_link_info
from uicore.templatetags.js_tags import jsonify

register = Library()


@register.inclusion_tag("classification/tags/condition_match.html")
def condition_match(condition_match: ConditionTextMatch, indent=0):
    return {
        "condition_match": condition_match,
        "indent": indent + 1,
        "indent_px": (indent + 1) * 16 + 8
    }


@register.inclusion_tag("classification/tags/classification_group_row.html")
def classification_group_row(group: ClassificationGroup, sub_row: Optional[int] = None, sub_index: Optional[int] = None, show_pending_changes: Optional[bool] = True):
    return {
        "group": group,
        "row_class": f"cc-{sub_row} collapse" if sub_row else "",
        "sub_index": sub_index,
        "show_username": settings.VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME,
        "show_allele_origin": settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN,
        "show_specimen_id": settings.VARIANT_CLASSIFICAITON_SHOW_SPECIMEN_ID
    }


@register.inclusion_tag("classification/tags/classification_groups.html", takes_context=True)
def classification_groups(
        context,
        classification_modifications: Iterable[ClassificationModification],
        show_diffs: bool = True,
        download_link: Optional[str] = None,
        history_link: Optional[str] = None,
        link_discordance_reports: bool = False,
        genome_build: Optional[GenomeBuild] = None,
        title: Optional[str] = None,
        context_object: Optional[Model] = None,
        group_utils: Optional[ClassificationGroupUtils] = None,
        default_sort: Optional[str] = 'c_hgvs'):
    """
    :param context: Auto included
    :param classification_modifications: The classification modifications to render
    :param show_diffs: Should a link to show diffs be shown
    :param download_link: URL to download this data
    :param history_link: URL to see the history of this data
    :param link_discordance_reports: Should link to discordance reports (if so will subdivide by clinical context)
    :param genome_build: Preferred genome build
    :param title: Heading to give the table
    :param context_object: If all these records are from an allele, provide "allele" if from a discordance report provide "discordance_report" etc
    :param old_classification_modifications: For showing what a discordance report used to be
    :param default_sort: The column to sort by default
    """
    sort_order_index = 1
    if default_sort == 'clinical_significance':
        sort_order_index = 2

    if not group_utils:
        group_utils = ClassificationGroupUtils(
            modifications=classification_modifications,
            calculate_pending=True
        )
    groups = ClassificationGroups(classification_modifications, genome_build=genome_build, group_utils=group_utils)

    tag_context = {
        "title": title,
        "classification_groups": groups,
        "user": context.request.user,
        "genome_build": groups.genome_build,
        "table_id": str(uuid.uuid4()).replace("-", "_"),
        "show_allele_origin": settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN,
        "sort_order_index": sort_order_index
    }
    ordered_classifications = list(groups.modifications)
    # classifications are sorted by group, display them so they're sorted by date
    ordered_classifications.sort(key=lambda cm: cm.curated_date_check, reverse=True)

    if groups and download_link:
        tag_context["download_link"] = download_link
    if groups and history_link and context.request.user.is_superuser:
        tag_context["history_link"] = history_link

    if show_diffs:
        if 1 < len(groups) <= 20 and len(groups) != len(ordered_classifications):
            diff_latest = ",".join([str(group.most_recent.classification.id) for group in groups])
            tag_context["diff_latest"] = diff_latest
        if 1 < len(ordered_classifications) <= 20:
            tag_context["diff_all"] = ",".join([str(cm.classification.id) for cm in ordered_classifications])

    tag_context["logging_key"] = ""
    if context_object:
        # in some contexts we get a string instead of Gene Symbol
        if isinstance(context_object, str):
            if gene_symbol := GeneSymbol.objects.filter(symbol=context_object).first():
                context_object = gene_symbol
        try:
            logging_key = context_object.metrics_logging_key
            tag_context["logging_key"] = f"&{logging_key[0]}={logging_key[1]}"
        except:
            raise ValueError(f"Context Object {context_object} does not have metrics_logging_key property")

    if link_discordance_reports:
        all_clinical_groupings = set()
        for cm in ordered_classifications:
            all_clinical_groupings.add(cm.classification.clinical_context)
        clinical_grouping_list = list(all_clinical_groupings)
        clinical_grouping_list.sort(key=lambda cg:(not cg.is_default if cg else False, cg.name if cg else 'No Allele'))
        tag_context["clinical_contexts"] = clinical_grouping_list

    tag_context["paging"] = len(groups) > 10

    return tag_context


def render_ekey(val, key: Optional[str] = None, value_if_none: Optional[str] = None):
    if isinstance(val, ClassificationModification):
        val = val.get(key)
    if not key:
        raise ValueError('ekey filter must have a key')
    e_key = EvidenceKeyMap.cached_key(key)
    if e_key.is_dummy:
        print(f"Warning, dummy evidence key {key}")
    pretty_val = e_key.pretty_value(val, dash_for_none=True)
    if val is None or val == '':
        if value_if_none is not None:
            return value_if_none
        return mark_safe(f'<span class="no-value">{escape(pretty_val)}</span>')
    elif (isinstance(val, list) or isinstance(val, set)) and len(val) == 0:
        if value_if_none is not None:
            return value_if_none
        return mark_safe(f'<span class="no-value">-</span>')
    return pretty_val


@register.filter
def ekey(val, key: Optional[str] = None):
    return render_ekey(val, key)


@register.filter
def ekey_raw(val, key: Optional[str] = None):
    """
    Renders blank if no value
    """
    return render_ekey(val, key, "")


@register.inclusion_tag("classification/tags/classification_history.html")
def classification_changes(changes):
    return {
        "changes": changes
    }


@register.inclusion_tag("classification/tags/clinical_significance.html")
def clinical_significance(value):
    key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    label = key.option_dictionary.get(value, value) or "Unclassified"
    if value == "withdrawn":
        label = "Withdrawn"

    return {
        "key": value.lower(),
        "label": label
    }

@register.inclusion_tag("classification/tags/clinical_significance_inline.html")
def clinical_significance_inline(value):
    key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    colors = {
        "B": "#44d",
        "LB": "#88d",
        "VUS": "#666",
        "LP": "#d88",
        "P": "#d44"
    }
    return {
        "key": value.lower(),
        "color": colors.get(value) or "#aaa",
        "label": key.option_dictionary.get(value, value) or "Unclassified"
    }

@register.inclusion_tag("classification/tags/lab.html")
def lab(lab: Lab, your_lab: Optional[Lab] = None):
    return {
        "lab": lab,
        "is_your_lab": your_lab is True or your_lab == lab
    }


@register.inclusion_tag("classification/tags/clinical_significance_select.html")
def clinical_significance_select(name, value):
    key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    return {
        "name": name,
        "options": key.virtual_options,
        "value": value
    }


@register.inclusion_tag("classification/tags/clinical_context.html", takes_context=True)
def clinical_context(context, cc: ClinicalContext, show_link: Optional[bool] = None):
    if show_link is None:
        show_link = context.request.user.is_superuser
    return {"cc": cc, "link": show_link}


@register.inclusion_tag("classification/tags/classification_quick.html", takes_context=True)
def classification_quick(context, vc: Union[Classification, ClassificationModification], show_clinical_grouping=True):
    user = context.request.user
    vcm = vc
    if isinstance(vc, Classification):
        vcm = ClassificationModification.latest_for_user(user=user, classification=vc, published=True, exclude_withdrawn=False).first()
    return {"vcm": vcm, "show_clinical_grouping": show_clinical_grouping}


class ClinicalGrouping:

    def __init__(self, cc: Optional[ClinicalContext]):
        self.cc = cc
        self.latest_report = DiscordanceReport.latest_report(cc)
        self.vcms = []

    @property
    def has_multiple(self):
        return len(self.vcms) > 1


@register.inclusion_tag("classification/tags/classification_table.html", takes_context=True)
def classification_table(
        context,
        records,
        genome_build=None,
        user=None,
        show_variant_link=False,
        show_clinical_context=False,
        variant: Variant = None,
        allele: Allele = None,
        edit_clinical_groupings=False):
    if user is None:
        user = context.request.user
    if not genome_build:
        if not user:
            raise ValueError('Must provide genome build or user')
        else:
            genome_build = UserSettings.get_for_user(user).default_genome_build

    if isinstance(records, QuerySet):
        records = list(records.all())
    mods = []
    for r in records:
        if isinstance(r, Classification):
            r = r.last_published_version
        if r:
            mods.append(r)

    # group by cc if we're doing clinical context
    if show_clinical_context:
        groupings = {}
        no_cc_grouping = ClinicalGrouping(cc=None)
        for vcm in mods:
            cc = vcm.classification.clinical_context
            if cc and vcm.share_level in ShareLevel.DISCORDANT_LEVEL_KEYS:
                grouping = groupings.get(cc.id)
                if not grouping:
                    grouping = ClinicalGrouping(cc=cc)
                    groupings[cc.id] = grouping
                grouping.vcms.append(vcm)
            else:
                no_cc_grouping.vcms.append(vcm)

        records = list(groupings.values())
        if no_cc_grouping.vcms:
            records.append(no_cc_grouping)
    else:
        not_grouped = ClinicalGrouping(cc=None)
        not_grouped.vcms = mods
        records = [not_grouped]

    # now to sort groups, and then sort vcms inside groups
    def clinical_group_sort_score(cg: ClinicalGrouping):
        cc = cg.cc
        if not cc:
            return 'zzzz'
        else:
            return cc.name

    cs_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    def classification_modification_sort_score(vcm: ClassificationModification):
        cs = vcm.evidence.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        options = cs_key.matched_options(cs)
        if options:
            return options[0].get('index', 0), vcm.id_str
        else:
            return None, vcm.id_str

    records.sort(key=clinical_group_sort_score)
    for record in records:
        record.vcms.sort(key=classification_modification_sort_score)

    return {
        "records": records,
        "variant": variant,
        "allele": allele,
        "genome_build": genome_build,
        "show_variant_link": show_variant_link,
        "show_clinical_context": show_clinical_context,
        "edit_clinical_groupings": edit_clinical_groupings,
        "show_allele_origin": settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN,
        "user": user,
        "discordance_enabled": settings.DISCORDANCE_ENABLED
    }


def _to_c_hgvs(c_hgvs: Any) -> CHGVS:
    if isinstance(c_hgvs, ClassificationModification):
        if c_hgvs := c_hgvs.classification.get_c_hgvs(GenomeBuildManager.get_current_genome_build()):
            c_hgvs = CHGVS(c_hgvs)
            c_hgvs.genome_build = GenomeBuildManager.get_current_genome_build()
    elif isinstance(c_hgvs, str):
        c_hgvs = CHGVS(c_hgvs)

    if c_hgvs is None:  # might have got a none c.hgvs from the ClassificationModification
        c_hgvs = CHGVS("")

    return c_hgvs


@register.inclusion_tag("classification/tags/c_hgvs.html")
def c_hgvs(c_hgvs: Union[CHGVS, ClassificationModification, str], show_genome_build: Optional[bool] = None):
    c_hgvs = _to_c_hgvs(c_hgvs)

    if show_genome_build is None:
        show_genome_build = c_hgvs.is_desired_build is False or c_hgvs.is_normalised is False

    return {"c_hgvs": c_hgvs, "show_genome_build": show_genome_build and c_hgvs.genome_build is not None}


@register.inclusion_tag("classification/tags/classification_row.html", takes_context=True)
def classification_row(
        context,
        record: Union[Classification, ClassificationModification],
        genome_build: GenomeBuild,
        user: User = None,
        show_variant_link=False,
        show_clinical_context=False,
        edit_clinical_context=False):
    if user is None:
        user = context.request.user

    vc = record
    vcm = None
    if isinstance(record, ClassificationModification):
        vc = record.classification
        vcm = record
    else:
        vcm = record.last_published_version
    icon = 'icons/share_level/' + vc.share_level_enum.key + '.png'

    try:
        curated = Classification.to_date(record.get(SpecialEKeys.CURATION_DATE, None))
    except ValueError:
        curated = None

    can_write = False
    if user:
        can_write = vc.can_write(user=user)

    c_hgvs = vc.c_hgvs_best(preferred_genome_build=genome_build)
    p_hgvs = None
    if settings.VARIANT_CLASSIFICATION_GRID_SHOW_PHGVS:
        p_hgvs = record.get(SpecialEKeys.P_HGVS)
        if p_hgvs:
            p_dot = p_hgvs.find('p.')
            if p_dot != -1:
                p_hgvs = p_hgvs[p_dot::]

    return {
        "evidence": record.evidence,
        "condition_obj": vc.condition_resolution_obj,
        "curated": curated,
        "c_hgvs": c_hgvs,
        "gene_symbol": vcm.get(SpecialEKeys.GENE_SYMBOL),
        "vc": vc,
        "vcm": vcm,
        "icon": icon,
        "p_hgvs": p_hgvs,
        "can_write": can_write,
        "show_variant_link": show_variant_link,
        "show_clinical_context": show_clinical_context,
        "edit_clinical_context": edit_clinical_context,
        "show_allele_origin": settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN,
        "show_specimen_id": settings.VARIANT_CLASSIFICAITON_SHOW_SPECIMEN_ID
    }


@register.inclusion_tag("classification/tags/classification.html")
def classification(classification, user):
    ref = ClassificationRef.init_from_obj(user=user, obj=classification)
    record = ref.record
    return {
        "lab_name": record.lab.name,
        "lab_record_id": record.lab_record_id
    }


@register.filter
def option_label(value):
    label = value.get('label')
    if label:
        return label
    else:
        return EvidenceKey.pretty_label_from_string(value.get('key'))


@register.filter
def classification_count(obj: Allele) -> int:
    if isinstance(obj, Allele):
        return Classification.objects.filter(variant__in=obj.variants).count()
    elif isinstance(obj, Variant):
        return Classification.objects.filter(variant=obj).count()
    else:
        return 0


@register.inclusion_tag("classification/tags/variant_card.html", takes_context=True)
def variant_card(context, allele: Allele, genome_build: GenomeBuild):
    request = context.request
    can_create_classification = Classification.can_create_via_web_form(request.user)
    va: VariantAllele = allele.variant_alleles().filter(genome_build=genome_build).first()
    liftover_error_qs = allele.liftovererror_set.filter(liftover__genome_build=genome_build)

    unfinished_liftover = None
    can_create_variant = False
    if va is None:
        unfinished_liftover = VariantAlleleSource.get_liftover_for_allele(allele, genome_build)
        if unfinished_liftover is None:
            try:
                check_can_create_variants(request.user)
                try:
                    # See if we can have data already to liftover
                    conversion_tool, _ = allele.get_liftover_tuple(genome_build)
                    can_create_variant = conversion_tool is not None
                except (Contig.ContigNotInBuildError, GenomeFasta.ContigNotInFastaError):
                    pass
            except CreateManualVariantForbidden:
                pass

    return {
        "user": request.user,
        "allele": allele,
        "unfinished_liftover": unfinished_liftover,
        "genome_build": genome_build,
        "variant_allele": va,
        "liftover_error_qs": liftover_error_qs,
        "can_create_classification": can_create_classification,
        "can_create_variant": can_create_variant,
    }


@register.filter
def quick_link_data(variant_allele: VariantAllele):
    """ Needs to be VariantAllele as we need genome_build too """
    data = variant_link_info(variant_allele.variant, variant_allele.genome_build)
    return jsonify(data)


@register.inclusion_tag("classification/tags/db_ref.html")
def db_ref(data: VCDbRefDict, css: Optional[str] = ''):
    context = dict(data)
    context['css'] = css
    return context


@register.inclusion_tag("classification/tags/condition.html")
def condition(condition_obj: ConditionResolved, limit: Optional[int] = 100):
    return {"condition": condition_obj, "limit": limit}


# look at removing this
@register.inclusion_tag("classification/tags/discordance_report.html")
def discordance_report(discordance_report: DiscordanceReport):
    return {"discordance_report": discordance_report}


@register.inclusion_tag("classification/tags/discordance_report_row.html")
def discordance_report_row(discordance_report_summary: DiscordanceReportRowData, selected: Optional[DiscordanceReport] = None, filter: bool = False):
    return {"summary": discordance_report_summary, "filter": filter, "is_selected": discordance_report_summary.discordance_report == selected}


@register.inclusion_tag("classification/tags/classification_lab_summaries.html")
def classification_lab_summaries(lab_classification_summaries: Iterable[ClassificationLabSummary], shared: bool = True, include_acmg: bool = False):
    return {
        "lab_classification_summaries": lab_classification_summaries,
        "shared": shared,
        "include_acmg": include_acmg
    }

@register.inclusion_tag("classification/tags/criteria_strength.html")
def criteria_strength(strength: CriteriaStrength):
    return {"strength": strength}


@register.inclusion_tag("classification/tags/criteria_strengths.html")
def criteria_strengths(strengths: Collection[CriteriaStrength]):
    # going to display NM, NS, NA all the same
    has_different_points = set(strength.strength if strength.is_met else 'NM' for strength in strengths)
    if len(has_different_points) == 1:
        strengths = [next(iter(strengths))]
    return {
        "strengths": strengths,
        "single_not_met": len(strengths) == 1 and strengths[0].strength_direction == "N"
    }