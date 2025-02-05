import json
import uuid
from datetime import timedelta
from html import escape
from typing import Union, Optional, Iterable, Any, Collection

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Model
from django.db.models.query import QuerySet
from django.template import Library
from django.utils.safestring import mark_safe
from django.utils.timezone import localtime

from classification.criteria_strengths import CriteriaStrength, AcmgPointScore
from classification.enums import SpecialEKeys
from classification.enums.classification_enums import ShareLevel
from classification.models import ConditionTextMatch, ConditionResolved, ClassificationLabSummary, ImportedAlleleInfo, \
    EvidenceMixin
from classification.models.classification import ClassificationModification, Classification
from classification.models.classification_groups import ClassificationGroup, ClassificationGroups, \
    ClassificationGroupUtils
from classification.models.classification_ref import ClassificationRef
from classification.models.clinical_context_models import ClinicalContext
from classification.models.discordance_models import DiscordanceReport
from classification.models.discordance_models_utils import DiscordanceReportRowData, DiscordanceReportTableData
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from classification.models.evidence_mixin import VCDbRefDict
from eventlog.models import ViewEvent
from genes.hgvs import CHGVS
from genes.models import GeneSymbol
from library.health_check import HealthCheckRequest
from ontology.models import OntologyTerm
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Lab
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.models.models_variant import Allele, Variant

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
        default_sort: Optional[str] = 'c_hgvs',
        allele_origin_filter_enabled: bool = True
    ):
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
    if isinstance(classification_modifications, QuerySet):
        classification_modifications = classification_modifications.select_related(
            'classification',
            'classification__clinical_context',
            'classification__lab',
            'classification__lab__organization'
        )

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
        "sort_order_index": sort_order_index,
        "allele_origin_filter_enabled": allele_origin_filter_enabled
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
        clinical_grouping_list = list({cm.classification.clinical_context for cm in ordered_classifications if cm.classification.clinical_context})
        clinical_grouping_list.sort(key=lambda cg:(not cg.is_default if cg else False, cg.name if cg else 'No Allele'))
        tag_context["clinical_contexts"] = clinical_grouping_list

    tag_context["paging"] = len(groups) > 10

    return tag_context


@register.inclusion_tag("classification/tags/classification_groupings.html", takes_context=True)
def classification_groupings(context, show_allele_origin_filter=True):
    """
    Shows the new database based classification grouping table. To filter the data implement a JavaScript method on the page
    <script>
        function classificationGroupingFilter(data) {
            data.ontology_term_id = {{ term.id | jsonify }};
        }
    </script>
    :param show_allele_origin_filter: True by default, set to False to hardcode the filtering to all records
    """
    return {"show_allele_origin_filter": show_allele_origin_filter, "genome_build": GenomeBuildManager.get_current_genome_build()}



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
def clinical_significance(value, evidence_key=SpecialEKeys.CLINICAL_SIGNIFICANCE, show_if_none=True):
    if isinstance(value, EvidenceMixin):
        value = value.get(evidence_key)
    if value is None and not show_if_none:
        return {"skip": True}

    key = EvidenceKeyMap.cached_key(evidence_key)
    label = key.option_dictionary.get(value, value) or "No Data"
    title = key.pretty_label
    if value == "withdrawn":
        label = "Withdrawn"

    prefix = "cs" if key.key == SpecialEKeys.CLINICAL_SIGNIFICANCE else "scs"
    css_value = value.lower() if value else "none"
    css_class = f"{prefix} {prefix}-{css_value}"

    return {
        "css_class": css_class,
        "label": label,
        "prefix": "cs" if key.key == SpecialEKeys.CLINICAL_SIGNIFICANCE else "scs",
        "title": title
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
def clinical_context(context, cc: ClinicalContext, orientation: str = 'horizontal'):
    return {"cc": cc, "orientation": orientation}


@register.inclusion_tag("classification/tags/classification_quick.html", takes_context=True)
def classification_quick(context,
                         vc: Union[Classification, ClassificationModification],
                         show_clinical_grouping=True,
                         show_lab=True,
                         show_condition=False,
                         show_flags=False,
                         record_count: Optional[int] = None, mode: Optional[str] = "detailed"):
    user = context.request.user
    vcm = vc
    if isinstance(vc, Classification):
        vcm = ClassificationModification.latest_for_user(user=user, classification=vc, published=True, exclude_withdrawn=False).first()
    return {
        "vcm": vcm,
        "show_clinical_grouping": show_clinical_grouping,
        "mode": mode,
        "show_lab": show_lab,
        "show_condition": show_condition,
        "show_flags": show_flags,
        "record_count": record_count
    }


class ClinicalGrouping:

    def __init__(self, cc: Optional[ClinicalContext]):
        self.cc = cc
        self.latest_report = DiscordanceReport.latest_report(cc) if cc else None
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


@register.inclusion_tag("classification/tags/allele.html")
def allele(allele: Allele):
    return {"allele": allele}


@register.inclusion_tag("classification/tags/allele_origin_toggle.html", takes_context=True)
def allele_origin_toggle(context, style: str = "", show_label: bool = True, prefill=None):
    """
    Embed the allele origin toggle on the page.
    It's expected this will only appear once per page and that there will be a JavaScript implementation of
    function alleleOriginToggle(value) {}
    on the page, where value can be "A", "G" or "S"
    :return:
    """
    user = context.request.user
    User_settings = UserSettings.get_for_user(user)
    value = prefill if prefill else User_settings.default_allele_origin
    return {"value": value, "style": style, "show_label": show_label}


@register.inclusion_tag("classification/tags/gene_symbol.html")
def gene_symbol(gene_symbol: GeneSymbol):
    return {"gene_symbol": gene_symbol}


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
        can_write = vc.can_write(user_or_group=user)

    c_hgvs = vc.c_hgvs_best(preferred_genome_build=genome_build)
    p_hgvs = None
    if settings.CLASSIFICATION_GRID_SHOW_PHGVS:
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
        "show_specimen_id": settings.CLASSIFICATION_SHOW_SPECIMEN_ID
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


@register.inclusion_tag("classification/tags/db_ref.html")
def db_ref(data: VCDbRefDict, css: Optional[str] = ''):
    context = dict(data)
    context['css'] = css
    return context


@register.inclusion_tag("classification/tags/condition.html")
def condition(condition_obj: Union[OntologyTerm, ConditionResolved],
              limit: Optional[int] = 100,
              show_link: Optional[bool] = True,
              no_condition_message: bool = False):
    if isinstance(condition_obj, OntologyTerm):
        condition_obj = ConditionResolved(terms=[condition_obj])
    return {"condition": condition_obj, "limit": limit, "show_link": show_link, "no_condition_message": no_condition_message}


# look at removing this
@register.inclusion_tag("classification/tags/discordance_report.html")
def discordance_report(discordance_report: DiscordanceReport):
    return {"discordance_report": discordance_report}


@register.inclusion_tag("classification/tags/discordance_report_row.html")
def discordance_report_row(
        discordance_report_summary: DiscordanceReportRowData,
        selected: Optional[DiscordanceReport] = None,
        filter: bool = False,
        discuss: bool = False,
        read_only: bool = False):
    return {
        "summary": discordance_report_summary,
        "filter": filter,
        "discuss": discuss,
        "is_selected": discordance_report_summary.discordance_report == selected,
        "read_only": read_only
    }


@register.inclusion_tag("classification/tags/discordance_report_table.html")
def discordance_report_table(table: DiscordanceReportTableData, filter: bool = False, discuss: bool = False, read_only: bool = False):
    return {
        "table": table,
        "filter": filter,
        "discuss": discuss,
        "read_only": read_only
    }


@register.inclusion_tag("classification/tags/classification_lab_summaries.html")
def classification_lab_summaries(
        lab_classification_summaries: Iterable[ClassificationLabSummary],
        shared: bool = True,
        include_acmg: bool = False,
        read_only: bool = False):

    return {
        "lab_classification_summaries": lab_classification_summaries,
        "shared": shared,
        "include_acmg": include_acmg,
        "read_only": read_only
    }

@register.inclusion_tag("classification/tags/criteria_strength.html")
def criteria_strength(strength: CriteriaStrength):
    return {"strength": strength}


@register.inclusion_tag("classification/tags/criteria_strength_td.html")
def criteria_strength_td(strength: Union[CriteriaStrength, Collection[CriteriaStrength]]):
    # going to display NM, NS, NA all the same
    if isinstance(strength, list):
        all_met_strengths = [s for s in strength if s.is_met]
        if len(all_met_strengths) > 1:
            return {
                "strengths": all_met_strengths
            }
        elif len(all_met_strengths) == 1:
            strength = all_met_strengths[0]
        else:
            strength = strength[0]

    return {
        "strength": strength
    }


@register.inclusion_tag("classification/tags/acmg_points.html")
def acmg_points(points: AcmgPointScore):
    return {"points": points}


@register.inclusion_tag("classification/tags/evidence_key_heading.html")
def evidence_key_heading(key: str, include_help: bool = True):
    e_key = EvidenceKeyMap.cached_key(key)
    return {
        "e_key": e_key,
        "help":  e_key.description if include_help else None
    }


@register.inclusion_tag("classification/tags/imported_allele_info.html")
def imported_allele_info(imported_allele_info: ImportedAlleleInfo, on_allele_page: bool = False):
    return {
        "imported_allele_info": imported_allele_info,
        "on_allele_page":  on_allele_page
    }


@register.simple_tag
def user_view_events(user: User, days: int = 1):
    if isinstance(user, str):
        try:
            user = User.objects.get(username=user)
        except User.DoesNotExist:
            return {}
    now = localtime()
    since = now - timedelta(days=days)
    health_request = HealthCheckRequest(since=since, now=now)
    view_events = ViewEvent.objects.filter(user=user, created__gte=health_request.since,
                                           created__lt=health_request.now).order_by('-created')
    view_event_data = []
    for event in view_events[:20]:
        view_event_data.append({
            'created': event.created,
            'view_name': event.view_name,
            'args': json.dumps(event.args)
        })
    return view_event_data