from html import escape

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models.query import QuerySet
from django.template import Library
from typing import Union, Optional

from django.utils.safestring import mark_safe

from snpdb.models import VariantAllele
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.models.models_variant import Allele, Variant
from snpdb.variant_links import variant_link_info
from classification.enums import SpecialEKeys
from classification.enums.variant_classification_enums import ShareLevel
from classification.models import BestHGVS
from classification.models.clinical_context_models import ClinicalContext
from classification.models.discordance_models import DiscordanceReport, \
    DiscordanceReportClassification
from classification.models.evidence_key import EvidenceKey, \
    EvidenceKeyMap
from classification.models.variant_classification import VariantClassificationModification, \
    VariantClassification
from classification.models.variant_classification_ref import VariantClassificationRef
from classification.templatetags.js_tags import jsonify

register = Library()


@register.filter
def ekey(val, key: str = None):
    if not key:
        raise ValueError('ekey filter must have a key')
    e_key = EvidenceKeyMap.cached_key(key)
    pretty_val = e_key.pretty_value(val, dash_for_none=True)
    if val is None or val == '':
        return mark_safe(f'<span class="no-value">{escape(pretty_val)}</span>')
    return pretty_val


@register.inclusion_tag("classification/tags/variant_classification_history.html")
def variant_classification_changes(changes):
    return {
        "changes": changes
    }

@register.inclusion_tag("classification/tags/clinical_significance.html")
def clinical_significance(value):
    key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    return {
        "key": value.lower(),
        "label": key.option_dictionary.get(value, value) or "Unclassified"
    }


@register.inclusion_tag("classification/tags/clinical_significance_select.html")
def clinical_significance_select(name, value):
    key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    return {
        "name": name,
        "options": key.virtual_options,
        "value": value
    }

@register.inclusion_tag("classification/tags/clinical_context.html")
def clinical_context(cc: ClinicalContext, user: User):
    return {"cc": cc, "link": user.is_superuser}


@register.inclusion_tag("classification/tags/classification_quick.html", takes_context=True)
def classification_quick(context, vc: Union[VariantClassification, VariantClassificationModification]):
    user = context.request.user
    vcm = vc
    if isinstance(vc, VariantClassification):
        vcm = VariantClassificationModification.latest_for_user(user=user, variant_classification=vc, published=True, exclude_withdrawn=False).first()
    return {"vcm": vcm}


class ClinicalGrouping:

    def __init__(self, cc: Optional[ClinicalContext]):
        self.cc = cc
        self.latest_report = DiscordanceReport.latest_report(cc)
        self.vcms = []

    @property
    def has_multiple(self):
        return len(self.vcms) > 1


@register.inclusion_tag("classification/tags/variant_classification_table.html", takes_context=True)
def variant_classification_table(
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
        if isinstance(r, VariantClassification):
            r = r.last_published_version
        if r:
            mods.append(r)

    # group by cc if we're doing clinical context
    if show_clinical_context:
        groupings = {}
        no_cc_grouping = ClinicalGrouping(cc=None)
        for vcm in mods:
            cc = vcm.variant_classification.clinical_context
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

    def variant_classification_modification_sort_score(vcm: VariantClassificationModification):
        cs = vcm.evidence.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        options = cs_key.matched_options(cs)
        if options:
            return options[0].get('index', 0), vcm.id_str
        else:
            return None, vcm.id_str

    records.sort(key=clinical_group_sort_score)
    for record in records:
        record.vcms.sort(key=variant_classification_modification_sort_score)

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


@register.inclusion_tag("classification/tags/hgvs.html", takes_context=True)
def hgvs(context, hgvs: BestHGVS, show_variant_link: bool = True):
    return {"hgvs": hgvs, "show_variant_link": show_variant_link}


@register.inclusion_tag("classification/tags/variant_classification_row.html", takes_context=True)
def variant_classification_row(
        context,
        record: Union[VariantClassification, VariantClassificationModification],
        genome_build: GenomeBuild,
        user: User = None,
        show_variant_link=False,
        show_clinical_context=False,
        edit_clinical_context=False):
    if user is None:
        user = context.request.user

    vc = record
    vcm = None
    if isinstance(record, VariantClassificationModification):
        vc = record.variant_classification
        vcm = record
    else:
        vcm = record.last_published_version
    icon = 'icons/share_level/' + vc.share_level_enum.key + '.png'

    try:
        curated = VariantClassification.to_date(record.get(SpecialEKeys.CURATION_DATE, None))
    except ValueError:
        curated = None

    can_write = False
    if user:
        can_write = vc.can_write(user=user)

    best_hgvs = vc.best_hgvs(genome_build)
    p_hgvs = None
    if settings.VARIANT_CLASSIFICATION_GRID_SHOW_PHGVS:
        p_hgvs = record.get(SpecialEKeys.P_HGVS)
        if p_hgvs:
            p_dot = p_hgvs.find('p.')
            if p_dot != -1:
                p_hgvs = p_hgvs[p_dot::]

    return {
        "evidence": record.evidence,
        "curated": curated,
        "best_hgvs": best_hgvs,
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


@register.inclusion_tag("classification/tags/variant_classification.html")
def variant_classification(classification, user):
    ref = VariantClassificationRef.init_from_obj(user=user, obj=classification)
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
def variant_classification_count(obj: Allele) -> int:
    if isinstance(obj, Allele):
        return VariantClassification.objects.filter(variant__in=obj.variants).count()
    elif isinstance(obj, Variant):
        return VariantClassification.objects.filter(variant=obj).count()
    else:
        return 0


@register.inclusion_tag("classification/tags/variant_classification_discordance_row.html")
def variant_classification_discordance_row(row: DiscordanceReportClassification, show_flags=False):
    vc = row.classification_original.variant_classification
    icon = 'icons/share_level/' + vc.share_level_enum.key + '.png'
    return {
        "vc": vc,
        "icon": icon,
        "action_log": row.action_log,
        "best_hgvs": row.classfication_effective.get(SpecialEKeys.C_HGVS, None),
        "starting": row.classification_original,
        "closing": row.classfication_effective,
        "starting_curated": row.classification_original.get(SpecialEKeys.CURATION_DATE, None),
        "closing_curated": row.classfication_effective.get(SpecialEKeys.CURATION_DATE, None),
        "show_flags": show_flags,
    }


@register.inclusion_tag("classification/tags/variant_card.html", takes_context=True)
def variant_card(context, allele: Allele, build: GenomeBuild):
    request = context.request
    can_create_classification = VariantClassification.can_create_via_web_form(request.user)
    va: VariantAllele = allele.variant_alleles().filter(genome_build=build).first()
    liftover_error_qs = allele.liftovererror_set.filter(liftover__genome_build=build)

    return {
        "user": request.user,
        "genome_build": build,
        "variant_allele": va,
        "can_create_classification": can_create_classification,
        "liftover_error_qs": liftover_error_qs,
    }


@register.filter
def quick_link_data(variant: Variant):
    data = variant_link_info(variant)
    return jsonify(data)
