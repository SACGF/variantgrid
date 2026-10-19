"""
Extra info for the FlagCollections attached to Classifications (label, variant, clinical context, write access and
any ongoing discordance report), answering flag_collection_extra_info_signal for the flags pages.
"""
from django.contrib.auth.models import User
from django.dispatch.dispatcher import receiver

from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models.classification import Classification
from classification.models.discordance_models import DiscordanceReportClassification
from flags.models.models import FlagCollection, FlagInfos, flag_collection_extra_info_signal


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs) -> None:  # pylint: disable=unused-argument
    """
    Allows us to provide extra info for FlagCollections attached to Classification.
    e.g. linking to the appropriate allele page, discordance report etc.
    :param flag_infos: Information on the flag collections being displayed to the user.
    Populates this with the extra info
    :param user: The current user
    :param kwargs: Required by @receiver
    """
    vcs = Classification.objects.filter(flag_collection__in=flag_infos.ids).select_related('lab')
    drcs = DiscordanceReportClassification.objects.filter(classification_original__classification__in=vcs,
                                                          report__resolution=DiscordanceReportResolution.ONGOING)
    drcs_dict = {}
    for drc in drcs.values_list('classification_original__classification', 'report'):
        drcs_dict[drc[0]] = drc[1]

    for vc in vcs:
        """
            # If we are to use sub-flags in their current incarnation, here's how we'd include
            # one if zygosity or condition are missing
            flags_for_vc = flag_infos.flags_for_collection[vc.flag_collection_id]
            unsubmitted_flag = find_unsubmitted(flags_for_vc)

            if unsubmitted_flag:
                if vc.get(SpecialEKeys.ZYGOSITY) is None:
                    flag_infos.record_sub_flag(unsubmitted_flag, 'Missing Zygosity', 'Z')
                if vc.get(SpecialEKeys.CONDITION) is None:
                    flag_infos.record_sub_flag(unsubmitted_flag, 'Missing Condition', 'C')
        """

        context = {
            'label': vc.friendly_label,
            'vc_id': vc.id
        }
        if vc.variant_id:
            context['variant'] = vc.variant_id
        if vc.clinical_context_id:
            context['clinical_context'] = vc.clinical_context_id
        if user and vc.can_write(user):
            context['can_write'] = True

        dc_report_id = drcs_dict.get(vc.id)
        if dc_report_id:
            context['discordance_report'] = dc_report_id

        flag_infos.set_extra_info(vc.flag_collection_id, context, source_object=vc)
