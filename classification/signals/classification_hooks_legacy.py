# from dataclasses import dataclass
# from typing import Dict, Any, Optional, Set
#
# import django
# from django.dispatch import receiver
# from classification.enums import SpecialEKeys
# from classification.models import Classification, ImportedAlleleInfoStatus, classification_flag_types, \
#     classification_variant_set_signal, variants_classification_changed_signal, ImportedAlleleInfo
# from flags.models import FlagStatus, FlagComment, Flag
# from genes.hgvs import HGVSMatcher, CHGVS, CHGVSDiff, chgvs_diff_description
# from library.guardian_utils import admin_bot
# from library.log_utils import report_exc_info
# from snpdb.models import Allele, GenomeBuild, allele_flag_types
#
# """
# IMPORTANT!!
#
# NONE OF THIS CODE SHOULD BE HOOKED UP TO ANYTHING OR RUN ANYMORE.
# Will delete this soon one everything has been migrated
# """
#
# allele_validate_signal = django.dispatch.Signal()  # args: "allele"
# # we no longer validate alleles
#
# ### LEGACY ALLELE STUFF
#
# class AlleleLegacy(Allele):
#
#     def validate(self, liftover_complete: bool):
#         # this used to be Allele.validate()
#         if liftover_complete:
#             v37 = self.variant_alleles().filter(genome_build=GenomeBuild.grch37()).first()
#             v38 = self.variant_alleles().filter(genome_build=GenomeBuild.grch38()).first()
#
#             if v37:
#                 self.close_open_flags_of_type(allele_flag_types.missing_37)
#             else:
#                 self.flag_collection_safe.get_or_create_open_flag_of_type(flag_type=allele_flag_types.missing_37,
#                                                                           only_if_new=True)
#
#             if v38:
#                 self.close_open_flags_of_type(allele_flag_types.missing_38)
#             else:
#                 self.flag_collection_safe.get_or_create_open_flag_of_type(flag_type=allele_flag_types.missing_38,
#                                                                           only_if_new=True)
#
#         allele_validate_signal.send(sender=Allele, allele=self)
#
#
# @dataclass(frozen=True)
# class TranscriptDifference:
#     transcript: str
#     chgvs37: str
#     chgvs38: str
#
#     @property
#     def comment(self) -> str:
#         return (
#             f'Attached classification with transcript {self.transcript} appears as the following in \n\n'
#             f'{self.chgvs37} (GRCh37)\n'
#             f'{self.chgvs38} (GRCh38)')
#
#     @property
#     def data(self) -> Dict[str, Any()]:
#         return {
#             "transcript": self.transcript,
#             'chgvs37': self.chgvs37,
#             'chgvs38': self.chgvs38
#         }
#
#     @staticmethod
#     def from_flag(flag: Flag) -> 'TranscriptDifference':
#         data = flag.data or {}
#         return TranscriptDifference(
#             transcript=data.get('transcript'),
#             chgvs37=data.get('chgvs37'),
#             chgvs38=data.get('chgvs38')
#         )
#
#
# @receiver(allele_validate_signal, sender=Allele)
def compare_chgvs(sender, allele: Allele, **kwargs):  # pylint: disable=unused-argument
    # FIXME remove this code, no longer required compared to AlleleInfo
    transcript_differences: Set[TranscriptDifference] = set()
    vcs = Classification.objects.filter(variant__in=allele.variants).order_by('id')
    v37 = allele.grch37
    v38 = allele.grch38
    if v37 and v38:
        matcher37 = HGVSMatcher(genome_build=GenomeBuild.grch37())
        matcher38 = HGVSMatcher(genome_build=GenomeBuild.grch38())
        classification_transcripts = set()
        vc: Classification
        for vc in vcs:
            if transcript := vc.transcript:
                classification_transcripts.add(transcript)

        for transcript in classification_transcripts:
            chgvs37 = None
            chgvs38 = None
            try:
                chgvs37 = matcher37.variant_to_c_hgvs(v37, transcript)
            except ValueError as ve:
                chgvs37 = f'Error: {str(ve)}'

            try:
                chgvs38 = matcher38.variant_to_c_hgvs(v38, transcript)
            except ValueError as ve:
                chgvs38 = f'Error: {str(ve)}'

            if chgvs37 != chgvs38:
                transcript_differences.add(TranscriptDifference(transcript=transcript, chgvs37=chgvs37, chgvs38=chgvs38))

        # loop through existing flags and close them if we don't have a mismatch for them anymore
        # or re-open them if we do (if flag was closed by admin_bot)
        existing_flag: Flag
        for existing_flag in list(allele.flag_collection_safe.flag_set.filter(flag_type=allele_flag_types.allele_37_not_38)):
            existing_td = TranscriptDifference.from_flag(existing_flag)

            if existing_td in transcript_differences:
                # current mismatch already has a flag
                transcript_differences.remove(existing_td)

                # flag is closed, but it was closed by a bot, so we can re-open it
                if existing_flag.resolution.status != FlagStatus.OPEN and FlagComment.last(existing_flag).user == admin_bot():
                    existing_flag.flag_action(comment="A classification with this transcript, GRCh37 and GRCh38 has been re-added or resolved", resolution=existing_flag.flag_type.resolution_for_status(FlagStatus.OPEN))
                # otherwise leave open flag (or closed by a human flag) alone

            elif existing_flag.resolution.status == FlagStatus.OPEN:
                existing_flag.flag_action(comment="No classifications with this transcript, GRCh37 and GRCh38 combination remaining", resolution=existing_flag.flag_type.resolution_for_status(FlagStatus.CLOSED))

        for new_issue_td in transcript_differences:
            # there shouldn't be any other flags with this data based on above checks
            allele.flag_collection_safe.get_or_create_open_flag_of_type(
                flag_type=allele_flag_types.allele_37_not_38,
                comment=new_issue_td.comment,
                data=new_issue_td.data)

    else:
        # if there's no 37 or no 38, close any flag comparing the two
        allele.close_open_flags_of_type(allele_flag_types.allele_37_not_38, comment="Lacking representation in both 37 and 38")

#
# ### LEGACY CLASSIFICATION STUFF
#
# class LegacyClassification(Classification):
#
#     def allele_classification_changed(self):
#         """ Notifies all variants linked to allele that the classification has changed """
#         if allele := self.allele_info.allele:
#             for variant_allele in allele.variant_alleles():
#                 variants_classification_changed_signal.send(sender=Classification,
#                                                             variants=[variant_allele.variant],
#                                                             genome_build=variant_allele.genome_build)
#
#
#     # was called on Classification when variant was set
#     def variant_set(self, allele_or_variant_changed: bool):
#         """
#         No longer called, this is how Classifications updated their flags when variants were set
#         but we don't use those flags anymore
#         :param allele_or_variant_changed:
#         :return:
#         """
#
#         allele_info = self.allele_info
#         variant = allele_info.matched_variant
#
#         failed = allele_info.status == ImportedAlleleInfoStatus.FAILED
#         message = allele_info.message
#
#         # Classifications no longer have classification_import, that's in ImportedAlleleInfo
#         # # don't want to be considered as part of the import anymore
#         # # as we've failed matching somewhere along the line
#         # if variant is None and failed:
#         #     self.classification_import = None
#
#         # see if the clinical context is still relevant (if we have one)
#         if not variant or ((cc := self.clinical_context) and (variant.allele != cc.allele)):
#             self.clinical_context = None
#
#         if not self.id:
#             # need to save to have a flag collection
#             self.save()
#
#         flag_collection = self.flag_collection_safe
#         try:
#             if failed:
#                 # failed matching
#                 if not message:
#                     c_hgvs = self.get(SpecialEKeys.C_HGVS)
#                     build_name = self.get(SpecialEKeys.GENOME_BUILD)
#                     message = f'Could not resolve {build_name} {c_hgvs}'
#
#                 flag_collection.ensure_resolution(classification_flag_types.matching_variant_flag,
#                                                   resolution='matching_failed',
#                                                   comment=message)
#             elif variant:
#                 # matching success
#                 flag_collection.close_open_flags_of_type(classification_flag_types.matching_variant_flag,
#                                                          comment='Variant Matched')
#                 classification_variant_set_signal.send(sender=Classification, classification=self, variant=variant)
#
#                 # FIXME, just move this to allele info
#                 if allele_or_variant_changed:
#                     self.allele_classification_changed()
#             else:
#                 # matching ongoing
#                 flag_collection.ensure_resolution(classification_flag_types.matching_variant_flag,
#                                                   resolution='open',
#                                                   comment=message)
#                 flag_collection.close_open_flags_of_type(classification_flag_types.matching_variant_warning_flag,
#                                                          comment='Variant Re-Matching')
#                 flag_collection.close_open_flags_of_type(classification_flag_types.transcript_version_change_flag,
#                                                          comment='Variant Re-Matching')
#
#             self._perform_c_hgvs_validation()
#
#         except:
#             report_exc_info()
#             flag_collection.ensure_resolution(classification_flag_types.matching_variant_flag,
#                                               resolution='matching_failed',
#                                               comment='Could not set variant for unexpected reason')
#
#     def _perform_c_hgvs_validation(self):
#         # if we had a previously opened flag match warning - don't re-open
#         if self.allele_info:
#             if not self.id:
#                 # need to save to have a flag collection
#                 self.save()
#
#             # record the fact that we did match
#             flag_collection = self.flag_collection_safe
#
#             self.save()  # we should be done changing values at this point
#             # though we may modify flags (which are external objects)
#
#             if not self.lab.external:
#                 # don't re-raise flags on external classifications,
#                 # they would have already been raised in the external system and resolved there
#
#                 # see if the match looks suspect
#                 c_hgvs = self.get(SpecialEKeys.C_HGVS)
#                 resolved_chgvs: Optional[CHGVS] = None
#                 transcript_comment: Optional[str] = None
#                 matching_warning_comment: Optional[str] = None
#                 compare_to: Optional[str] = None
#                 diff: Optional[CHGVSDiff] = None
#                 if self.get(SpecialEKeys.GENOME_BUILD):
#                     genome_build = self.get_genome_build()
#                     if genome_build == GenomeBuild.grch37():
#                         compare_to = self.chgvs_grch37
#                     elif genome_build == GenomeBuild.grch38():
#                         compare_to = self.chgvs_grch38
#
#                     if c_hgvs and compare_to:
#                         original_chgvs = CHGVS(c_hgvs)
#                         resolved_chgvs = CHGVS(compare_to)
#
#                         diff = original_chgvs.diff(resolved_chgvs)
#
#                         if diff:
#                             # DIFF_RAW_CGVS_EXPANDED is minor and expected process
#                             diff = diff & ~CHGVSDiff.DIFF_RAW_CGVS_EXPANDED
#
#                             if diff == CHGVSDiff.DIFF_TRANSCRIPT_VER:
#                                 transcript_comment = \
#                                     (f'For c.hgvs {c_hgvs} the transcripts are:\n\n'
#                                      f'{original_chgvs.transcript} (imported)\n{resolved_chgvs.transcript} (resolved)')
#                             elif diff:
#                                 important_diffs = chgvs_diff_description(diff)
#                                 diff_desc = ' - ' + '\n - '.join(important_diffs)
#                                 matching_warning_comment = (
#                                     f'Imported c.hgvs and matched c.hgvs differ in the following ways:\n\n{diff_desc}\n\n'
#                                     f'{original_chgvs.full_c_hgvs} (imported)\n{resolved_chgvs.full_c_hgvs} (resolved)')
#
#                 if transcript_comment:
#                     flag_collection.get_or_create_open_flag_of_type(
#                         flag_type=classification_flag_types.transcript_version_change_flag,
#                         comment=transcript_comment,
#                         only_if_new=True,
#                         reopen_if_bot_closed=True,
#                         data={'resolved': resolved_chgvs.transcript},
#                         close_other_data=True
#                     )
#                 else:
#                     flag_collection.close_open_flags_of_type(classification_flag_types.transcript_version_change_flag)
#
#                 if matching_warning_comment:
#                     flag_collection.get_or_create_open_flag_of_type(
#                         flag_type=classification_flag_types.matching_variant_warning_flag,
#                         comment=matching_warning_comment,
#                         only_if_new=True,
#                         reopen_if_bot_closed=True,
#                         data={'resolved': resolved_chgvs.full_c_hgvs},
#                         close_other_data=True
#                     )
#                 else:
#                     flag_collection.close_open_flags_of_type(classification_flag_types.matching_variant_warning_flag)
#
#         if allele := self.allele:
#             # now that everything's saved - see if 37 rep != 38 rep
#             # but note that the liftover might not be complete (if there is a liftover happening
#             # validation will be called again anyway)
#             allele.validate(liftover_complete=False)
