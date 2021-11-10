from dataclasses import dataclass
from typing import Set, Dict, Any

from django.dispatch.dispatcher import receiver

from classification.models.classification import Classification
from flags.models import FlagStatus, FlagComment, Flag
from genes.hgvs import HGVSMatcher
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild
from snpdb.models.flag_types import allele_flag_types
from snpdb.models.models_variant import Allele, allele_validate_signal


@dataclass(frozen=True)
class TranscriptDifference:
    transcript: str
    chgvs37: str
    chgvs38: str

    @property
    def comment(self) -> str:
        return (
            f'Attached classification with transcript {self.transcript} appears as the following in \n\n'
            f'{self.chgvs37} (GRCh37)\n'
            f'{self.chgvs38} (GRCh38)')

    @property
    def data(self) -> Dict[str, Any]:
        return {
            "transcript": self.transcript,
            'chgvs37': self.chgvs37,
            'chgvs38': self.chgvs38
        }

    @staticmethod
    def from_flag(flag: Flag) -> 'TranscriptDifference':
        data = flag.data or {}
        return TranscriptDifference(
            transcript=data.get('transcript'),
            chgvs37=data.get('chgvs37'),
            chgvs38=data.get('chgvs38')
        )


@receiver(allele_validate_signal, sender=Allele)
def compare_chgvs(sender, allele: Allele, **kwargs):  # pylint: disable=unused-argument
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
