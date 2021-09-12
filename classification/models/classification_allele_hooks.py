from django.dispatch.dispatcher import receiver

from classification.models.classification import Classification
from flags.models import FlagStatus, FlagCollection, Flag, FlagResolution
from genes.hgvs import HGVSMatcher
from snpdb.models import GenomeBuild
from snpdb.models.flag_types import allele_flag_types
from snpdb.models.models_variant import Allele, allele_validate_signal


@receiver(allele_validate_signal, sender=Allele)
def compare_chgvs(sender, allele: Allele, **kwargs):  # pylint: disable=unused-argument
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

        # There may be open flags against transcripts no longer used in classifications - remove them
        for flag in allele.flags_of_type(allele_flag_types.allele_37_not_38).filter(resolution__status=FlagStatus.OPEN):
            if transcript := flag.data.get("transcript"):
                if transcript not in classification_transcripts:
                    allele.close_open_flags_of_type(
                        allele_flag_types.allele_37_not_38,
                        data={'transcript': transcript},
                        comment="Transcript (version) no longer used",
                    )

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

            are_same = chgvs37 == chgvs38

            if are_same:
                allele.close_open_flags_of_type(
                    allele_flag_types.allele_37_not_38,
                    data={'transcript': transcript},
                    comment="37 and 38 representations are now the same"
                )
            else:
                # have to be careful as we only want one flag per transcript, but there might be multiple transcripts per allele
                # then for each transcript, we want to make sure the flag has the correct 37 and 38 representation (so can close and re-create if representation changes)
                raise_new = True
                if existing := allele.flag_collection_safe.get_flag_of_type(flag_type=allele_flag_types.allele_37_not_38, data={'transcript': transcript}):
                    if existing.data.get('chgvs37') == chgvs37 and existing.data.get('chgvs38') == chgvs38:
                        # record with the same transcript has the same 37 and 38 representation, no need to create a new flag
                        raise_new = False
                    elif existing.resolution.status == FlagStatus.OPEN:
                        # close old flag (with wrong 37 or 38 representation) so we can create a new one
                        existing.flag_action(comment="Resolution of 37, 38 changed", resolution=FlagResolution.objects.get(pk="closed"))

                if raise_new:
                    comment = (
                        f'Attached classification with transcript {transcript} appears as the following in \n\n'
                        f'{chgvs37} (GRCh37)\n'
                        f'{chgvs38} (GRCh38)')
                    allele.flag_collection_safe.get_or_create_open_flag_of_type(
                        flag_type=allele_flag_types.allele_37_not_38,
                        comment=comment,
                        data={'transcript': transcript, 'chgvs37': chgvs37, 'chgvs38': chgvs38},
                        only_if_new=True)
    else:
        # if there's no 37 or no 38, close any flag comparing the two
        allele.close_open_flags_of_type(allele_flag_types.allele_37_not_38, comment="Lacking representation in both 37 and 38")

