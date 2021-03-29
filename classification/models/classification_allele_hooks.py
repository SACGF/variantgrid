from django.dispatch.dispatcher import receiver

from genes.hgvs import HGVSMatcher
from snpdb.models import GenomeBuild
from snpdb.models.flag_types import allele_flag_types
from snpdb.models.models_variant import Allele, allele_validate_signal
from classification.models.classification import Classification


@receiver(allele_validate_signal, sender=Allele)
def compare_chgvs(sender, allele: Allele, **kwargs):
    vcs = Classification.objects.filter(variant__in=allele.variants).order_by('id')
    v37 = allele.grch37
    v38 = allele.grch38
    if v37 and v38:
        matcher37 = HGVSMatcher(genome_build=GenomeBuild.grch37())
        matcher38 = HGVSMatcher(genome_build=GenomeBuild.grch38())
        transcripts = set()
        existing_flags = allele.flags_of_type(allele_flag_types.allele_37_not_38).values_list('data__transcript',
                                                                                              flat=True)
        for transcript in existing_flags:
            if transcript:
                transcripts.add(transcript)
        vc: Classification
        for vc in vcs:
            transcript = vc.transcript
            if transcript:
                transcripts.add(transcript)

        for transcript in transcripts:
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
                    data={'transcript': transcript}
                )
            else:
                comment = (
                    f'Attached classification with transcript {transcript} appears as the following in \n\n'
                    f'{chgvs37} (GRCh37)\n'
                    f'{chgvs38} (GRCh38)')
                allele.flag_collection_safe.get_or_create_open_flag_of_type(
                    flag_type=allele_flag_types.allele_37_not_38,
                    comment=comment,
                    data={'transcript': transcript},
                    only_if_new=True)
    else:
        # if there's no 37 or no 38, close any flag comparing the two
        allele.close_open_flags_of_type(allele_flag_types.allele_37_not_38)
