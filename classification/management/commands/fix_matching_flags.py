from typing import Optional

from django.core.management import BaseCommand

from classification.models import classification_flag_types, Classification
from flags.models import Flag, FlagComment
import re

from genes.hgvs import CHGVS
from snpdb.models import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--apply', action='store_true', default=False)

    def handle(self, *args, **options):

        apply = options.get('apply')
        if apply:
            print("Will be applying changes")

        flags_of_interest = Flag.objects.filter(flag_type__in=[
            classification_flag_types.transcript_version_change_flag,
            classification_flag_types.matching_variant_warning_flag
        ]).filter(data__isnull=True)

        re_find_in_comment = [
            re.compile(r'^(.*) \(resolved\)$', re.MULTILINE),
            re.compile(r'^(.*) \(matched\)$', re.MULTILINE)
        ]

        for flag in flags_of_interest:

            resolved: Optional[str] = None

            opening_comment: FlagComment
            if opening_comment := flag.flagcomment_set.order_by('created').first():
                if comment_text := opening_comment.text:
                    for pattern in re_find_in_comment:
                        if m := pattern.search(comment_text):
                            resolved = m[1]
                            break
                    if not resolved:
                        print(f'** Couldnt find the resolved value in flag {flag.id}')
                        print('---')
                        print(comment_text)
                        print('---')

            #  find the classification, in the case of not being able to work it out, let's just leave things None
            #  and flag will be closed and re-opened since the data doesn't match
            """
            if not resolved:
                c: Classification
                if c := Classification.objects.filter(flag_collection_id=flag.collection_id).first():
                    if genome_build := c.get_genome_build():
                        compare_to: Optional[str] = None
                        if genome_build == GenomeBuild.grch37():
                            compare_to = c.chgvs_grch37
                        elif genome_build == GenomeBuild.grch38():
                            compare_to = c.chgvs_grch38

                        if compare_to:
                            compare_to_chgvs = CHGVS(compare_to)
                            if flag.flag_type == classification_flag_types.matching_variant_warning_flag:
                                resolved = compare_to_chgvs.full_c_hgvs
                            elif flag.flag_type == classification_flag_types.transcript_version_change_flag:
                                resolved = compare_to_chgvs.transcript
            """
            if resolved:
                print(f"Resolved flag {flag.id} to {resolved}")
                if apply:
                    flag.data = {'resolved': resolved}
                    flag.save()

        if apply:
            c: Classification
            for index, c in enumerate(Classification.objects.all()):
                if index % 100 == 0:
                    print(f"Processed {index} classifications for updating flags")
                c.update_cached_c_hgvs()
