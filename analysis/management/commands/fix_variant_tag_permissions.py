from django.core.management.base import BaseCommand
from guardian.shortcuts import get_groups_with_perms, get_group_perms, assign_perm

from analysis.models.models_variant_tag import VariantTag
from analysis.models.models_analysis import Analysis

from library.guardian_utils import assign_permission_to_user_and_groups


class Command(BaseCommand):
    """ VariantTags used to always have an Analysis and thus used that permission
        we now have to handle tag existing w/o analyses and thus they have their own permissions """
    def handle(self, *args, **options):
        analysis_read_perm = Analysis.get_read_perm()
        analysis_write_perm = Analysis.get_write_perm()
        variant_tag_read_perm = VariantTag.get_read_perm()
        variant_tag_write_perm = VariantTag.get_write_perm()
        PERM_CONVERTER = {
            analysis_read_perm: variant_tag_read_perm,
            analysis_write_perm: variant_tag_write_perm,
        }

        num_tags = 0
        for variant_tag in VariantTag.objects.filter(analysis__isnull=False):
            # Convert existing analysis permissions
            for group in get_groups_with_perms(variant_tag.analysis):
                for old_perm in get_group_perms(group, variant_tag.analysis):
                    perm = PERM_CONVERTER[old_perm]
                    assign_perm(perm, group, variant_tag)

            # Default permissions for user
            assign_permission_to_user_and_groups(variant_tag.user, variant_tag)
            num_tags += 1

        print(f"Set permission on {num_tags} tags")
