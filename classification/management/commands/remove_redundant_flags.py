from django.core.management import BaseCommand

from flags.models import FlagType, Flag


class Command(BaseCommand):

    def handle(self, *args, **options):

        def delete_flag_of_type(type_str: str):
            if type_obj := FlagType.objects.filter(pk=type_str).first():
                type_obj_count = Flag.objects.filter(flag_type=type_obj).count()
                print(f"{type_str} has {type_obj_count} flags... deleting")
                type_obj.delete()
                # cascade delete will delete all the flags of that type
            else:
                print(f"{type_str} flag already deleted")

        delete_flag_of_type('classification_matching_variant')
        delete_flag_of_type('classification_matching_variant_warning')
        delete_flag_of_type('classification_transcript_version_change')
        delete_flag_of_type('allele_37_not_38')
        delete_flag_of_type('allele_missing_38')
        delete_flag_of_type('allele_missing_37')
