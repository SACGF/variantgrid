from django.contrib.auth.models import User, Group
from django.core.management import BaseCommand

from classification.models import Classification
from snpdb.models import Lab, Organization


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--delete', action='store_true', default=False)

    def handle(self, *args, **options):
        org = Organization.objects.get(group_name="sa_pathology")
        legacy, _ = Lab.objects.get_or_create(group_name="sa_pathology/legacy", defaults={
            "name": "LEGACY",
            "organization": org,
            "city": "Adelaide",
            "state": "SA",
            "country": "Australia",
            "lat": -34.91203,
            "long": 138.600349,
            "url": "http://www.sapathology.sa.gov.au",
            "css_class": "sa-pathology",
            "classification_config": {"zygosity": {"mandatory": False}, "condition": {"mandatory": False}, "clinical_significance": {"mandatory": False}}
        })

        old_legacy_labs = Lab.objects.filter(group_name__startswith="sa_pathology/legacy_")
        old_classifications = Classification.objects.filter(lab__in=old_legacy_labs)
        print(f"Classifications to update = {old_classifications.count()}")
        legacy_class: Classification
        count = 0
        for legacy_class in old_classifications:
            if count % 100 == 0:
                print(f"Processed {count} classifications")
            legacy_class.lab = legacy
            legacy_class.save(update_fields=["lab"])
            legacy_class.fix_permissions()
            count += 1
        print(f"Processed {count} classifications")

        old_group_names = []
        for lab in old_legacy_labs:
            old_group_names.append(lab.group_name)

        for user in User.objects.all():
            if user.groups.filter(name__in=old_group_names):
                for group_name in old_group_names:
                    if user.groups.filter(name=group_name):
                        user.groups.remove(Group.objects.get(name=group_name))
                if not user.groups.filter(name='sa_pathology/legacy'):
                    user.groups.add(Group.objects.get(name="sa_pathology/legacy"))
                print(f"Updated user {user.username} to belong to just one legacy")

        if options["delete"]:
            print("Deleting labs")
            print(old_group_names)
            for lab in old_legacy_labs:
                Group.objects.filter(name=lab.group_name).delete()
                lab.delete()
