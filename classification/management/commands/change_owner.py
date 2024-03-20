from django.contrib.auth.models import User
from django.core.management import BaseCommand
from django.db.models import Q

from classification.models import Classification, ClassificationModification
from snpdb.models import Lab


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('lab_id', type='str', help="The lab_id with classifications that we want to change")
        parser.add_argument('old_user_id', type='str', help="The user_id we don't want anymore")
        parser.add_argument('new_user_id', type='str', help="The user_id we don't want anymore")
        parser.add_argument('test', action='store_true')

    @staticmethod
    def get_user(user_id:str) -> User:
        if user_id.isnumeric():
            return User.objects.get(pk=int(user_id))
        else:
            return User.objects.get(Q(username__iexact=user_id) | Q(email__iexact=user_id))

    @staticmethod
    def get_lab(lab_id: str) -> Lab:
        if lab_id.isnumeric():
            return Lab.objects.get(pk=int(lab_id))
        else:
            return Lab.objects.get(group_name=lab_id)

    def handle(self, *args, **options):
        lab = Command.get_lab(options["lab_id"])
        old_user = Command.get_user(options["old_user_id"])
        new_user = Command.get_user(options["new_user_id"])

        print(f"Lab = {lab}")
        print(f"Old User = {old_user}")
        print(f"New User = {new_user}")

        records_to_update_qs = ClassificationModification.objects.filter(
            classification__lab=lab,
            is_last_edited=True,
            is_last_published=True,
            published_evidence__owner__value=old_user.email,
            delta__owner__value=old_user.email
        )

        print(f"Records to update = {records_to_update_qs.count()}")

        if not options["test"]:
            for cm in records_to_update_qs:
                cm.published_evidence["owner"]["value"] = new_user.email
                cm.delta["owner"]["value"] = new_user.email
                cm.save(update_fields=["published_evidence", "delta"])
                cm.classification.evidence["owner"]["value"] = new_user.email
                cm.save(update_fields=["evidence"])
            print("Done")
