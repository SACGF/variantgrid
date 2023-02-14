from django.core.management import BaseCommand


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--file', type=str, required=True)

    def handle(self, *args, **options):
        """
        Take a file of internal classification IDs, and rematch them
        """
        raise Exception("classification_re_matching has been removed (rematching is now done against allele infos")
        # filename = options["file"]
        # qs: QuerySet[Classification] = Classification.objects.none()
        #
        # with open(filename) as file:
        #     all_ids = set(int(cid) for cid in re.compile(r"\d+").findall(file.read()))
        #     print(f"Found {len(all_ids)} IDs in file")
        #     qs = Classification.objects.filter(pk__in=all_ids)
        #     if qs.count() != len(all_ids):
        #         for c in qs:
        #             all_ids.remove(c.pk)
        #
        #         print(f"Number of classification IDs not present in the database = {len(all_ids)}")
        #         print(f"Example missing IDs = {list(all_ids)[0:5]}")
        #         return
        #
        # if qs:
        #     user = admin_bot()
        #     reattempt_variant_matching(user, qs)
        #     print("Re-matching has been queued")
