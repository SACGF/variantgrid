from django.core.management import BaseCommand


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--file', type=str, required=True)

    def handle(self, *args, **options):
        """
        Take a file of internal classification IDs, and rematch them
        """
        raise Exception("classification_re_matching has been removed (rematching is now done against allele infos")
