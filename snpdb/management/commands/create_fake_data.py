"""
Fake data for development, so pages have something realistic to render.

One subcommand per kind of data, each owned by a module in the app the data belongs to. To add another,
write a class with HELP / add_arguments / create / delete (@see analysis.fake_variant_tags.FakeVariantTags)
and register it below.

Everything created is obviously fake, and every subcommand takes '--delete' to remove what it made.
"""
from django.core.management.base import BaseCommand

from analysis.fake_variant_tags import FakeVariantTags

FAKE_DATA = {
    "tags": FakeVariantTags,
}


class Command(BaseCommand):
    help = "Create obviously fake data for development"

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="subcommand", required=True)
        for name, fake_data_klass in FAKE_DATA.items():
            subparser = subparsers.add_parser(name, help=fake_data_klass.HELP)
            subparser.add_argument("--delete", action="store_true", help="Remove everything this created")
            fake_data_klass.add_arguments(subparser)

    def handle(self, *args, **options):
        fake_data = FAKE_DATA[options["subcommand"]](self.stdout)
        if options["delete"]:
            fake_data.delete(**options)
        else:
            fake_data.create(**options)
