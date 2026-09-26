"""
Obviously fake data for a development database: 'create_fake_data all' on an empty database gives users, labs,
genes, annotated variants, a trio with genotypes, an analysis, tags and classifications.

    manage.py create_fake_data list                                  # the steps, what each makes, what it needs
    manage.py create_fake_data all [--genome-build GRCh37] [--seed N] [--delete]
    manage.py create_fake_data <step> [--delete] [step options]      # runs the steps it requires first

Each app owns the steps for its own models, in '<app>/fake_data.py' or a 'fake_data' package. To add one, subclass
library.fake_data.FakeData there and decorate it with @register - this command finds it without importing the app.
'--delete' removes only what the named step made; 'all --delete' removes every step's data, last step first.
"""
import argparse

from django.core.management.base import BaseCommand
from django.db import transaction

from library.fake_data import FakeDataContext, discover, in_dependency_order
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    category = "dev"
    help = "Create obviously fake data for development"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fake_data = discover()

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="step", required=True)
        subparsers.add_parser("list", help="The steps, what each one makes and the steps it requires")

        # 'all' also takes every step's options, and one given applies to each step that has it. Two steps may
        # share an option name, hence "resolve", and defaults are suppressed so only what was given is passed on
        all_parser = subparsers.add_parser("all", help="Every step", conflict_handler="resolve")
        self._add_shared_arguments(all_parser)
        for fake_data_klass in self.fake_data.values():
            step_parser = subparsers.add_parser(fake_data_klass.name, help=fake_data_klass.help)
            self._add_shared_arguments(step_parser)
            fake_data_klass.add_arguments(step_parser)

            existing_actions = set(all_parser._actions)
            fake_data_klass.add_arguments(all_parser)
            for action in all_parser._actions:
                if action not in existing_actions:
                    action.default = argparse.SUPPRESS

    @staticmethod
    def _add_shared_arguments(parser):
        parser.add_argument("--genome-build", default="GRCh37", help="Build to make the data in")
        parser.add_argument("--seed", type=int, default=1816, help="Same seed, same data")
        parser.add_argument("--delete", action="store_true",
                            help="Remove what this step made (not what it requires) - for 'all', everything")

    def handle(self, *args, **options):
        step = options["step"]
        if step == "list":
            self._list()
            return

        context = FakeDataContext(genome_build=GenomeBuild.get_name_or_alias(options["genome_build"]),
                                  seed=options["seed"], stdout=self.stdout)
        if step == "all":
            step_names = in_dependency_order(self.fake_data, self.fake_data)
            given_options = options
        else:
            step_names = in_dependency_order([step], self.fake_data)
            given_options = {}  # the named step's options are its own - requirements run on their defaults

        def step_options(name: str) -> dict:
            return options if name == step else self._step_options(name, given_options)

        if options["delete"]:
            for name in reversed(step_names if step == "all" else [step]):
                self.stdout.write(self.style.MIGRATE_HEADING(f"Deleting fake {name}"))
                with transaction.atomic():
                    self.fake_data[name]().delete(context, **step_options(name))
            return

        for name in step_names:
            self.stdout.write(self.style.MIGRATE_HEADING(f"Fake {name}"))
            with transaction.atomic():
                self.fake_data[name]().create(context, **step_options(name))

    def _step_options(self, name: str, given_options: dict) -> dict:
        """ The step's defaults, with whichever of its options were given """
        parser = argparse.ArgumentParser()
        self.fake_data[name].add_arguments(parser)
        step_options = vars(parser.parse_args([]))
        step_options.update({key: value for key, value in given_options.items() if key in step_options})
        return step_options

    def _list(self):
        for name in in_dependency_order(self.fake_data, self.fake_data):
            fake_data_klass = self.fake_data[name]
            requires = f" (requires {', '.join(fake_data_klass.requires)})" if fake_data_klass.requires else ""
            self.stdout.write(f"{name}{requires}: {fake_data_klass.help}")
