import logging

from django.core.management import BaseCommand
from django.db.models import QuerySet

from classification.models import ImportedAlleleInfo
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        At the moment this is only used for liftover but I'm using subcommands so we can have multiple fixes
    """
    def add_arguments(self, parser):
        subparsers = parser.add_subparsers(dest="command", help="Subcommands for different fixes")
        _liftover_parser = subparsers.add_parser('liftover', help='Re-link to any missing variants in other builds')

    def handle(self, *args, **options):
        sub_commands = {
            "liftover": self._handle_liftover,
        }

        if cmd := options["command"]:
            commands = [cmd]
        else:
            # None given run all
            commands = sub_commands.keys()

        for cmd in commands:
            logging.info("Running %s", cmd)
            sub_commands[cmd](*args, **options)

    @staticmethod
    def get_imported_allele_info_liftover_builds_and_qs() -> list[tuple[GenomeBuild, QuerySet[ImportedAlleleInfo]]]:
        builds_and_liftover = []
        for genome_build in ImportedAlleleInfo.supported_genome_builds:
            field = ImportedAlleleInfo._genome_build_to_attr(genome_build)
            kwargs = {
                f"{field}__isnull": True,
                "allele__variantallele__genome_build": genome_build,
            }
            iai_qs = ImportedAlleleInfo.objects.filter(**kwargs)
            builds_and_liftover.append((genome_build, iai_qs))
        return builds_and_liftover

    def _handle_liftover(self, *args, **options):
        """ Historically some ImportedAlleleInfos were lifted over didn't have resolved variants linked
            See https://github.com/SACGF/variantgrid_private/issues/3761#issuecomment-2777293132
        """
        for genome_build, iai_qs in self.get_imported_allele_info_liftover_builds_and_qs():
            if num_missing := iai_qs.count():
                logging.info("ImportedAlleleInfo: %d missing resolved variants for %s", num_missing, genome_build)
                for allele_info in iai_qs:
                    logging.info("Refreshing %s", allele_info)
                    allele_info.refresh_and_save(force_update=True)

                still_missing = iai_qs.count()
                fixed = num_missing - still_missing
                logging.info("ImportedAlleleInfo - refresh matched %d, %d still missing resolved variants for %s",
                             fixed, still_missing, genome_build)

