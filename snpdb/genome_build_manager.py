from typing import Optional

from threadlocals.threadlocals import get_current_user, get_request_variable, set_request_variable

from snpdb.models import UserSettings, GenomeBuild


class GenomeBuildManager:

    @staticmethod
    def get_current_genome_build():
        # FIXME, allow overriding of current genome build in a consistent way
        # e.g. "genome_build" in the request
        genome_build: Optional[GenomeBuild] = None
        if genome_build := get_request_variable("building_manager_genome_build"):
            return genome_build

        if user := get_current_user():
            genome_build = UserSettings.get_genome_build_or_default(user)
        else:
            # TODO is there a better default than this?
            genome_build = GenomeBuild.builds_with_annotation()[0]

        GenomeBuildManager.set_current_genome_build(genome_build)
        return genome_build

    @staticmethod
    def set_current_genome_build(genome_build: GenomeBuild):
        """
        Only needs to be called if the genome build wont be the user's default
        """
        set_request_variable("building_manager_genome_build", genome_build)
