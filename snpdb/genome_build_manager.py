from typing import Optional
from django.core.handlers.wsgi import WSGIRequest
from threadlocals.threadlocals import get_current_user, get_request_variable, set_request_variable, get_current_request
from snpdb.models import UserSettings, GenomeBuild
import re


# hardcoded to look for 37 or 38, needs support for other genome builds
GENOME_BUILD_RE = re.compile(r"(?:^|/)(GRCh3[78])(?:\?|$|/)")


class GenomeBuildManager:

    @staticmethod
    def get_current_genome_build():
        # See if we've set a request variable already with the genome build
        # after we calculate the genome build we cache it as a request variable
        genome_build: Optional[GenomeBuild] = None
        if genome_build := get_request_variable("building_manager_genome_build"):
            return genome_build

        # Now try to get the genome build, try fetching it from the following in priority order
        # 1: a URL GET parameter of genome_build, e.g. /someurl?genome_build=GRCh37
        # 2: as part of the URL path, e.g. /classifications/GRCh38
        # 3: if there's an authenticated user, from their UserSettings default_genome_build
        # 4: if all else fails, the first annotated genome build

        request: WSGIRequest
        if request := get_current_request():
            if gb := request.GET.get("genome_build"):
                try:
                    genome_build = GenomeBuild.get_name_or_alias(gb)
                except GenomeBuild.DoesNotExist:
                    pass
            elif genome_build_match := GENOME_BUILD_RE.search(request.get_full_path()):
                try:
                    genome_build = GenomeBuild.get_name_or_alias(genome_build_match.group(1))
                except GenomeBuild.DoesNotExist:
                    pass

        if genome_build is None:
            if user := get_current_user():
                genome_build = UserSettings.get_genome_build_or_default(user)

        if genome_build is None:
            genome_build = GenomeBuild.builds_with_annotation()[0]

        GenomeBuildManager.set_current_genome_build(genome_build)
        return genome_build

    @staticmethod
    def set_current_genome_build(genome_build: GenomeBuild):
        """
        Only needs to be called if the genome build wont be the user's default
        """
        set_request_variable("building_manager_genome_build", genome_build)
