"""
Git: the repo's current hash, version, branch, last modified date and a GitHub link, read by the site footer,
`vg status`, deployment checks and the VARIANTGRID_VERSION setting.
"""
import re
import subprocess
from functools import cached_property

from dateutil.parser import parse

from library.utils import is_url


class Git:
    ALLOWED_SUBCOMMANDS = {"rev-parse", "log", "config", "describe"}

    def __init__(self, directory=None):
        self.directory = directory

    def git_cmd(self, subcommand: str, *params) -> str:
        if subcommand not in self.ALLOWED_SUBCOMMANDS:
            raise ValueError(f"git subcommand '{subcommand}' is not in the allow-list {self.ALLOWED_SUBCOMMANDS}")
        output = subprocess.check_output(["git", subcommand, *params], cwd=self.directory, stderr=subprocess.PIPE)
        return output.decode().strip()

    @cached_property
    def hash(self) -> str:
        return self.git_cmd("rev-parse", "HEAD")

    def version(self, major_version: int) -> str:
        """ Counted from the latest vg<major>.* release tag, e.g. 'vg4.0-12-gc174556',
            or 'vg4-gc174556' while that major has no release tag yet """
        try:
            return self.git_cmd("describe", "--tags", "--match", f"vg{major_version}.*")
        except subprocess.CalledProcessError:
            short_hash = self.git_cmd("rev-parse", "--short", "HEAD")
            return f"vg{major_version}-g{short_hash}"

    @cached_property
    def last_modified_date(self):
        date_string = self.git_cmd("log", "-1", "--format=%cd")
        return parse(date_string)

    @cached_property
    def branch(self) -> str:
        return self.git_cmd("rev-parse", "--abbrev-ref", "HEAD")

    @cached_property
    def site(self):
        git_site = self.git_cmd("config", "--get", "remote.origin.url")
        return re.sub(r"([^/]+@|\.git)", "", git_site)

    @cached_property
    def branch_link(self):
        git_branch_link = None
        if is_url(self.site) and self.branch:
            git_branch_link = f"{self.site}/commits/{self.branch}"
        return git_branch_link
