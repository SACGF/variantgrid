import re
import subprocess
from functools import cached_property

from dateutil.parser import parse

from library.utils import is_url


class Git:
    def __init__(self, directory=None):
        self.directory = directory

    def git_cmd(self, params) -> str:
        output = subprocess.check_output(["git"] + params, cwd=self.directory)
        return output.decode().strip()

    @cached_property
    def hash(self) -> str:
        return self.git_cmd(["rev-parse", "HEAD"])

    @cached_property
    def last_modified_date(self):
        date_string = self.git_cmd(["log", "-1", "--format=%cd"])
        return parse(date_string)

    @cached_property
    def branch(self) -> str:
        return self.git_cmd(["rev-parse", "--abbrev-ref", "HEAD"])

    @cached_property
    def site(self):
        git_site = self.git_cmd(["config", "--get", "remote.origin.url"])
        return re.sub(r"([^/]+@|\.git)", "", git_site)

    @cached_property
    def branch_link(self):
        git_branch_link = None
        if is_url(self.site) and self.branch:
            git_branch_link = f"{self.site}/commits/{self.branch}"
        return git_branch_link
