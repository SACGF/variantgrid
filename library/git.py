import re
import subprocess

from dateutil.parser import parse

from lazy import lazy

from library.utils import is_url


class Git:
    def __init__(self, directory=None):
        self.directory = directory

    def git_cmd(self, params):
        output = subprocess.check_output(["git"] + params, cwd=self.directory)
        return output.decode().strip()

    @lazy
    def hash(self):
        return self.git_cmd(["rev-parse", "HEAD"])

    @lazy
    def last_modified_date(self):
        date_string = self.git_cmd(["log", "-1", "--format=%cd"])
        return parse(date_string)

    @lazy
    def branch(self):
        return self.git_cmd(["rev-parse", "--abbrev-ref", "HEAD"])

    @lazy
    def site(self):
        git_site = self.git_cmd(["config", "--get", "remote.origin.url"])
        return re.sub(r"([^/]+@|\.git)", "", git_site)

    @lazy
    def branch_link(self):
        git_branch_link = None
        if is_url(self.site) and self.branch:
            git_branch_link = f"{self.site}/commits/{self.branch}"
        return git_branch_link

