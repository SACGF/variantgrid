import subprocess
from tempfile import TemporaryDirectory

from django.test import SimpleTestCase

from library.git import Git


class GitVersionTest(SimpleTestCase):

    def test_version_counts_from_major_release_tag(self):
        with TemporaryDirectory() as repo_dir:
            def git(*args):
                subprocess.check_output(["git", "-c", "user.name=test", "-c", "user.email=test@example.com",
                                         "-c", "commit.gpgsign=false", "-c", "tag.gpgsign=false", *args],
                                        cwd=repo_dir)

            git("init", "--quiet")
            git("commit", "--quiet", "--allow-empty", "-m", "first")
            git("tag", "vg3.0")
            git("commit", "--quiet", "--allow-empty", "-m", "second")

            repo = Git(repo_dir)
            short_hash = repo.git_cmd("rev-parse", "--short", "HEAD")
            self.assertEqual(repo.version(3), f"vg3.0-1-g{short_hash}")
            self.assertEqual(repo.version(4), f"vg4-g{short_hash}")
