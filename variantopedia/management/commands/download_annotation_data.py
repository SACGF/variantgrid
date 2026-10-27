import os
import shlex
import tarfile
from typing import NamedTuple
from urllib.parse import urlparse

import requests
from django.core.management.base import BaseCommand
from tqdm import tqdm

from variantgrid.deployment_validation.annotation_files_check import annotation_data_exists

# vep_config entries that point at a directory of files rather than a single file.
# These are hosted as '<dir>.tar.gz' on variantgrid.com and extracted in place.
_DIRECTORY_ARCHIVE_KEYS = {"maxentscan"}


# Entries that are programs, made executable once downloaded
_EXECUTABLE_KEYS = {"somalier"}


def _base_key(key):
    # key may be plain ("maxentscan") or flattened ("GRCh37/maxentscan")
    return key.rsplit("/", 1)[-1]


def _is_directory_archive(key):
    return _base_key(key) in _DIRECTORY_ARCHIVE_KEYS


def _is_executable(key):
    return _base_key(key) in _EXECUTABLE_KEYS


class _Download(NamedTuple):
    key: str
    url: str
    dirname: str
    filename: str
    is_archive: bool
    optional: bool
    refresh: bool


class Command(BaseCommand):
    category = "ops"
    help = "Check for missing annotation data files and download them from variantgrid.com"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="List what would be downloaded without downloading")
        parser.add_argument('--force', action='store_true',
                            help="Re-download even if file already exists")
        parser.add_argument('--bash-script', nargs='?', const='-', metavar='PATH',
                            help="Write a bash script that downloads the missing and outdated annotation files "
                                 "(to PATH, or stdout if no path given) instead of downloading now")
        parser.add_argument('--download-all', action='store_true',
                            help="With --bash-script, write every annotation file into the script, not just the "
                                 "missing and outdated ones (same as --force)")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        force = options["force"] or options["download_all"]
        bash_script = options["bash_script"]

        checks = annotation_data_exists(flat=True, include_tbi_for_gz=True, include_csi=True)
        downloads = self._select_downloads(checks, force)

        if bash_script is not None:
            self._write_bash_script(downloads, bash_script)
            return

        missing = []
        not_hosted = []
        for download in downloads:
            if download.optional and not _is_hosted(download.url):
                not_hosted.append(download)
            else:
                missing.append(download)

        if not_hosted:
            self.stdout.write(f"{len(not_hosted)} optional index file(s) are not hosted on variantgrid.com, skipped "
                              f"(build them with annotation/annotation_data/generate_annotation/build_tabix_csi.sh):")
            for download in not_hosted:
                self.stdout.write(f"  {download.key}: {download.url}")

        if not missing:
            self.stdout.write("All annotation data files present.")
            return

        self.stdout.write(f"Found {len(missing)} missing or outdated annotation file(s):")
        for download in missing:
            self.stdout.write(f"  {download.key}: {download.url}")

        if dry_run:
            self.stdout.write("Dry run; not downloading.")
            return

        failures = []
        for download in missing:
            self.stdout.write(f"\nDownloading {download.key}: {download.url}")
            try:
                os.makedirs(download.dirname, exist_ok=True)
                _download(download.url, download.filename)
                if download.is_archive:
                    _extract_tarball(download.filename, download.dirname)
                    os.remove(download.filename)
                elif _is_executable(download.key):
                    os.chmod(download.filename, 0o755)
                self.stdout.write(f"  -> {download.filename}")
            except Exception as e:
                self.stderr.write(f"  FAILED: {e}")
                failures.append((download.key, download.url, str(e)))

        if failures:
            self.stdout.write(f"\n{len(failures)} download(s) failed:")
            for key, url, err in failures:
                self.stdout.write(f"  {key}: {url}  ({err})")
        else:
            self.stdout.write("\nAll downloads complete.")

    def _select_downloads(self, checks, force) -> list[_Download]:
        """The files to fetch: missing, stale (a plugin we ship whose hosted copy differs) or all of them with force.
        A file whose fix is not a wget ships with the repo and is never downloaded."""
        downloads = []
        for key, data in checks.items():
            fix = data.get("fix", "")
            if "wget" not in fix:
                continue
            url, dirname = _parse_fix(fix)
            if url is None:
                self.stderr.write(f"Could not parse fix for {key}: {fix}")
                continue
            filename = os.path.join(dirname, os.path.basename(urlparse(url).path))
            stale = data.get("refresh") and data.get("valid") and _differs_from_hosted(url, filename)
            if not (force or stale or not data.get("valid")):
                continue
            is_archive = _is_directory_archive(key)
            if is_archive:
                # The config value is a directory; download the hosted tarball and extract it
                url = url + ".tar.gz"
                filename = os.path.join(dirname, os.path.basename(urlparse(url).path))
            downloads.append(_Download(key, url, dirname, filename, is_archive,
                                       optional=bool(data.get("optional")), refresh=bool(data.get("refresh"))))
        return downloads

    def _write_bash_script(self, downloads, path):
        lines = [
            "#!/bin/bash",
            "# Generated by 'manage.py download_annotation_data --bash-script'",
            "# Downloads VariantGrid annotation data files.",
            "set -euo pipefail",
            "",
        ]
        if not downloads:
            lines.append("# All annotation data files present.")
        for download in downloads:
            quoted_dir = shlex.quote(download.dirname)
            quoted_url = shlex.quote(download.url)
            quoted_filename = shlex.quote(download.filename)
            lines.append(f"# {download.key}")
            lines.append(f"mkdir -p {quoted_dir}")
            if download.is_archive:
                lines.append(f"wget -c -P {quoted_dir} {quoted_url}")
                lines.append(f"tar -xzf {quoted_filename} -C {quoted_dir}")
                lines.append(f"rm {quoted_filename}")
            elif download.optional:
                lines.append(f"wget -c -P {quoted_dir} {quoted_url} || echo 'optional, not hosted: {download.key}'")
            elif download.refresh:
                # ours; the hosted copy replaces whatever the plugins checkout holds
                lines.append(f"wget -O {quoted_filename} {quoted_url}")
            else:
                lines.append(f"wget -c -P {quoted_dir} {quoted_url}")
                if _is_executable(download.key):
                    lines.append(f"chmod a+x {quoted_filename}")
            lines.append("")

        script = "\n".join(lines)
        if path == "-":
            self.stdout.write(script)
        else:
            with open(path, "w", encoding="utf-8") as f:
                f.write(script + "\n")
            os.chmod(path, 0o755)
            self.stdout.write(f"Wrote download script to {path}")


def _parse_fix(fix):
    """Parse a fix instruction of the form 'cd <dir>;wget <url>'."""
    if not fix or ";" not in fix:
        return None, None
    cd_part, _, wget_part = fix.partition(";")
    cd_part = cd_part.strip()
    wget_part = wget_part.strip()
    if not cd_part.startswith("cd ") or not wget_part.startswith("wget "):
        return None, None
    return wget_part[len("wget "):].strip(), cd_part[len("cd "):].strip()


def _differs_from_hosted(url, filename) -> bool:
    """For a "refresh" entry (a plugin we ship), whether the installed file differs from the hosted one, which is
    authoritative. Unreachable counts as unchanged so an offline run does not re-download."""
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException:
        return False
    with open(filename, "rb") as f:
        return f.read() != resp.content


def _is_hosted(url) -> bool:
    """Whether variantgrid.com serves this URL; optional files (CSI indexes) are only fetched when it does."""
    try:
        resp = requests.head(url, allow_redirects=True, timeout=30)
    except requests.RequestException:
        return False
    return resp.ok


def _extract_tarball(tarball, dest):
    with tarfile.open(tarball, "r:gz") as tf:
        try:
            tf.extractall(dest, filter="data")  # guards against path traversal (Python 3.12+)
        except TypeError:
            tf.extractall(dest)


def _download(url, dest):
    tmp = dest + ".part"
    resume_from = os.path.getsize(tmp) if os.path.exists(tmp) else 0

    headers = {"Range": f"bytes={resume_from}-"} if resume_from else {}
    with requests.get(url, stream=True, timeout=60, headers=headers) as resp:
        if resume_from and resp.status_code == 416:
            # Range not satisfiable - the partial file is already the complete download.
            os.replace(tmp, dest)
            return
        resp.raise_for_status()

        content_length = int(resp.headers.get("content-length", 0)) or None
        if resume_from and resp.status_code == 206:
            # Server honoured the range request; append to the existing partial.
            mode = "ab"
            total = content_length + resume_from if content_length else None
            initial = resume_from
        else:
            # No partial, or the server ignored the range (200 OK) - start fresh.
            mode = "wb"
            total = content_length
            initial = 0

        with open(tmp, mode) as f, tqdm(
            total=total, initial=initial, unit="B", unit_scale=True, desc=os.path.basename(dest)
        ) as pbar:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
    os.replace(tmp, dest)
