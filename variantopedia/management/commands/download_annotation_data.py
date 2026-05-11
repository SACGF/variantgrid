import os
from urllib.parse import urlparse

import requests
from django.core.management.base import BaseCommand
from tqdm import tqdm

from variantgrid.deployment_validation.annotation_files_check import annotation_data_exists


class Command(BaseCommand):
    help = "Check for missing annotation data files and download them from variantgrid.com"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="List what would be downloaded without downloading")
        parser.add_argument('--force', action='store_true',
                            help="Re-download even if file already exists")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        force = options["force"]

        checks = annotation_data_exists(flat=True, include_tbi_for_gz=True)

        missing = []
        for key, data in checks.items():
            if force or not data.get("valid"):
                url, dirname = _parse_fix(data.get("fix", ""))
                if url is None:
                    self.stderr.write(f"Could not parse fix for {key}: {data.get('fix')}")
                    continue
                filename = os.path.join(dirname, os.path.basename(urlparse(url).path))
                missing.append((key, url, dirname, filename))

        if not missing:
            self.stdout.write("All annotation data files present.")
            return

        self.stdout.write(f"Found {len(missing)} missing annotation file(s):")
        for key, url, _dirname, _filename in missing:
            self.stdout.write(f"  {key}: {url}")

        if dry_run:
            self.stdout.write("Dry run; not downloading.")
            return

        failures = []
        for key, url, dirname, filename in missing:
            self.stdout.write(f"\nDownloading {key}: {url}")
            try:
                os.makedirs(dirname, exist_ok=True)
                _download(url, filename)
                self.stdout.write(f"  -> {filename}")
            except Exception as e:
                self.stderr.write(f"  FAILED: {e}")
                failures.append((key, url, str(e)))

        if failures:
            self.stdout.write(f"\n{len(failures)} download(s) failed:")
            for key, url, err in failures:
                self.stdout.write(f"  {key}: {url}  ({err})")
        else:
            self.stdout.write("\nAll downloads complete.")


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


def _download(url, dest):
    tmp = dest + ".part"
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0)) or None
        with open(tmp, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=os.path.basename(dest)
        ) as pbar:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
    os.replace(tmp, dest)
