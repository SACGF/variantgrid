"""
Per-GFF assets from cdot's GitHub data releases.

cdot publishes two kinds of JSON.gz per data release:

    combo files, one per genome build, holding the latest copy of every transcript
        eg 'cdot-0.2.33.refseq.GRCh38.json.gz' - @see cdot.data_release.get_latest_combo_file_urls
    per-GFF files, one per gene set VEP can be built against
        eg 'cdot-0.2.33.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz'

The per-GFF files are what a GeneAnnotationRelease needs, as a release has to match the
exact gene/transcript versions VEP used.
"""
import gzip
import logging
import re
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional, Iterator

import requests
from cdot.data_release import get_latest_data_release

from library.constants import MINUTE_SECS

# The release segment can itself contain underscores (eg 'RS_2025_08'), so anchor on the
# annotation consortium rather than splitting on '_'
CDOT_RELEASE_ASSET_PATTERN = re.compile(
    r"^cdot-(?P<cdot_version>\d+\.\d+\.\d+)\.Homo_sapiens_"
    r"(?P<genome_build>.+?)_(?P<annotation_consortium>RefSeq|Ensembl)_"
    r"(?P<release>.+)\.(?P<file_type>gff|gtf)\.json\.gz$"
)


@dataclass(frozen=True)
class CdotReleaseAsset:
    cdot_version: str
    genome_build: str
    annotation_consortium: str
    release: str
    file_type: str
    url: str

    @property
    def filename(self) -> str:
        return self.url.rsplit("/", maxsplit=1)[-1]


def parse_release_asset_url(url: str) -> Optional[CdotReleaseAsset]:
    """ None for combo files and anything else that isn't a per-GFF release asset """
    filename = url.rsplit("/", maxsplit=1)[-1]
    if m := CDOT_RELEASE_ASSET_PATTERN.match(filename):
        return CdotReleaseAsset(url=url, **m.groupdict())
    return None


def get_latest_release_assets() -> list[CdotReleaseAsset]:
    """ Per-GFF assets from the most recent cdot data release whose schema matches our
        installed cdot library """
    assets = []
    for asset in get_latest_data_release().get("assets", []):
        if release_asset := parse_release_asset_url(asset["browser_download_url"]):
            assets.append(release_asset)
    return assets


def find_latest_release_asset(genome_build_name: str, annotation_consortium_label: str,
                              release: str) -> tuple[Optional[CdotReleaseAsset], list[CdotReleaseAsset]]:
    """ Returns (matching asset or None, all assets for that build/consortium) so callers can
        report what is available when there's no match """
    available = [a for a in get_latest_release_assets()
                 if a.genome_build.lower() == genome_build_name.lower()
                 and a.annotation_consortium.lower() == annotation_consortium_label.lower()]
    for asset in available:
        if asset.release.lower() == release.lower():
            return asset, available
    return None, available


@contextmanager
def download_cdot_json(url: str) -> Iterator[gzip.GzipFile]:
    """ Stream a cdot JSON.gz to a temp file on disk and yield an open handle on it.

        On disk rather than in memory because the importer reads the file twice (genes, then
        transcripts) so it has to be seekable, but it doesn't need to live in RAM - that would
        OOM on the larger builds. """
    logging.info("Downloading %s", url)
    with requests.get(url, stream=True, timeout=2 * MINUTE_SECS) as response:
        response.raise_for_status()
        with tempfile.NamedTemporaryFile(suffix=".json.gz") as tmp:
            for chunk in response.iter_content(chunk_size=1 << 20):
                tmp.write(chunk)
            tmp.flush()
            tmp.seek(0)
            with gzip.GzipFile(fileobj=tmp) as fz:
                yield fz
