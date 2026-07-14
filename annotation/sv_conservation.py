"""
Structural-variant conservation scoring with pyBigWig (#1657).

VEP's `--custom ...,summary_stats=max,format=bigwig` computes the max bigWig value over a variant's
reference span. On 1 bp small variants that's cheap, but on multi-Mb SVs it forces VEP's Perl bigWig
reader to walk conservation values across the entire span (O(span)) for each of the 4 conservation
tracks, which makes large-SV annotation never finish.

This module reproduces the same 4 conservation `_max` columns (phastCons/phyloP) with pyBigWig, using
an aligned-window decomposition that is bit-exact vs VEP but ~ms/SV even on multi-Mb spans. The 4
conservation `--custom` args are dropped from the SV VEP command (their ColumnVEPField rows are
STANDARD-only, see annotation.vep_columns), and this stage populates the same DB columns on the import
path (see annotation.tasks.annotate_variants + BulkVEPVCFAnnotationInserter).
"""
import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Optional

import pyBigWig
from django.conf import settings

from annotation import vep_columns
from annotation.models.models_enums import VEPCustom
from annotation.vep_config import VEPConfig
from library.utils.file_utils import open_handle_gzip
from snpdb.models.models_genome import GenomeBuild

# The 4-per-build conservation bigWig tracks whose _max columns VEP computes via --custom.
# 100way applies to both builds; 30way is GRCh38-only, 46way is GRCh37-only (see vep_columns).
CONSERVATION_VEP_CUSTOMS = (
    VEPCustom.PHASTCONS_100_WAY,
    VEPCustom.PHASTCONS_30_WAY,
    VEPCustom.PHASTCONS_46_WAY,
    VEPCustom.PHYLOP_100_WAY,
    VEPCustom.PHYLOP_30_WAY,
    VEPCustom.PHYLOP_46_WAY,
)

# Aligned-window decomposition constants (validated in claude/sv_pybigwig_plan.md / PoC).
# Interior is read from zoom summaries on a fine grid of bin width G (so its zoom reduction is <= G),
# with exact per-base scans on the two ragged ends of width E. The invariant E > G guarantees any
# interior zoom bin that straddles the query boundary reaches at most ~G past the interior - which
# lands inside an exact-scanned end and is capped correctly. Net: bit-exact MAX, but fast.
_G = 16384       # interior zoom bin width
_E = 2 * _G      # exact-end width (must exceed _G so straddles are capped)


@dataclass(frozen=True)
class ConservationTrack:
    """ One conservation bigWig track + the VariantAnnotation column its _max value populates. """
    name: str          # VEP custom short_name, e.g. 'phastCons100way'
    path: str          # bigWig (.bw) file path
    db_column: str     # VariantGrid column, e.g. 'phastcons_100_way_vertebrate'


def span_max_aligned(bw, chrom: str, s: int, e: int) -> Optional[float]:
    """ MAX bigWig value over 0-based half-open [s, e) - bit-exact, fast on huge spans.

        Exact per-base on the two ragged ends (width _E); interior read from zoom summaries on a fine
        aligned grid (bin width _G) in a single stats() call. Invariant _E > _G caps any zoom straddle
        inside the exact-scanned ends. Spans <= 2*_E just do a full exact scan. """
    if e - s <= 2 * _E:
        return bw.stats(chrom, s, e, type="max", exact=True)[0]
    lo, hi = s + _E, e - _E
    vals = [
        bw.stats(chrom, s, lo, type="max", exact=True)[0],
        bw.stats(chrom, hi, e, type="max", exact=True)[0],
    ]
    nb = max(1, (hi - lo) // _G)
    vals.extend(v for v in bw.stats(chrom, lo, hi, type="max", exact=False, nBins=nb) if v is not None)
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


def sv_conservation_window(pos: int, end: Optional[int], svlen: Optional[int]) -> tuple[int, int]:
    """ VEP's summary_stats reference footprint for an SV, as 0-based half-open [s, e).

        VEP takes summary over 1-based inclusive [POS, max(END, POS+|SVLEN|)]:
        - DEL/DUP/INV: END (== POS+|SVLEN|)
        - INS: END is only POS+1, so the footprint is POS+|SVLEN| """
    s = pos - 1                                   # 0-based inclusive start
    e = max(end or pos, pos + abs(svlen or 0))    # 0-based exclusive end
    return s, e


def get_sv_conservation_tracks(genome_build: GenomeBuild) -> list[ConservationTrack]:
    """ Conservation bigWig tracks configured for this build, paired with the DB column each populates.

        Derived from the ColumnVEPField registry (the summary_stats='max' conservation rows) so the
        track->column mapping stays in one place. Drops any track whose bigWig data file isn't
        configured for the build. """
    vc = VEPConfig(genome_build)
    tracks: list[ConservationTrack] = []
    for vep_custom in CONSERVATION_VEP_CUSTOMS:
        # pipeline_type is intentionally left unfiltered - the _max rows are STANDARD-only so that
        # get_vep_command drops them on the SV command, but we still want their track->column mapping.
        cvf_list = [c for c in vep_columns.filter_for(vep_config=vc, vep_custom=vep_custom)
                    if c.summary_stats == 'max']
        if not cvf_list:
            continue  # not configured / no data file / wrong build or VEP version
        try:
            path = vc[vep_custom.label.lower()]
        except KeyError:
            continue
        db_column = cvf_list[0].variant_grid_columns[0]
        tracks.append(ConservationTrack(name=vep_custom.label, path=path, db_column=db_column))
    return tracks


def read_sv_vcf_variants(vcf_filename: str) -> list[tuple[int, str, int, Optional[int], Optional[int]]]:
    """ Read (variant_id, chrom, pos, end, svlen) for each record in a dumped/annotated SV VCF.
        variant_id / END / SVLEN come from the INFO column (written by snpdb.variants_to_vcf). """
    variants = []
    with open_handle_gzip(vcf_filename, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, pos, info = cols[0], int(cols[1]), cols[7]
            kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            variant_id = kv.get("variant_id")
            if variant_id is None:
                continue
            end = int(kv["END"]) if "END" in kv else None
            svlen = int(kv["SVLEN"]) if "SVLEN" in kv else None
            variants.append((int(variant_id), chrom, pos, end, svlen))
    return variants


def _score_track(track: ConservationTrack, variants) -> dict[int, float]:
    """ Score every variant against a single bigWig track, opening one handle for this worker. """
    results: dict[int, float] = {}
    bw = pyBigWig.open(track.path)
    try:
        chroms = bw.chroms()
        prefix = "chr" if "chr1" in chroms else ""    # bigWig chrom naming differs per file
        for variant_id, chrom, pos, end, svlen in variants:
            ck = prefix + chrom
            chrom_len = chroms.get(ck)
            if chrom_len is None:
                continue
            s, e = sv_conservation_window(pos, end, svlen)
            s = max(s, 0)
            e = min(e, chrom_len)
            if e <= s:
                e = s + 1
            try:
                value = span_max_aligned(bw, ck, s, e)
            except Exception:
                logging.warning("pyBigWig %s failed for %s:%d-%d", track.name, ck, s, e)
                value = None
            if value is not None:
                results[variant_id] = value
    finally:
        bw.close()
    return results


def conservation_sidecar_filename(vcf_annotated_filename: str) -> str:
    """ Deterministic sidecar path next to the annotated VCF (mirrors the VEP `_warnings.txt`
        pattern), so the import lane can find it without a DB field. """
    return vcf_annotated_filename + ".conservation.tsv"


def write_conservation_sidecar(filename: str, results: dict[int, dict[str, float]],
                               tracks: list[ConservationTrack]) -> None:
    """ Write {variant_id: {db_column: value}} as a TSV keyed by variant_id, one column per track. """
    columns = [t.db_column for t in tracks]
    with open(filename, "wt", encoding="utf-8") as f:
        f.write("\t".join(["variant_id", *columns]) + "\n")
        for variant_id, values in results.items():
            cells = ["" if values.get(c) is None else repr(values[c]) for c in columns]
            f.write("\t".join([str(variant_id), *cells]) + "\n")


def read_conservation_sidecar(filename: str) -> dict[int, dict[str, float]]:
    """ Read a sidecar written by write_conservation_sidecar back into {variant_id: {db_column: value}}. """
    results: dict[int, dict[str, float]] = {}
    with open(filename, "rt", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        columns = header[1:]
        for line in f:
            cells = line.rstrip("\n").split("\t")
            variant_id = int(cells[0])
            values = {c: float(v) for c, v in zip(columns, cells[1:]) if v != ""}
            if values:
                results[variant_id] = values
    return results


def score_sv_vcf(vcf_filename: str, genome_build: GenomeBuild,
                 threads: Optional[int] = None) -> dict[int, dict[str, float]]:
    """ Compute the 4 conservation _max columns for every SV in `vcf_filename`.

        Returns {variant_id: {db_column: value}} - keyed by variant, ready to merge onto the import
        path. The bigWig tracks are scored in parallel across a thread pool (one pyBigWig handle per
        worker); the pool size defaults to settings.ANNOTATION_VEP_FORK (the setting that also drives
        VEP --fork). """
    tracks = get_sv_conservation_tracks(genome_build)
    variants = read_sv_vcf_variants(vcf_filename)
    results: dict[int, dict[str, float]] = {}
    if not tracks or not variants:
        return results

    if threads is None:
        threads = settings.ANNOTATION_VEP_FORK or 1
    max_workers = max(1, min(threads, len(tracks)))

    logging.info("Scoring %d SVs across %d conservation tracks (%d threads)",
                 len(variants), len(tracks), max_workers)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        track_results = executor.map(lambda t: (t, _score_track(t, variants)), tracks)
        for track, track_values in track_results:
            for variant_id, value in track_values.items():
                results.setdefault(variant_id, {})[track.db_column] = value
    return results
