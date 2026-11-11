#!/usr/bin/env python3
"""
Adds the sites gnomAD's T2T liftover is missing (#2788): where the CHM13 reference carries GRCh38's ALT allele.

gnomAD's v4.1 T2T liftover (gnomad.*.v4.1.sites.GCA_009914755.4.trimmed_liftover.vcf.gz) only lifts records whose REF
matches CHM13. Where CHM13 carries GRCh38's ALT allele instead, the record is dropped - so a T2T variant there (which
is GRCh38's reference allele, often common) gets no gnomAD AF. This lifts those records from the GRCh38 joint file
with REF/ALT swapped, for gnomad_t2t_ref_swaps.sh to merge into the T2T file:

    1. sites   - positions in the hg38ToHs1 chain blocks where GRCh38 and CHM13 have different bases
    2. (bcftools view -T the sites file over the GRCh38 gnomAD file)
    3. swap    - the GRCh38 records whose ALT is the CHM13 base, as T2T records with REF/ALT swapped (skipping any
                 the T2T liftover already has)

The GRCh38 REF allele's record has AF = 1 - (every SNV ALT's AF), overall, per genetic ancestry group and for XY,
with grpmax recomputed from the groups - for a biallelic site AF = 1 - AF and AC = AN - AC. nhomalt and the filtering allele
frequencies can't be derived from the GRCh38 records so are left missing. The site's other SNV ALT alleles (also
dropped by the liftover, their REF no longer matching) are lifted with their own counts. SNVs only - an indel's swap
needs re-normalising on CHM13.

Needs numpy + pysam (both in VariantGrid's virtualenv).
"""
import gzip
import sys
from argparse import ArgumentParser
from collections import Counter, defaultdict
from collections.abc import Iterator

import numpy as np
import pysam

SUB_POPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"]
# grpmax excludes bottlenecked groups - as gnomad_data.py:BOTTLENECKED_SUB_POPS
GRPMAX_POPS = ["afr", "amr", "eas", "nfe", "sas"]
# INFO order of gnomad_data.py's T2T combined file, so swapped records read the same as lifted ones
OUTPUT_INFO = (["AF", "AC", "AN"] + [f"{p}_{pop}" for pop in
                                     ["afr", "amr", "asj", "eas", "fin", "nfe", "sas", "remaining", "mid"]
                                     for p in ("AF", "AC", "AN")]
               + ["nhomalt", "faf95", "faf99", "fafmax_faf95_max", "fafmax_faf99_max", "non_par",
                  "AC_XY", "AN_XY", "AF_XY", "gnomad_filtered", "AF_grpmax", "AC_grpmax", "AN_grpmax", "grpmax"])
PRIMARY_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
COMPLEMENT_BYTES = np.arange(256, dtype=np.uint8)
COMPLEMENT_BYTES[ACGT] = np.frombuffer(b"TGCA", dtype=np.uint8)


def read_chrom_mapping(filename) -> dict[str, str]:
    """ snpdb/genome/chrom_mapping_<build>.map: chrom name -> contig accession """
    with open(filename) as f:
        return dict(line.split() for line in f if line.strip())


def read_chains(chain_filename) -> Iterator[tuple]:
    """ (tName, tStart, qName, qSize, qStrand, [(tStart, qStart, size)]) per chain """
    with gzip.open(chain_filename, "rt") as f:
        chain = None
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "chain":
                if chain:
                    yield chain
                t, q = int(parts[5]), int(parts[10])
                chain = (parts[2], int(parts[5]), parts[7], int(parts[8]), parts[9], [])
                continue
            size = int(parts[0])
            chain[5].append((t, q, size))
            if len(parts) == 3:
                t += size + int(parts[1])
                q += size + int(parts[2])
        if chain:
            yield chain


def _sequence(fasta: pysam.FastaFile, contig: str, start=None, end=None) -> np.ndarray:
    """ Upper case bases as bytes - a chromosome is ~250MB, so a couple of copies at most """
    return np.frombuffer(fasta.fetch(contig, start, end).upper().encode(), dtype=np.uint8)


def write_sites(args):
    """ GRCh38 contig, pos, ref, T2T contig, pos, the CHM13 base on GRCh38's strand, strand - sorted by GRCh38 """
    grch38_contigs = read_chrom_mapping(args.grch38_chrom_mapping)
    t2t_contigs = read_chrom_mapping(args.t2t_chrom_mapping)
    chroms = args.chrom or PRIMARY_CHROMS
    chains_by_target = defaultdict(list)
    for chain in read_chains(args.chain):
        t_chrom, q_chrom = chain[0].removeprefix("chr"), chain[2].removeprefix("chr")
        if t_chrom in chroms and q_chrom in PRIMARY_CHROMS:
            chains_by_target[t_chrom].append(chain)

    grch38_fasta = pysam.FastaFile(args.grch38_fasta)
    t2t_fasta = pysam.FastaFile(args.t2t_fasta)
    with open(args.output, "w") as out:
        for t_chrom in chroms:
            t_seq = _sequence(grch38_fasta, grch38_contigs[t_chrom])
            rows = []
            for _, _, q_name, q_size, q_strand, blocks in chains_by_target[t_chrom]:
                q_chrom = q_name.removeprefix("chr")
                # Only the span the chain covers (in its strand's coordinates), so a small chain onto another
                # chromosome doesn't load that whole chromosome
                q_lo = blocks[0][1]
                q_hi = blocks[-1][1] + blocks[-1][2]
                if q_strand == "-":
                    # Blocks are on the reverse strand: compare against the reverse complement
                    span = _sequence(t2t_fasta, t2t_contigs[q_chrom], q_size - q_hi, q_size - q_lo)
                    span = COMPLEMENT_BYTES[span[::-1]]
                else:
                    span = _sequence(t2t_fasta, t2t_contigs[q_chrom], q_lo, q_hi)
                for t_start, q_start, size in blocks:
                    t_block = t_seq[t_start:t_start + size]
                    q_block = span[q_start - q_lo:q_start - q_lo + size]
                    differ = (t_block != q_block) & np.isin(t_block, ACGT) & np.isin(q_block, ACGT)
                    for offset in np.flatnonzero(differ):
                        q_pos = q_start + offset
                        if q_strand == "-":
                            q_pos = q_size - q_pos - 1
                        rows.append((t_start + offset + 1, chr(t_block[offset]), q_chrom, q_pos + 1,
                                     chr(q_block[offset]), q_strand))
            rows.sort()
            for t_pos, t_base, q_chrom, q_pos, q_base, q_strand in rows:
                out.write("\t".join([grch38_contigs[t_chrom], str(t_pos), t_base,
                                     t2t_contigs[q_chrom], str(q_pos), q_base, q_strand]) + "\n")
            print(f"{t_chrom}: {len(rows):,} differing sites", file=sys.stderr)


def _to_t2t(allele: str, strand: str) -> str:
    """ GRCh38-strand allele on the T2T forward strand """
    return allele.translate(COMPLEMENT) if strand == "-" else allele


def _info_dict(info: str) -> dict[str, str]:
    d = {}
    for kv in info.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v if _ else None  # None = flag
    return d


def _int(v):
    return None if v in (None, "", ".") else int(v)


COUNT_SUFFIXES = [""] + [f"_{pop}" for pop in SUB_POPS] + ["_XY"]


def _t2t_info(filter_value: str, grch38_info: dict) -> dict:
    """ OUTPUT_INFO keys, missing unless set """
    out = dict.fromkeys(OUTPUT_INFO, ".")
    if "non_par" in grch38_info:
        out["non_par"] = None
    else:
        del out["non_par"]
    out["gnomad_filtered"] = "0" if filter_value in ("PASS", ".") else "1"
    return out


def _set_grpmax(out: dict):
    grpmax = [(int(out[f"AC_{pop}"]) / int(out[f"AN_{pop}"]), pop) for pop in GRPMAX_POPS
              if out[f"AN_{pop}"] != "." and int(out[f"AN_{pop}"])]
    if grpmax:
        af, pop = max(grpmax)
        out.update({"AF_grpmax": repr(af), "AC_grpmax": out[f"AC_{pop}"], "AN_grpmax": out[f"AN_{pop}"],
                    "grpmax": pop})


def swapped_info(swap_info: dict, other_infos: list[dict], filter_value: str) -> dict:
    """ INFO for the GRCh38 REF allele: 1 - the frequencies of every SNV ALT allele, over the CHM13 allele's AN.
        gnomAD's joint counts are per allele (an allele seen only in genomes has the genomes' AN), so the
        alleles' frequencies are summed rather than their counts """
    out = _t2t_info(filter_value, swap_info)
    for suffix in COUNT_SUFFIXES:
        an = _int(swap_info.get(f"AN{suffix}"))
        if not an:
            continue
        alt_af = 0.0
        for info in [swap_info, *other_infos]:
            ac, other_an = _int(info.get(f"AC{suffix}")), _int(info.get(f"AN{suffix}"))
            if ac and other_an:
                alt_af += ac / other_an
        ref_af = min(max(1.0 - alt_af, 0.0), 1.0)
        out[f"AC{suffix}"] = str(round(ref_af * an))
        out[f"AN{suffix}"] = str(an)
        out[f"AF{suffix}"] = f"{ref_af:.6f}"
    _set_grpmax(out)
    return out


def lifted_info(info: dict, filter_value: str) -> dict:
    """ INFO for another GRCh38 ALT allele at the site - its own counts, now against CHM13's base """
    out = _t2t_info(filter_value, info)
    for key in out:
        if key not in ("non_par", "gnomad_filtered") and info.get(key) not in (None, ""):
            out[key] = info[key]
    return out


def _format_info(info: dict) -> str:
    return ";".join(k if v is None else f"{k}={v}" for k, v in info.items())


def write_swaps(args):
    """ Reads the sites file and the GRCh38 gnomAD records at those sites (VCF text on stdin), writes the swapped
        T2T records (unsorted, no header) and prints counts """
    sites = defaultdict(list)  # (contig, pos) -> [(ref, T2T contig, T2T pos, CHM13 base, strand)]
    with open(args.sites) as f:
        for line in f:
            contig, pos, ref, q_contig, q_pos, q_base, strand = line.rstrip("\n").split("\t")
            sites[(contig, int(pos))].append((ref, q_contig, int(q_pos), q_base, strand))

    counts = Counter()
    existing = pysam.TabixFile(args.existing_vcf)

    def in_existing(q_contig, q_pos, t2t_ref, t2t_alt) -> bool:
        return any(line.split("\t", 5)[3:5] == [t2t_ref, t2t_alt]
                   for line in existing.fetch(q_contig, q_pos - 1, q_pos))

    with open(args.output, "w") as out:
        for (contig, pos), records in _records_by_position(sys.stdin):
            for ref, q_contig, q_pos, q_base, strand in sites.get((contig, pos), []):
                snvs = [r for r in records if r[3] == ref and len(r[4]) == 1]
                swap = [r for r in snvs if r[4] == q_base]
                if not swap:
                    counts["no GRCh38 record of the CHM13 allele"] += 1
                    continue
                record = swap[0]
                others = [r for r in snvs if r is not record]
                t2t_ref = _to_t2t(q_base, strand)
                if not _int(_info_dict(record[7]).get("AN")):
                    counts["skipped - AN=0"] += 1
                    continue
                new_info = swapped_info(_info_dict(record[7]), [_info_dict(r[7]) for r in others], record[6])
                if in_existing(q_contig, q_pos, t2t_ref, _to_t2t(ref, strand)):
                    counts["already in the T2T liftover"] += 1
                    continue
                out.write("\t".join([q_contig, str(q_pos), record[2], t2t_ref, _to_t2t(ref, strand), ".", ".",
                                     _format_info(new_info)]) + "\n")
                counts["swapped (GRCh38 REF allele)"] += 1
                counts[f"site with {len(others)} other SNV alleles"] += 1
                for other in others:
                    t2t_alt = _to_t2t(other[4], strand)
                    if in_existing(q_contig, q_pos, t2t_ref, t2t_alt):
                        counts["other allele already in the T2T liftover"] += 1
                        continue
                    out.write("\t".join([q_contig, str(q_pos), other[2], t2t_ref, t2t_alt, ".", ".",
                                         _format_info(lifted_info(_info_dict(other[7]), other[6]))]) + "\n")
                    counts["other SNV allele lifted"] += 1
    for k, v in sorted(counts.items()):
        print(f"{k}: {v:,}", file=sys.stderr)


def _records_by_position(lines) -> Iterator[tuple[tuple[str, int], list[list[str]]]]:
    current, records = None, []
    for line in lines:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        key = (cols[0], int(cols[1]))
        if key != current:
            if records:
                yield current, records
            current, records = key, []
        records.append(cols)
    if records:
        yield current, records


def main():
    parser = ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sites = sub.add_parser("sites", help="Positions where GRCh38 and CHM13 bases differ")
    sites.add_argument("--chain", required=True, help="hg38ToHs1.over.chain.gz")
    sites.add_argument("--grch38-fasta", required=True)
    sites.add_argument("--t2t-fasta", required=True)
    sites.add_argument("--grch38-chrom-mapping", required=True)
    sites.add_argument("--t2t-chrom-mapping", required=True)
    sites.add_argument("--chrom", action="append", help="GRCh38 chromosome (repeatable, default all)")
    sites.add_argument("--output", required=True)
    swap = sub.add_parser("swap", help="Swapped T2T records from GRCh38 records (VCF on stdin)")
    swap.add_argument("--sites", required=True)
    swap.add_argument("--existing-vcf", required=True, help="gnomAD's T2T liftover (combined) - records it has are left alone")
    swap.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "sites":
        write_sites(args)
    else:
        write_swaps(args)


if __name__ == "__main__":
    main()
