#!/usr/bin/env python3
"""
Convert denovo-db TSV file(s) to a VCF suitable for VEP --custom annotation.

The upstream denovo-db VCF only ships STUDIES, SAMPLE_CT and a bunch of redundant
score / AF fields — it lacks the StudyName / PubmedID / PrimaryPhenotype that
make the data clinically useful. The TSV has all three, so we build our own VCF
from the TSV.

Multiple TSV rows for the same variant are merged into a single VCF record:
  - StudyName / PubmedID / PrimaryPhenotype are '&'-joined parallel arrays
    (one element per contributing denovo-db record), preserving duplicates so
    downstream consumers can see "3 of these were autism".
  - CASE_CT / CONTROL_CT = counts of contributing records by case vs. control
    (control is identified by literal PrimaryPhenotype="control"). Filter
    downstream with e.g. `bcftools view -i 'INFO/CASE_CT>0'`.

Output is bare VCF on stdout. Pipe through `bcftools sort -Oz` to get a
properly-sorted compressed VCF.

Usage:
    denovo_db_tsv_to_vcf.py TSV [TSV ...] > out.vcf

Indels: skipped if REF/ALT contain anything other than ACGT (denovo-db
sometimes uses '-' for empty alleles, which is not VCF-legal — fixing that
requires a reference fasta lookup, deferred). A summary is printed to stderr.
"""

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

CHROM_ALIASES = ["Chr", "Chromosome", "CHROM", "#Chr"]
POS_ALIASES = ["Position", "Pos", "POS"]
VARIANT_ALIASES = ["Variant"]  # denovo-db format: "G>C"
REF_ALIASES = ["Reference", "Ref", "REF"]
ALT_ALIASES = ["Alternate", "Alt", "ALT"]
STUDY_ALIASES = ["StudyName"]
PUBMED_ALIASES = ["PubmedID", "PubMedID", "PMID"]
PHENO_ALIASES = ["PrimaryPhenotype"]

VALID_BASES = set("ACGT")
VERSION_RE = re.compile(r"^##version=denovo-db\.v\.(?P<version>\S+)\s*$")
CONTROL_PHENOTYPES = {"control"}  # denovo-db uses the literal lowercase 'control'


def find_col(header, aliases, required=True):
    for a in aliases:
        if a in header:
            return a
    if required:
        raise ValueError(f"None of {aliases} found in TSV header: {header}")
    return None


def smart_open(path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "rt", encoding="utf-8", errors="replace")


def normalize_chrom(s: str) -> str:
    s = s.strip()
    if not s:
        return ""
    if not s.startswith("chr"):
        s = "chr" + s
    if s == "chrMT":
        s = "chrM"
    return s


def vcf_safe(s: str) -> str:
    """Strip characters that break VCF INFO encoding."""
    out = (s or "").replace(";", "_").replace(",", "_").replace("=", "_").replace(" ", "_")
    return out or "."


def _contig_sort_key(c: str):
    """Order chr1..chr22, chrX, chrY, chrM, then anything else alphabetically."""
    base = c[3:] if c.startswith("chr") else c
    if base.isdigit():
        return (0, int(base))
    order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    if base in order:
        return (1, order[base])
    return (2, base)


def parse_version_line(line: str, source: str) -> str:
    m = VERSION_RE.match(line)
    if not m:
        raise ValueError(
            f"{source}: expected first line to be '##version=denovo-db.v.<VERSION>', got: {line!r}")
    return m.group("version")


def split_variant(value: str) -> tuple[str, str]:
    """denovo-db 'Variant' column is 'REF>ALT' (e.g. 'G>C')."""
    if not value or ">" not in value:
        return "", ""
    ref, _, alt = value.partition(">")
    return ref.strip().upper(), alt.strip().upper()


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    records = defaultdict(lambda: {"studies": [], "pubmeds": [], "phenotypes": []})
    skipped_indel = 0
    skipped_malformed = 0
    total_rows = 0
    version = None

    for tsv_path in argv[1:]:
        sys.stderr.write(f"Reading {tsv_path}\n")
        with smart_open(tsv_path) as fh:
            file_version = parse_version_line(fh.readline(), tsv_path)
            if version is None:
                version = file_version
            elif file_version != version:
                raise ValueError(
                    f"Version mismatch: {tsv_path} is {file_version}, earlier file was {version}")

            reader = csv.DictReader(fh, delimiter="\t")
            # Header line is '#SampleID\tStudyName\t...' — strip leading '#' from first field.
            if reader.fieldnames and reader.fieldnames[0].startswith("#"):
                reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
            chrom_c = find_col(reader.fieldnames, CHROM_ALIASES)
            pos_c = find_col(reader.fieldnames, POS_ALIASES)
            variant_c = find_col(reader.fieldnames, VARIANT_ALIASES, required=False)
            ref_c = find_col(reader.fieldnames, REF_ALIASES, required=False)
            alt_c = find_col(reader.fieldnames, ALT_ALIASES, required=False)
            if not variant_c and not (ref_c and alt_c):
                raise ValueError(
                    f"{tsv_path}: need either a 'Variant' column or both 'Ref' and 'Alt' columns. "
                    f"Got: {reader.fieldnames}")
            study_c = find_col(reader.fieldnames, STUDY_ALIASES, required=False)
            pubmed_c = find_col(reader.fieldnames, PUBMED_ALIASES, required=False)
            pheno_c = find_col(reader.fieldnames, PHENO_ALIASES, required=False)

            for row in reader:
                total_rows += 1
                chrom = normalize_chrom(row.get(chrom_c) or "")
                pos = (row.get(pos_c) or "").strip()
                if variant_c:
                    ref, alt = split_variant(row.get(variant_c) or "")
                else:
                    ref = (row.get(ref_c) or "").strip().upper()
                    alt = (row.get(alt_c) or "").strip().upper()

                if not (chrom and pos and ref and alt) or not pos.isdigit():
                    skipped_malformed += 1
                    continue
                if not (set(ref) <= VALID_BASES and set(alt) <= VALID_BASES):
                    skipped_indel += 1
                    continue

                key = (chrom, int(pos), ref, alt)
                bucket = records[key]
                if study_c:
                    bucket["studies"].append((row.get(study_c) or "").strip())
                if pubmed_c:
                    bucket["pubmeds"].append((row.get(pubmed_c) or "").strip())
                if pheno_c:
                    bucket["phenotypes"].append((row.get(pheno_c) or "").strip())

    out = sys.stdout
    out.write("##fileformat=VCFv4.2\n")
    out.write("##source=denovo-db_tsv_to_vcf.py\n")
    out.write("##commandline=" + " ".join(argv) + "\n")
    if version:
        out.write(f"##denovo_db_version={version}\n")
    for contig in sorted({k[0] for k in records}, key=_contig_sort_key):
        out.write(f"##contig=<ID={contig}>\n")
    out.write('##INFO=<ID=StudyName,Number=.,Type=String,'
              'Description="Contributing denovo-db studies, &-separated parallel array (one element per contributing record)">\n')
    out.write('##INFO=<ID=PubmedID,Number=.,Type=String,'
              'Description="PubMed IDs of contributing records, &-separated parallel array">\n')
    out.write('##INFO=<ID=PrimaryPhenotype,Number=.,Type=String,'
              'Description="Primary phenotype per contributing record, &-separated parallel array">\n')
    out.write('##INFO=<ID=CASE_CT,Number=1,Type=Integer,'
              'Description="Number of contributing records with a non-control PrimaryPhenotype">\n')
    out.write('##INFO=<ID=CONTROL_CT,Number=1,Type=Integer,'
              'Description="Number of contributing records with PrimaryPhenotype=control">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for key in sorted(records.keys(), key=lambda k: (k[0], k[1])):
        chrom, pos, ref, alt = key
        bucket = records[key]
        info_parts = []
        if any(bucket["studies"]):
            info_parts.append("StudyName=" + "&".join(vcf_safe(s) for s in bucket["studies"]))
        if any(bucket["pubmeds"]):
            info_parts.append("PubmedID=" + "&".join(vcf_safe(s) for s in bucket["pubmeds"]))
        if any(bucket["phenotypes"]):
            info_parts.append("PrimaryPhenotype=" + "&".join(vcf_safe(s) for s in bucket["phenotypes"]))
        control_ct = sum(1 for p in bucket["phenotypes"] if p.strip().lower() in CONTROL_PHENOTYPES)
        case_ct = len(bucket["phenotypes"]) - control_ct
        info_parts.append(f"CASE_CT={case_ct}")
        info_parts.append(f"CONTROL_CT={control_ct}")
        out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{';'.join(info_parts)}\n")

    sys.stderr.write(
        f"Done: {total_rows} TSV rows -> {len(records)} unique variants "
        f"(skipped {skipped_indel} non-ACGT, {skipped_malformed} malformed)\n"
    )


if __name__ == "__main__":
    main(sys.argv)
