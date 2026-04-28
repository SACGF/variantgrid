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
  - SAMPLE_CT = number of contributing records.

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
import sys
from collections import defaultdict
from pathlib import Path

CHROM_ALIASES = ["Chr", "Chromosome", "CHROM", "#Chr"]
POS_ALIASES = ["Position", "Pos", "POS"]
REF_ALIASES = ["Reference", "Ref", "REF"]
ALT_ALIASES = ["Alternate", "Alt", "ALT"]
STUDY_ALIASES = ["StudyName"]
PUBMED_ALIASES = ["PubmedID", "PubMedID", "PMID"]
PHENO_ALIASES = ["PrimaryPhenotype"]

VALID_BASES = set("ACGT")


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


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    records = defaultdict(lambda: {"studies": [], "pubmeds": [], "phenotypes": []})
    skipped_indel = 0
    skipped_malformed = 0
    total_rows = 0

    for tsv_path in argv[1:]:
        sys.stderr.write(f"Reading {tsv_path}\n")
        with smart_open(tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            chrom_c = find_col(reader.fieldnames, CHROM_ALIASES)
            pos_c = find_col(reader.fieldnames, POS_ALIASES)
            ref_c = find_col(reader.fieldnames, REF_ALIASES)
            alt_c = find_col(reader.fieldnames, ALT_ALIASES)
            study_c = find_col(reader.fieldnames, STUDY_ALIASES, required=False)
            pubmed_c = find_col(reader.fieldnames, PUBMED_ALIASES, required=False)
            pheno_c = find_col(reader.fieldnames, PHENO_ALIASES, required=False)

            for row in reader:
                total_rows += 1
                chrom = normalize_chrom(row.get(chrom_c) or "")
                pos = (row.get(pos_c) or "").strip()
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
    out.write('##INFO=<ID=StudyName,Number=.,Type=String,'
              'Description="Contributing denovo-db studies, &-separated parallel array (one element per contributing record)">\n')
    out.write('##INFO=<ID=PubmedID,Number=.,Type=String,'
              'Description="PubMed IDs of contributing records, &-separated parallel array">\n')
    out.write('##INFO=<ID=PrimaryPhenotype,Number=.,Type=String,'
              'Description="Primary phenotype per contributing record, &-separated parallel array">\n')
    out.write('##INFO=<ID=SAMPLE_CT,Number=1,Type=Integer,'
              'Description="Number of denovo-db records aggregated at this site">\n')
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
        info_parts.append(f"SAMPLE_CT={len(bucket['studies'])}")
        out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{';'.join(info_parts)}\n")

    sys.stderr.write(
        f"Done: {total_rows} TSV rows -> {len(records)} unique variants "
        f"(skipped {skipped_indel} non-ACGT, {skipped_malformed} malformed)\n"
    )


if __name__ == "__main__":
    main(sys.argv)
