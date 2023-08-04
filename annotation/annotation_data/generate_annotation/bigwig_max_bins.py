#!/usr/bin/env python3
"""
Conservation BigWigs take ages to run on large variants (>10kb) - and we just want the max anyway

So we'll process the file into bins
"""

import os
import pyBigWig
import numpy as np

from argparse import ArgumentParser

def name_from_filename(filename):
    """ removes path and extension"""
    return os.path.splitext(os.path.basename(filename))[0]


SUMMARY_STATS = {
    "min": np.min,
    "mean": np.mean,
    "max": np.max,
    "std_dev": np.std,
    "median": np.median,
}

def get_args():
    valid_stats_csv = ', '.join(sorted(SUMMARY_STATS))
    parser = ArgumentParser(description="Bin BigWig files")
    parser.add_argument("--summary-stat", default='max', type=str,
                        help=f"Summary statistic - must be one of '{valid_stats_csv}'")
    parser.add_argument("--binsize", default=500, type=int, help="BinSize")
    parser.add_argument("bigwig")
    args = parser.parse_args()
    if not args.bigwig.endswith(".bw"):
        parser.error("bigwig file must end with '.bw'")
    if args.summary_stat not in SUMMARY_STATS:
        parser.error(f'summary-stat must be one of {valid_stats_csv}')

    return args

def main(args):
    name = name_from_filename(args.bigwig)
    output_filename = f"{name}.{args.summary_stat}.binsize.{args.binsize}.bw"
    print(f"writing to {output_filename}")

    summary_stat = SUMMARY_STATS[args.summary_stat]

    with pyBigWig.open(args.bigwig) as input_bw:
        # Get the list of chromosomes and their lengths from the header
        chroms = input_bw.chroms()

        with pyBigWig.open(output_filename, "wb") as output_bw:
            # Write the header with the same chromosome lengths
            output_bw.addHeader(list(chroms.items()))

            # Process each chromosome
            for chrom, length in chroms.items():
                chroms = []
                starts = []
                ends = []
                summary_stats = []

                for start in range(0, length, args.binsize):
                    end = min(start + args.binsize, length)

                    values = input_bw.values(chrom, start, end)

                    # None represent regions missing data
                    values = [v for v in values if v is not None]

                    if values:
                        chroms.append(chrom)
                        starts.append(start)
                        ends.append(end)
                        summary_stats.append(summary_stat(values))

                if chroms:
                    output_bw.addEntries(chroms, starts, ends=ends, values=summary_stats)


if __name__ == "__main__":
    args = get_args()
    main(args)
