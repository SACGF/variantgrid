from django.conf import settings
import codecs
import logging
import os
import subprocess

from library import genomics
from library.file_utils import mk_path
from library.genomics.bed_file import BedFileReader


def process_bed_file(bed_file, processed_file, has_chr):
    """ formats chromosome, sorts and merges (so it can be used for bed tools intersections) """
    if not os.path.exists(bed_file):
        msg = f"'{bed_file}' does not exist"
        raise IOError(msg)

    mk_path(settings.PROCESSED_BED_FILES_DIR)

    script_name = os.path.join(settings.BASE_DIR, 'scripts', 'bed_sort_and_merge.sh')
    args = [script_name, bed_file]
    bed_pipe = subprocess.Popen(args, stdout=subprocess.PIPE)
    pipe_in = codecs.getreader('utf-8')(bed_pipe.stdout)

    num_records = 0
    with open(processed_file, "w") as out_f:
        for feature in BedFileReader(pipe_in):
            num_records += 1
            chrom = genomics.format_chrom(feature.iv.chrom, has_chr)
            out_f.write('\t'.join((chrom, str(feature.iv.start), str(feature.iv.end))) + '\n')

    logging.info("Wrote %d records to %s", num_records, processed_file)
    return num_records
