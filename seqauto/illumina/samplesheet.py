"""

Code taken from Jinghua (Frank) Feng's bwa_gatk_freebayes_annotate_pipeline

"""

from io import StringIO
import logging
import os

import pandas as pd

ILLUMINA = 'ILLUMINA'
ILLEGAL_CHARACTERS = "?()[]/\=+<>:;\"',*^|&. "  # From Illumina docs
NO_SPACES_COLUMNS = ['sample_id', 'sample_name']


def is_illumina_convention(fastq):
    if fastq.endswith('_R1_001.fastq.gz') or fastq.endswith('_R2_001.fastq.gz'):
        return True
    return False


def is_simplified_single_end_convention(fastq):
    if not is_illumina_convention(fastq):
        if not (fastq.endswith('_R1.fastq.gz') or fastq.endswith('_R2.fastq.gz')):
            return True
    return False


def is_new_sheet(sheet):
    with open(sheet) as f:
        first_line = f.readline().strip()
        return first_line.startswith('[Header]')

def csv_file_from_new_samplesheet(sheet):
    mem_file = StringIO()
    with open(sheet) as f:
        found_data = False
        for line in f:
            if found_data:
                mem_file.write(line)
            else:
                if line.strip().startswith('[Data]'):
                    found_data = True

    mem_file.seek(0)
    return mem_file


def is_standard_valid_flowcell_dir(sequencing_run_path, date_on_file: str = None):
    filename_parts = sequencing_run_path.split('_')
    if date_on_file is None:
        date_on_file = filename_parts[0]

    if len(date_on_file) != 6:
        logging.warning('SequencingRun: %s - date from file name section "%s" not in format YY-MM-DD: ',
                        sequencing_run_path, date_on_file)
        return False

    # There are 2 types of naming conventions:
    # 140924_M02027_0090_000000000-A8N7P
    # 150731_NB501009_0003_AHF7YMBGXX
    VALID_FLOWCELL_NAMES = ['A', 'B']
    flow_cell_id = filename_parts[-1]
    p = flow_cell_id.find("-")
    if p < 0:
        start = flow_cell_id[0]
        flow_cell_id = flow_cell_id[1:]
    else:
        start = flow_cell_id[p + 1]

    if start not in VALID_FLOWCELL_NAMES:
        logging.warning("SequencingRun: %s doesn't have flowcell ID '%s' start with %s", sequencing_run_path, flow_cell_id, VALID_FLOWCELL_NAMES)
        return False

    return True


def convert_sheet_to_df(sheet, date_on_file: str = None):
    """
    columns in returned df:
        sample_id, sample_name, sample_project, flowcell_id, lane, barcode, date, platform
    """
    platform = ILLUMINA  # for now, it is always ILLUMINA
    flow_cell_id = None  # default None = use column from sample sheet

    sequencing_run_path = os.path.basename(os.path.dirname(sheet))
    is_standard_valid_flowcell_dir(sequencing_run_path, date_on_file=date_on_file)

    filename_parts = sequencing_run_path.split('_')
    if date_on_file is None:
        date_on_file = filename_parts[0]

    if is_new_sheet(sheet):
        # 150828_NB501009_0004_AHFCC3BGXX
        #TODO: Remember in new sheet, lane column is available only for HiSeq and not for MiSeq and NextSeq (confirm this!!),
        #So make the lane all None in df for MiSeq and NextSeq, and use the lane number presented in the fastq file name.
        #Also, flowcell_id is not in the sheet for HiSeq, MiSeq or NextSeq yet -- It will either in
        #file name or inside the sheet, and I am still waiting for Andreas
        #to make a convention for the new sheet.
        csv_file = csv_file_from_new_samplesheet(sheet)
        sample_id = 'Sample_ID'  # Goddamn you Illumina this is quality scores all over again.

        flow_cell_id = filename_parts[-1]
    else:  # old version sheet e.g.
        # ${BIOINFO}/scripts/pipelines/mapping_and_variant_calling/test/sample_sheets/150521_SN1101_0170_AC69B8ACXX.csv
        csv_file = sheet
        sample_id = 'SampleID'

    if len(date_on_file) != 6:
        raise ValueError(f"{sheet}: date indicated in the file name ({date_on_file}) not in format YY-MM-DD")

    yy, mm, dd = date_on_file[:2], date_on_file[2:4], date_on_file[4:]
    if int(mm) > 12 or int(dd) > 31:
        raise ValueError(f"{sheet}: date indicated in the file name ({date_on_file}) not in format YY-MM-DD")

    date = '20{yy}-{mm}-{dd}T00:00:00+0930'.format(yy=yy, mm=mm, dd=dd)

    df = pd.read_csv(csv_file, index_col=None, comment='#', skip_blank_lines=True, skipinitialspace=True, dtype=str)
    # Remove empty rows - sometimes left by people editing SampleSheets.csv by hand
    df = df.dropna(axis=0, how='all')
    df['sample_id'] = df[sample_id].str.strip()
    df['sample_name'] = None
    if 'Sample_Name' in df.columns:
        df['sample_name'] = df['Sample_Name'].str.strip()

    no_sample_name_mask = df['sample_name'].isnull()
    df.loc[no_sample_name_mask, 'sample_name'] = df.loc[no_sample_name_mask, 'sample_id']

    if 'Sample_Project' in df.columns:
        sample_project = df['Sample_Project'].str.strip()
    else:
        sample_project = None
    df['sample_project'] = sample_project

    if 'Failed' in df.columns:
        failed = df['Failed'].str.strip().str.upper().isin(["TRUE", "Y", "YES", "FAIL", "FAILED"])
    else:
        failed = False
    df['failed'] = failed

    if not flow_cell_id:
        flow_cell_id = df['FCID'].str.strip()

    df['flowcell_id'] = flow_cell_id
    if 'Lane' in df.columns:
        lane = df['Lane'].str.strip()
    else:
        lane = None
    df['lane'] = lane
    index_names = ["Index", 'index']
    barcode = None
    for i in index_names:
        if i in df.columns:
            barcode = df[i]
            break
    df['barcode'] = barcode
    df['date'] = date
    df['platform'] = platform

    return df


def samplesheet_is_valid(sample_sheet_df):
    """ returns is_valid, error message """
    for column in NO_SPACES_COLUMNS:
        col = sample_sheet_df[column].dropna()  # Handle crazy NaN spreadsheets (saved as CSV with blank rows/cols)

        for row in col:
            for char in row:
                if char in ILLEGAL_CHARACTERS:
                    msg = "Cannot handle illegal character or whitespace in the sample sheet\n"
                    msg += "Column: %s, row value: '%s', illegal char = '%s' \n" % (column, row, char)
                    return False, msg

    barcode_lengths = sample_sheet_df["barcode"].str.len().value_counts()
    if len(barcode_lengths) != 1:
        msg = f"Barcodes must be all the same size, was: {barcode_lengths}"
        return False, msg

    return True, None
