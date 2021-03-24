from collections import defaultdict
import re

from library.utils import FLOAT_REGEX
from seqauto.models.models_enums import SequencerRead, PairedEnd
from snpdb.models import DataState


# TODO: Return a hash, instead of having knowledge of object (should be done in models)
def process_tiles(illumina_flowcell_qc, lines):
    REGEXES = [(r'Mean Cluster Density:\s+(\d+)', 'mean_cluster_density'),
               (r'Mean PF Cluster Density:\s+(\d+)', 'mean_pf_cluster_density'),
               (r'Total Clusters:\s+(\d+)', 'total_clusters'),
               (r'Total PF Clusters:\s+(\d+)', 'total_pf_clusters'),
               (r'Percentage of Clusters PF:\s+' + FLOAT_REGEX, 'percentage_of_clusters_pf'),
               (r'Aligned to PhiX:\s+' + FLOAT_REGEX, 'aligned_to_phix')]

    for line in lines:
        for regex, f in REGEXES:
            if m := re.match(regex, line):
                value = m.group(1)
                setattr(illumina_flowcell_qc, f, value)


def process_quality(illumina_flowcell_qc, lines):
    # Looks like:
    #Read 1: 96.604115
    #Read 2: 96.571730 (Index)
    REGEX = r"Read (\d+): " + FLOAT_REGEX + r"(.*)"
    read_num = 0
    index_num = 0

    READS = [PairedEnd.R1, PairedEnd.R2]
    INDEXES = [SequencerRead.I1, SequencerRead.I2]

    for line in lines:
        if m := re.match(REGEX, line):
            sequencer_read_id = int(m.group(1))
            percent = m.group(2)
            rest_of_line = m.group(3)
            is_index = 'Index' in rest_of_line

            if is_index:
                read = INDEXES[index_num]
                index_num += 1
            else:
                read = READS[read_num]
                read_num += 1

            illumina_flowcell_qc.readq30_set.create(sequencer_read_id=sequencer_read_id,
                                                    percent=percent,
                                                    read=read,
                                                    is_index=is_index)


def process_indexing(illumina_flowcell_qc, lines):
    # Looks like:
    #index_str          project_str  name_str
    #AAGAGGCA+ATAGAGAG  default      XXX        33522694
    #AGGCAGAA+GCGATCTA  default      YYY           19999024

    JUNK = ['index_str',
            'Name']

    for line in lines:
        skip = False
        for j in JUNK:
            if line.startswith(j):
                skip = True
                break
        if skip:
            continue

        cols = line.split()
        if len(cols) == 4:
            illumina_flowcell_qc.illuminaindexqc_set.create(index=cols[0],
                                                            project=cols[1],
                                                            name=cols[2],
                                                            reads=int(cols[3]))


SECTIONS = {'TILE METRICS': process_tiles,
            'QUALITY METRICS': process_quality,
            'INDEXING METRICS': process_indexing,
            'ERROR METRICS': None,
            'CORRECTED INTENSITY': None,
            'EXTRACTION METRICS': None,
            'CONTROL METRICS': None}

ERROR_MESSAGES = [
    "Data file incomplete or unparseable",
    "File not found"
]


def load_from_file(_seqauto_run, illumina_flowcell_qc):
    illumina_flowcell_qc.data_state = DataState.RUNNING
    illumina_flowcell_qc.save()
    try:
        with open(illumina_flowcell_qc.path) as f:
            sections = break_into_sections(f)

        for section, processor in SECTIONS.items():
            if processor:
                lines = sections[section]
                # check lines for errors
                for line in lines:
                    for error_message in ERROR_MESSAGES:
                        if line.endswith(error_message):
                            raise ValueError(error_message)

                processor(illumina_flowcell_qc, lines)

        illumina_flowcell_qc.error_exception = None
        illumina_flowcell_qc.data_state = DataState.COMPLETE
    except Exception as e:
        illumina_flowcell_qc.error_exception = str(e)
        illumina_flowcell_qc.data_state = DataState.ERROR

    illumina_flowcell_qc.save()


def break_into_sections(lines):
    current_section = 'HEADER'
    sections = defaultdict(list)

    for line in lines:
        line = line.strip()
        if line.startswith('---------'):
            continue
        if line in SECTIONS:
            current_section = line

        sections[current_section].append(line)

    return sections
