import os

import xml.etree.ElementTree as ET


def get_single_element(run_info_file, root, element_name):
    values = list(root.iter(element_name))
    if len(values) != 1:
        msg = f"Not exactly 1  <{element_name}> tag in {run_info_file}"
        raise ValueError(msg)
    return values[0]


def get_instrument_from_run_info(flowcell_dir):
    run_info_file = os.path.join(flowcell_dir, "RunInfo.xml")
    tree = ET.parse(run_info_file)
    root = tree.getroot()
    return get_single_element(run_info_file, root, "Instrument")


def get_run_parameters(run_parameters_dir):
    """ returns instrument, experiment_name """

    # Damn you Illumina! MiSeq has lower-case 'r'
    RUN_PARAMETERS_FILES = ["RunParameters.xml", "runParameters.xml"]
    run_info_file = None
    for f in RUN_PARAMETERS_FILES:
        potential_run_info_file = os.path.join(run_parameters_dir, f)
        if os.path.exists(potential_run_info_file):
            run_info_file = potential_run_info_file
            break

    if run_info_file is None:
        msg = f"Couldn't find {' or '.join(RUN_PARAMETERS_FILES)} in {run_parameters_dir}"
        raise ValueError(msg)

    tree = ET.parse(run_info_file)
    root = tree.getroot()
    experiment = get_single_element(run_info_file, root, "ExperimentName")
    try:
        instrument = get_single_element(run_info_file, root, "InstrumentID")
    except:
        instrument = get_instrument_from_run_info(run_parameters_dir)

    return instrument.text, experiment.text
