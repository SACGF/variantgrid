
def get_sequencer_model_from_name(name):

    # This code was taken from Illuminate library - https://bitbucket.org/invitae/illuminate
    model = None
    if name.startswith("NS") or name.startswith("NB"):
        model = "NextSeq 500"
    elif name.startswith("M"):
        model = "MiSeq"
    elif name.startswith("D"):
        model = "HiSeq 2500"
    elif name.startswith("SN"):
        model = "HiSeq 2000"
    # elif machine_id.startswith("??"):
    # model = "Hiseq 3000"
    elif name.startswith("K"):
        model = "HiSeq 4000"
    elif name.startswith("ST"):
        model = "HiSeq X"

    return model
