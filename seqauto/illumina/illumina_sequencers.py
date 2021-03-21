
# Regex tested on :
#     MiSeq: 150130_M01761_0114_000000000-ACUR0
#     NextSeq: 150202_NS500318_0047_AH3KLMBGXX
#     HiSeq 2000: 130919_SN792_0281_BD2CHRACXX
#     HiSeq 2500: 150203_D00535_0052_AC66RWANXX
#     HiSeq 4000: 150210_K00111_0013_AH2372BBXX
#     HiSeq X: 141121_ST-E00107_0356_AH00C3CCXX
SEQUENCING_RUN_REGEX = r"\d{6}[_-](NS|NB|M|D|SN|K|ST)(.{3,7})_\d{4}_(0{9}-.{5}|.{10})"


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
