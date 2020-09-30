
class MolecularConsequenceColors:
    STOPGAIN = "Stopgain"
    FRAMESHIFT = "Frameshift"
    MISSENSE = "Missense"
    SPLICE_SITE = "Splice site"
    OTHER = "Other"
    NOT_TESTED = "Not tested"

    CONSEQUENCE_LOOKUPS = {
        "missense_variant": MISSENSE,
        "stop_gained": STOPGAIN,
        "frameshift_variant": FRAMESHIFT,
        "splice_acceptor_variant": SPLICE_SITE,
        "splice_donor_variant": SPLICE_SITE,
    }

    VARIANT_TYPES = [
        STOPGAIN,
        FRAMESHIFT,
        MISSENSE,
        SPLICE_SITE,
        OTHER,
        NOT_TESTED,
    ]

    # In order of rank, most important 1st
    ONCOPLOT_COLORS = [
        (STOPGAIN, "#ff0000"),
        (FRAMESHIFT, "#000000"),
        (MISSENSE, "#0000ff"),
        (SPLICE_SITE, "#7030a0"),
        (OTHER, "#ffffff"),
        (NOT_TESTED, "#d9d9d9"),
    ]

    # White works as background but not as hotspot lollypop
    HOTSPOT_COLORS = [
        (STOPGAIN, "#ff0000"),
        (FRAMESHIFT, "#000000"),
        (MISSENSE, "#0000ff"),
        (SPLICE_SITE, "#7030a0"),
        (OTHER, "#d9d9d9"),
    ]
