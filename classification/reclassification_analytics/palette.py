"""
Colours for the reclassification analytics tables and charts. Shared by the result objects, the significance
matrix and the plotly JSON, so they live apart from all three.
"""
from snpdb.models.models_enums import ClinicalSignificance

BENIGN_DIRECTION_COLOUR = ClinicalSignificance.chart_colour(ClinicalSignificance.BENIGN)
PATHOGENIC_DIRECTION_COLOUR = ClinicalSignificance.chart_colour(ClinicalSignificance.PATHOGENIC)
# counts that carry no benign/pathogenic meaning stay outside the significance language
NEUTRAL_COUNT_COLOUR = "#6c757d"
# the matrix shades by how busy a cell is, which is a different thing from significance, so it gets
# its own single hue rather than borrowing the pill colours
MATRIX_HEAT_COLOUR = "#4a6f8a"

EVIDENCE_APPLIED_COLOUR = "#4a8f5b"
EVIDENCE_UNAPPLIED_COLOUR = "#b5564e"
EVIDENCE_STRENGTHENED_COLOUR = "#84b590"
EVIDENCE_WEAKENED_COLOUR = "#d59a95"
EVIDENCE_CHANGED_COLOUR = "#4f7ca6"

# curation effort and calls moving are counted side by side all over the page, and neither of them is
# a benign or pathogenic direction, so they stay out of the significance palette
REVIEW_COLOUR = "#5b8f8a"
CHANGE_COLOUR = "#a8743f"


def translucent(colour: str, alpha: float = 0.45) -> str:
    """ Expands the 4 digit chart colours so the sankey bands can be laid over each other """
    digits = colour.lstrip("#")
    if len(digits) == 3:
        digits = "".join(digit * 2 for digit in digits)
    r, g, b = (int(digits[index:index + 2], 16) for index in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
