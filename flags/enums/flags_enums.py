"""
FlagStatus: the status a FlagResolution carries (open, closed, rejected); a flag is open while its
resolution's status is OPEN. Stored codes: add values, never change one.
"""


class FlagStatus:
    OPEN = 'O'
    CLOSED = 'C'
    REJECTED = 'R'

    CHOICES = (
        (OPEN, 'Open'),
        (CLOSED, 'Closed'),
        (REJECTED, 'REJECTED'),
    )
