"""
review's enums as one namespace, outside the models package so importing one doesn't run review/models/__init__.
Callers write `from review.enums import ReviewMedium, ReviewParticipants`.
"""
from review.enums.review_enums import *
