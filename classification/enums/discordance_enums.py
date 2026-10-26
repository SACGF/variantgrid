"""
Discordance report vocabulary: resolution, continued discordance reasons, triage status,
next step, and NotifyLevel for discordance notifications.
"""
from enum import Enum

from django.db.models import TextChoices


class DiscordanceReportResolution:
    ONGOING = None  # FIX ME, rename to Active Discordance
    CONCORDANT = 'C'
    CONTINUED_DISCORDANCE = 'D'

    CHOICES = (
        (CONCORDANT, 'Concordant'),
        (CONTINUED_DISCORDANCE, 'Continued Discordance')
    )


class ContinuedDiscordanceReason:
    UNRESPONSIVE = 'U'
    DIFFERENT_CURATION_METHODS = 'D'
    NOT_DEFINED = 'X'

    CHOICES = (
        (UNRESPONSIVE, 'One or more labs was unresponsive'),
        (DIFFERENT_CURATION_METHODS, 'Different curation methods are being employed'),
        (NOT_DEFINED, "See notes")
    )


class DiscordanceReportTriageStatus(TextChoices):
    PENDING = "P", "Pending Triage"
    REVIEWED_WILL_FIX = "F", "Will Amend"
    REVIEWED_WILL_DISCUSS = "D", "For Joint Discussion"
    REVIEWED_SATISFACTORY = "R", "Confident in Classification"
    COMPLEX = "X", "Low Penetrance/Risk Allele etc"


class NotifyLevel(str, Enum):
    NEVER_NOTIFY = "never-notify"
    NOTIFY_IF_CHANGE = "notify-if-change"
    ALWAYS_NOTIFY = "always-notify"


class DiscordanceReportNextStep(str, Enum):
    UNANIMOUSLy_COMPLEX = "C"
    AWAITING_YOUR_TRIAGE = "T"
    AWAITING_YOUR_AMEND = "A"
    AWAITING_OTHER_LAB = "O"
    TO_DISCUSS = "D"
