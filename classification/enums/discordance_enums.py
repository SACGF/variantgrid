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
