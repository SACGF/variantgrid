class DiscordanceReportResolution:
    ONGOING = None
    CONCORDANT = 'C'
    CONTINUED_DISCORDANCE = 'D'

    CHOICES = (
        (CONCORDANT, 'Concordant'),
        (CONTINUED_DISCORDANCE, 'Continued Discordance')
    )

class ContinuedDiscordanceReason:
    UNRESPONSIVE = 'U'
    DIFFERENT_CURATION_METHODS = 'D'

    CHOICES = (
        (UNRESPONSIVE, 'One or more labs was unresponsive'),
        (DIFFERENT_CURATION_METHODS, 'Different curation methods are being employed')
    )
