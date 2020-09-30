class SyncStatus:
    IN_PROGRESS = 'P'
    SUCCESS = 'S'
    # With no records, we had nothing to upload or download
    # there was no failure, but no success
    NO_RECORDS = 'N'
    FAILED = 'F'
    CHOICES = [
        (IN_PROGRESS, "P"),
        (SUCCESS, "S"),
        (FAILED, "F"),
        (NO_RECORDS, "N")
    ]
