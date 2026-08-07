class NoTranscript(ValueError):
    """
    Extends ValueError for backwards compatibility.
    Indicates the transcript we are looking for is not in our database
    """


class NoTranscriptVersion(NoTranscript):
    pass


class MissingTranscript(NoTranscript):
    """
    Transcript exists in RefSeq/Ensembl, so c.hgvs (or otherwise) might be okay.
    """


class BadTranscript(NoTranscript):
    """
    Transcript not found in Ensembl or RefSeq (User error)
    """
