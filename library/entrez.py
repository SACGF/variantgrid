"""
NCBI E-utilities (efetch, esearch, esummary, epost) over requests, with a timeout: Bio.Entrez calls urlopen with
none, so a stalled NCBI hangs the web or celery worker. Entry point: eutils_request - parse what it returns with Bio
(Entrez.read on .content, Medline / SeqIO on .text).

Credentials are Bio.Entrez's module settings (Entrez.email, Entrez.api_key), set from settings in annotation/apps.py.
"""
import threading
import time

from Bio import Entrez
from requests import Response, Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_TIMEOUT_SECS = 20

# NCBI allows 3 requests/second, or 10 with an API key
_MIN_INTERVAL_SECS = 0.37
_MIN_INTERVAL_WITH_API_KEY_SECS = 0.1

# NCBI answers 429 under load even when the rate limit is kept, so retry it with the server errors
# (raise_on_status=False so the caller's raise_for_status reports the actual status, not a RetryError)
_retry = Retry(total=2, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET", "POST"],
               backoff_factor=1, raise_on_status=False)
_session = Session()
_session.mount("https://", HTTPAdapter(max_retries=_retry))

_rate_limit_lock = threading.Lock()
_last_request_time = 0.0


def _wait_for_rate_limit(has_api_key: bool):
    global _last_request_time  # pylint: disable=global-statement
    min_interval = _MIN_INTERVAL_WITH_API_KEY_SECS if has_api_key else _MIN_INTERVAL_SECS
    with _rate_limit_lock:
        wait = _last_request_time + min_interval - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        _last_request_time = time.monotonic()


def eutils_request(utility: str, timeout: float = EUTILS_TIMEOUT_SECS, **params) -> Response:
    """ POSTs (no URL length limit on long id lists) to the E-utility e.g. "efetch", raising for an HTTP error.
        A list/tuple param (e.g. id) is sent comma separated """
    data = {"tool": Entrez.tool}
    if Entrez.email:
        data["email"] = Entrez.email
    if Entrez.api_key:
        data["api_key"] = Entrez.api_key
    for key, value in params.items():
        if isinstance(value, (list, tuple)):
            value = ",".join(str(v) for v in value)
        data[key] = value

    _wait_for_rate_limit(bool(Entrez.api_key))
    response = _session.post(f"{EUTILS_URL}/{utility}.fcgi", data=data, timeout=timeout)
    response.raise_for_status()
    return response
