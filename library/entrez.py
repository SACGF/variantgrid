"""
NCBI E-utilities over requests, with a timeout: Bio.Entrez calls urlopen with none, so a stalled NCBI hangs the
web or celery worker. Entry points: efetch (returns text for Medline / SeqIO / XML parsers), esearch, esummary and
epost (return what Entrez.read parses).

Credentials are Bio.Entrez's module settings (Entrez.email, Entrez.api_key), set from settings in annotation/apps.py.
"""
import threading
import time
from collections.abc import Sequence
from io import BytesIO
from typing import Any, Optional, Union

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


IdsType = Union[str, Sequence[str]]


def _eutils_post(utility: str, timeout: float, **params) -> Response:
    """ POSTs (no URL length limit on long id lists) to the E-utility, raising for an HTTP error.
        None params are left out, a list/tuple (e.g. id) is sent comma separated """
    data = {"tool": Entrez.tool}
    if Entrez.email:
        data["email"] = Entrez.email
    if Entrez.api_key:
        data["api_key"] = Entrez.api_key
    for key, value in params.items():
        if value is None:
            continue
        if isinstance(value, (list, tuple)):
            value = ",".join(str(v) for v in value)
        data[key] = value

    _wait_for_rate_limit(bool(Entrez.api_key))
    response = _session.post(f"{EUTILS_URL}/{utility}.fcgi", data=data, timeout=timeout)
    response.raise_for_status()
    # E-utilities is UTF-8, but text/plain without a charset would otherwise decode as ISO-8859-1
    response.encoding = "utf-8"
    return response


def efetch(db: str, ids: Optional[IdsType] = None, *, rettype: str, retmode: Optional[str] = None,
           webenv: Optional[str] = None, query_key: Optional[str] = None, idtype: Optional[str] = None,
           timeout: float = EUTILS_TIMEOUT_SECS) -> str:
    """ Records by ids, or by a history server (webenv, query_key) from epost """
    response = _eutils_post("efetch", timeout, db=db, id=ids, rettype=rettype, retmode=retmode,
                            WebEnv=webenv, query_key=query_key, idtype=idtype)
    return response.text


def esearch(db: str, term: str, *, retmax: int = 20, timeout: float = EUTILS_TIMEOUT_SECS) -> dict[str, Any]:
    """ Entrez.read result, matching ids are in ["IdList"] """
    response = _eutils_post("esearch", timeout, db=db, term=term, retmax=retmax)
    return Entrez.read(BytesIO(response.content))


def esummary(db: str, ids: IdsType, *, timeout: float = EUTILS_TIMEOUT_SECS) -> list[dict[str, Any]]:
    """ Entrez.read result, one summary per id """
    response = _eutils_post("esummary", timeout, db=db, id=ids)
    return Entrez.read(BytesIO(response.content))


def epost(db: str, ids: IdsType, *, timeout: float = EUTILS_TIMEOUT_SECS) -> dict[str, str]:
    """ Stores ids on the history server, returning {"WebEnv": ..., "QueryKey": ...} to pass to efetch """
    response = _eutils_post("epost", timeout, db=db, id=ids)
    return Entrez.read(BytesIO(response.content))
