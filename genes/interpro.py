"""
Pfam protein domains, fetched per-gene from the InterPro REST API and cached in PfamDomains.

@see https://www.ebi.ac.uk/interpro/api/

EBI stopped publishing the per-organism Pfam file the bulk ingest used, and the remaining bulk sources
cost hours to fetch 90k proteins to serve the few thousand anyone looks at. A gene is a median of 3
UniProt sequences, so we fetch a gene's worth at a time, the first time someone views it.
@see https://github.com/SACGF/variantgrid/issues/1554
"""
import logging
import time
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from http import HTTPStatus

from django.conf import settings
from django.db import transaction
from django.utils import timezone
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from genes.models import Pfam, PfamDomains, PfamSequence, PfamSequenceIdentifier

# "count" is how many Pfam entries match the protein, and the API default page size is 20 - ask for
# enough that a protein matching many families doesn't need paging
INTERPRO_PAGE_SIZE = 200

_retry = Retry(total=5, read=5, backoff_factor=0, allowed_methods=["GET"])
_session = Session()
_session.mount("http://", HTTPAdapter(max_retries=_retry))
_session.mount("https://", HTTPAdapter(max_retries=_retry))


def _lazy_domains_enabled() -> bool:
    return settings.PFAM_INTERPRO_LAZY_DOMAINS and not settings.UNIT_TEST


def _entry_url(seq_id: str) -> str:
    return f"{settings.PFAM_INTERPRO_API_URL}/entry/pfam/protein/uniprot/{seq_id}/?page_size={INTERPRO_PAGE_SIZE}"


def _retrieve_pfam_entries(seq_id: str) -> list[dict]:
    """ Runs in a worker thread (no DB access) - returns [] for a protein with no Pfam match """
    entries = []
    url = _entry_url(seq_id)
    while url:
        r = _session.get(url, timeout=settings.PFAM_INTERPRO_TIMEOUT_SECONDS)
        if r.status_code == HTTPStatus.NO_CONTENT:
            break
        r.raise_for_status()
        data = r.json()
        entries.extend(data.get("results") or [])
        url = data.get("next")
    return entries


def _pfam_domains(seq_id: str, entries: list[dict], pfam_ids: set[int]) -> tuple[list[PfamDomains], int]:
    """ :returns domains (1 per fragment, so a discontinuous domain draws as separate boxes),
                 number of Pfam accessions we have no Pfam record for """
    domains = []
    num_missing_pfam = 0
    for entry in entries:
        accession = entry["metadata"]["accession"]
        try:
            pfam_id = Pfam.get_pk_from_accession(accession)
        except ValueError:
            pfam_id = None
        if pfam_id is None or pfam_id not in pfam_ids:
            num_missing_pfam += 1
            continue

        for protein in entry.get("proteins") or []:
            for location in protein.get("entry_protein_locations") or []:
                for fragment in location.get("fragments") or []:
                    domains.append(PfamDomains(pfam_sequence_id=seq_id, pfam_id=pfam_id,
                                               start=fragment["start"], end=fragment["end"]))
    return domains, num_missing_pfam


@transaction.atomic
def _store_domains(seq_id: str, domains: list[PfamDomains]):
    """ Replace rather than add to, so a sequence holding stale rows converges """
    PfamDomains.objects.filter(pfam_sequence_id=seq_id).delete()
    if domains:
        PfamDomains.objects.bulk_create(domains)
    PfamSequence.objects.filter(pk=seq_id).update(domains_imported=timezone.now())


def store_domains_for_sequences(seq_ids: Iterable[str]) -> int:
    """ Fetches domains for any of these PfamSequences not already imported.

        Never raises - a sequence whose fetch failed is left unstamped and retried on the next view,
        so a graph renders without a domain track rather than erroring.

        :returns number of PfamDomains created """
    if not _lazy_domains_enabled():
        return 0

    unimported = list(PfamSequence.objects.filter(pk__in=list(seq_ids),
                                                  domains_imported__isnull=True).values_list("pk", flat=True))
    if not unimported:
        return 0

    pfam_ids = set(Pfam.objects.values_list("pk", flat=True))
    num_domains = 0
    num_missing_pfam = 0
    num_imported = 0
    start = time.monotonic()

    executor = ThreadPoolExecutor(max_workers=settings.PFAM_INTERPRO_MAX_WORKERS)
    try:
        futures = {executor.submit(_retrieve_pfam_entries, seq_id): seq_id for seq_id in unimported}
        try:
            for future in as_completed(futures, timeout=settings.PFAM_INTERPRO_DEADLINE_SECONDS):
                seq_id = futures[future]
                try:
                    entries = future.result()
                except Exception as e:
                    logging.warning("InterPro: couldn't retrieve domains for '%s': %s", seq_id, e)
                    continue
                domains, missing_pfam = _pfam_domains(seq_id, entries, pfam_ids)
                _store_domains(seq_id, domains)
                num_domains += len(domains)
                num_missing_pfam += missing_pfam
                num_imported += 1
        except TimeoutError:
            # Whatever hasn't come back stays unstamped for the next view to pick up
            logging.warning("InterPro: %d of %d sequences done in %ds, leaving the rest for next time",
                            num_imported, len(unimported), settings.PFAM_INTERPRO_DEADLINE_SECONDS)
    finally:
        executor.shutdown(wait=False, cancel_futures=True)

    if num_missing_pfam:
        logging.warning("InterPro: skipped %d domains with no Pfam record", num_missing_pfam)
    logging.debug("InterPro: %d domains for %d sequences in %.1fs", num_domains, num_imported,
                  time.monotonic() - start)
    return num_domains


def store_domains_for_transcripts(transcript_ids: Iterable[str]) -> int:
    """ Fetches domains for every UniProt sequence these transcripts map to """
    if not _lazy_domains_enabled():
        return 0

    seq_ids = PfamSequenceIdentifier.objects.filter(transcript__in=list(transcript_ids)) \
        .values_list("pfam_sequence_id", flat=True)
    return store_domains_for_sequences(set(seq_ids))
