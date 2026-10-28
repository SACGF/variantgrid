"""
`vg css unused --rendered`: settles the dynamic bucket of `vg css unused` (names like `cs-{{ status }}` whose
prefix is built up at runtime) by rendering pages and collecting every class that appears in a `class=`
attribute. A dynamic name seen in rendered HTML is live; one never seen needs a reader (or a canary), since the
box may just lack the data that produces it, or JS may add it after load.

Owns the crawl. Seeds are the first-party URL patterns that reverse with no arguments and whose view names a
template; each rendered page's links (LINK_ATTRS) add pages whose URL name has been rendered fewer than
`per_name` times, so detail pages (an allele, a classification, a sample) are reached through real data, a few
of each. A URL name that reads as an action (ACTION_RE) is never requested. Every request runs as the
`claude_agent` user inside a transaction that is rolled back, under a statement timeout, with celery's
send_task, Rollbar, outgoing POSTs (Slack), email and subprocesses stubbed, so a view that writes, queues,
reports or shells out (annotation_build_detail runs VEP on a test VCF) leaves nothing behind.

Entry points: `crawl_rendered_classes` (returns a CrawlResult) and `render_rendered_report`.
"""
import re
from collections import Counter, deque
from dataclasses import dataclass, field
from urllib.parse import urlsplit

from bs4 import BeautifulSoup
from django.db import transaction
from django.test import Client
from django.urls import NoReverseMatch, Resolver404, get_resolver, resolve, reverse

from library.django_utils.database_utils import pg_settings
from library.vg.css import UnusedReport
from library.vg.maps.urls import _is_api, _templates, _view_callable, _walk
from library.vg.page import get_agent_user, quiet_side_effects
from library.vg.repo import first_party_packages

# URL names that do something rather than show something; a GET to one is never made
ACTION_RE = re.compile(r"delete|remove|create|clone|activate|load|cancel|rerun|retry|reload|upload|download|export"
                       r"|logout|sync|send|trigger|save|set_|toggle|reset|restart|import|refresh|queue|run_|graph"
                       r"|wizard|redirect|liftover|annotate|populate|merge|batch|bulk", re.IGNORECASE)
LINK_ATTRS = ("href", "data-url", "data-href")
STATEMENT_TIMEOUT = "30s"


@dataclass
class CrawlResult:
    classes: set[str] = field(default_factory=set)
    rendered: Counter = field(default_factory=Counter)  # URL name -> pages rendered
    failed: dict[str, str] = field(default_factory=dict)  # url -> status or exception


def page_url_names() -> set[str]:
    """ First-party, non-API URL names whose view renders a template: the pages worth requesting """
    names = set()
    for name, path, pattern in _walk(get_resolver()):
        view = _view_callable(pattern)
        owner = (getattr(view, "__module__", "") or "").split(".")[0]
        if name and owner in first_party_packages() and not _is_api(view, path) and _templates(view) \
                and not ACTION_RE.search(name):
            names.add(name)
    return names


def _seeds(names: set[str]) -> list[str]:
    seeds = []
    for name in sorted(names):
        try:
            seeds.append(reverse(name))
        except NoReverseMatch:  # needs arguments
            continue
    return seeds


def _links(soup: BeautifulSoup) -> list[str]:
    links = []
    for attr in LINK_ATTRS:
        for element in soup.find_all(attrs={attr: True}):
            parts = urlsplit(element[attr])
            if parts.scheme in ("", "http", "https") and parts.netloc in ("", "localhost") and parts.path.startswith("/"):
                links.append(parts.path)
    return links


def crawl_rendered_classes(username: str | None = None, max_pages: int = 500, per_name: int = 2,
                           progress=None) -> CrawlResult:
    result = CrawlResult()
    names = page_url_names()
    client = Client(SERVER_NAME="localhost", raise_request_exception=False)
    queue = deque(_seeds(names))
    seen_paths = set(queue)
    with quiet_side_effects():
        with transaction.atomic():
            client.force_login(get_agent_user(username) if username else get_agent_user())
            transaction.set_rollback(True)  # the session lives in the cache; last_login need not be kept
        while queue and sum(result.rendered.values()) < max_pages:
            path = queue.popleft()
            try:
                url_name = resolve(path).view_name
            except Resolver404:
                continue
            if url_name not in names or result.rendered[url_name] >= per_name:
                continue
            result.rendered[url_name] += 1
            if progress:
                progress(path)
            try:
                with transaction.atomic(), pg_settings(local=True, statement_timeout=STATEMENT_TIMEOUT):
                    response = client.get(path)
                    transaction.set_rollback(True)
            except Exception as e:  # pylint: disable=broad-exception-caught
                # a page that errors tells us nothing about classes; keep crawling
                result.failed[path] = f"{type(e).__name__}: {e}"[:200]
                continue
            if response.status_code != 200 or "text/html" not in response.get("Content-Type", ""):
                result.failed[path] = str(response.status_code)
                continue
            soup = BeautifulSoup(response.content, "html.parser")
            for element in soup.find_all(class_=True):
                result.classes.update(element["class"])
            for link in _links(soup):
                if link not in seen_paths:
                    seen_paths.add(link)
                    queue.append(link)
    return result


def render_rendered_report(report: UnusedReport, crawl: CrawlResult) -> str:
    seen = [s for s in report.dynamic if s.name in crawl.classes]
    unseen = [s for s in report.dynamic if s.name not in crawl.classes]
    out = [f"Rendered {sum(crawl.rendered.values())} pages across {len(crawl.rendered)} URL names "
           f"({len(crawl.failed)} failed or not HTML), {len(crawl.classes)} distinct classes",
           f"\n{len(seen)} dynamic names seen in rendered HTML (live): " + " ".join(s.spelled for s in seen),
           f"\n{len(unseen)} dynamic names never rendered - JS-added, or need data this box lacks; read the code:"]
    for selector in unseen:
        where = " ".join(f"{sheet}:{line}" for sheet, line in selector.locations)
        out.append(f"  {selector.spelled:40} {where}")
    return "\n".join(out)
