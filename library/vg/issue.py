"""
`vg issue <ref>`: the whole state of a GitHub issue (or PR) before any work starts. The issue, its comments, and every
piece of existing work linked to it: PRs that reference it with their branch and state, remote branches named for it
(`<issue>_<slug>`), a plan under claude/plans/ and commits on master that cite it. The PostToolUse hook on Bash
(`.claude/hooks/post_bash.py`) appends the linked-work line to any `gh issue view`, so an open PR is seen even when
this command is skipped.

Everything goes through `gh` and `git`; nothing here needs Django. Only SACGF/variantgrid is public: for an issue in a
private repo (handed over by the user, so reading that one is allowed) linked issues and PRs are listed by number and
state only, since their titles are other private issues' content (AGENTS.md, "Which repo an issue goes in").
"""
import json
import re
import subprocess
from dataclasses import asdict, dataclass, field
from typing import Optional

from library.vg.repo import REPO_ROOT, git

PUBLIC_REPOS = frozenset({"SACGF/variantgrid"})
DEFAULT_REPO = "SACGF/variantgrid"
GH_TIMEOUT_SECONDS = 30
MAX_COMMITS = 20
MAX_UNREFERENCED_PRS = 100

_URL_RE = re.compile(r"(?:github\.com|\brepos)/(?P<repo>[\w.-]+/[\w.-]+)/(?:issues|pull)/(?P<number>\d+)")
_REPO_HASH_RE = re.compile(r"^(?P<repo>[\w.-]+/[\w.-]+)#(?P<number>\d+)$")
_NUMBER_RE = re.compile(r"^#?(?P<number>\d+)$")


@dataclass(frozen=True)
class IssueRef:
    repo: str
    number: int

    @property
    def public(self) -> bool:
        return self.repo in PUBLIC_REPOS

    @property
    def url(self) -> str:
        return f"https://github.com/{self.repo}/issues/{self.number}"

    def __str__(self):
        return f"{self.repo}#{self.number}"


@dataclass
class LinkedWork:
    """ A PR or issue whose body or comments mention the issue (GitHub's cross-referenced timeline event) """
    kind: str  # "pr" or "issue"
    repo: str
    number: int
    state: str  # open, closed, merged, draft
    title: Optional[str] = None  # withheld for a private repo
    branch: Optional[str] = None

    def __str__(self):
        label = "PR" if self.kind == "pr" else "issue"
        parts = [f"{label} #{self.number}" if self.repo == DEFAULT_REPO else f"{label} {self.repo}#{self.number}",
                 self.state]
        if self.branch:
            parts.append(f"({self.branch})")
        if self.title:
            parts.append(self.title)
        return " ".join(parts)


@dataclass
class Comment:
    author: str
    created: str
    body: str
    url: str


@dataclass
class IssueWork:
    ref: IssueRef
    title: str
    state: str
    author: str
    created: str
    labels: list[str]
    body: str
    is_pull_request: bool
    comments: list[Comment] = field(default_factory=list)
    linked: list[LinkedWork] = field(default_factory=list)
    branches: list[str] = field(default_factory=list)
    plans: list[str] = field(default_factory=list)
    commits: list[str] = field(default_factory=list)
    unreferenced_prs: list[LinkedWork] = field(default_factory=list)
    """ Open PRs citing no issue at all, gathered only when nothing links to the issue: one of them may be its fix """

    @property
    def has_existing_work(self) -> bool:
        return bool(self.linked or self.branches or self.plans or self.commits)


def parse_issue_ref(text: str, default_repo: str = DEFAULT_REPO) -> IssueRef:
    """ '2317', '#2317', 'SACGF/variantgrid_private#12', an issue / PR / comment URL """
    text = text.strip()
    for pattern in (_URL_RE, _REPO_HASH_RE, _NUMBER_RE):
        if match := pattern.search(text):
            return IssueRef(repo=match.groupdict().get("repo") or default_repo, number=int(match["number"]))
    raise ValueError(f"Not an issue reference: {text!r}")


def origin_repo() -> str:
    """ owner/name of this checkout's origin remote, or the public repo when there is none """
    try:
        url = git("remote", "get-url", "origin").strip()
    except subprocess.CalledProcessError:
        return DEFAULT_REPO
    match = re.search(r"github\.com[:/](?P<repo>[\w.-]+/[\w.-]+?)(?:\.git)?$", url)
    return match["repo"] if match else DEFAULT_REPO


def gh_api(endpoint: str):
    """ One page of up to 100 items is plenty for an issue's comments or timeline """
    result = subprocess.run(["gh", "api", endpoint], cwd=REPO_ROOT, check=True, capture_output=True, text=True,
                            timeout=GH_TIMEOUT_SECONDS)
    return json.loads(result.stdout)


def gather_issue(ref: IssueRef, *, comments: bool = True) -> IssueWork:
    issue = gh_api(f"repos/{ref.repo}/issues/{ref.number}")
    work = IssueWork(
        ref=ref,
        title=issue["title"],
        state=_issue_state(issue),
        author=issue["user"]["login"],
        created=issue["created_at"][:10],
        labels=[label["name"] for label in issue.get("labels") or []],
        body=issue.get("body") or "",
        is_pull_request="pull_request" in issue,
        linked=linked_work(ref),
        branches=issue_branches(ref.number),
        plans=issue_plans(ref.number),
        commits=issue_commits(ref.number),
    )
    if not work.has_existing_work:
        work.unreferenced_prs = unreferenced_open_prs(ref.repo if ref.public else DEFAULT_REPO)
    if comments and issue.get("comments"):
        work.comments = [
            Comment(author=c["user"]["login"], created=c["created_at"][:10], body=c.get("body") or "",
                    url=c["html_url"])
            for c in gh_api(f"repos/{ref.repo}/issues/{ref.number}/comments?per_page=100")
        ]
    return work


def linked_work(ref: IssueRef) -> list[LinkedWork]:
    """ Cross-references from the issue's timeline, PRs first, each once; a PR gets its branch from the pulls API """
    seen: dict[tuple[str, int], LinkedWork] = {}
    for event in gh_api(f"repos/{ref.repo}/issues/{ref.number}/timeline?per_page=100"):
        if event.get("event") != "cross-referenced":
            continue
        source = (event.get("source") or {}).get("issue") or {}
        if not source.get("number"):
            continue
        source_repo = (source.get("repository") or {}).get("full_name") or ref.repo
        key = (source_repo, source["number"])
        if key in seen:
            continue
        is_pr = "pull_request" in source
        linked = LinkedWork(kind="pr" if is_pr else "issue", repo=source_repo, number=source["number"],
                            state=_issue_state(source))
        if source_repo in PUBLIC_REPOS:
            linked.title = source.get("title")
        if is_pr:
            _add_pull_details(linked)
        seen[key] = linked
    return sorted(seen.values(), key=lambda w: (w.kind != "pr", -w.number))


_ISSUE_CITATION_RE = re.compile(
    r"(?<![\w/])#\d+|(?<![\w/])[\w.-]+/[\w.-]+#\d+|github\.com/[\w.-]+/[\w.-]+/(?:issues|pull)/\d+")


def unreferenced_open_prs(repo: str) -> list[LinkedWork]:
    """ Open PRs in the public repo whose title and body cite no issue in any repo (AGENTS.md says every PR must) """
    try:
        pulls = gh_api(f"repos/{repo}/pulls?state=open&per_page={MAX_UNREFERENCED_PRS}")
    except (subprocess.SubprocessError, ValueError):
        return []
    found = []
    for pull in pulls:
        if _ISSUE_CITATION_RE.search(f"{pull.get('title') or ''} {pull.get('body') or ''}"):
            continue
        found.append(LinkedWork(kind="pr", repo=repo, number=pull["number"], state="draft" if pull.get("draft") else "open",
                                title=pull.get("title"), branch=(pull.get("head") or {}).get("ref")))
    return found


def _add_pull_details(linked: LinkedWork) -> None:
    try:
        pull = gh_api(f"repos/{linked.repo}/pulls/{linked.number}")
    except (subprocess.SubprocessError, ValueError, KeyError):
        return
    linked.branch = (pull.get("head") or {}).get("ref")
    if pull.get("merged_at"):
        linked.state = "merged"
    elif pull.get("draft"):
        linked.state = "draft"


def _issue_state(issue: dict) -> str:
    if (issue.get("pull_request") or {}).get("merged_at"):
        return "merged"
    return (issue.get("state") or "").lower()


def issue_branches(number: int) -> list[str]:
    """ Remote branches named for the issue by the `<issue>_<slug>` convention """
    pattern = re.compile(rf"^\s*origin/{number}[_-]")
    try:
        lines = git("branch", "-r").splitlines()
    except subprocess.CalledProcessError:
        return []
    return sorted(line.strip() for line in lines if pattern.match(line))


def issue_plans(number: int) -> list[str]:
    plans_dir = REPO_ROOT / "claude" / "plans"
    if not plans_dir.is_dir():
        return []
    return sorted(str(path.relative_to(REPO_ROOT)) for path in plans_dir.glob(f"{number}_*.md"))


def issue_commits(number: int) -> list[str]:
    """ Commits on origin/master (or HEAD without a remote) whose message cites #<number> """
    for ref in ("origin/master", "HEAD"):
        try:
            log = git("log", ref, "--oneline", "--extended-regexp", f"--grep=#{number}([^0-9]|$)",
                      f"--max-count={MAX_COMMITS}")
        except subprocess.CalledProcessError:
            continue
        return log.splitlines()
    return []


def render_issue(work: IssueWork, *, body: bool = True, comments: bool = True) -> str:
    ref = work.ref
    kind = "PR" if work.is_pull_request else "issue"
    lines = [f"{ref}  {kind}  {work.state}  " + (", ".join(work.labels) if work.labels else "no labels"),
             work.title,
             f"by {work.author}, {work.created}  {ref.url}"]
    if not ref.public:
        lines.append("Private repo: this issue only. Linked items are listed by number and state, without titles.")
    lines += ["", "Existing work"]
    lines += [f"  {line}" for line in render_existing_work(work)]
    if work.unreferenced_prs:
        lines += ["", "Open PRs citing no issue at all - one may be this issue's fix, check by hand:"]
        lines += [f"  {pull}" for pull in work.unreferenced_prs]
    if body:
        lines += ["", "Body", work.body.rstrip() or "(empty)"]
    if comments and work.comments:
        lines += ["", f"Comments ({len(work.comments)})"]
        for comment in work.comments:
            lines += [f"--- {comment.author}  {comment.created}  {comment.url}", comment.body.rstrip()]
    return "\n".join(lines) + "\n"


def render_existing_work(work: IssueWork) -> list[str]:
    lines = [str(linked) for linked in work.linked]
    lines += [f"branch {branch}" for branch in work.branches]
    lines += [f"plan {plan}" for plan in work.plans]
    if work.commits:
        lines.append("commits on master citing it:")
        lines += [f"  {commit}" for commit in work.commits]
    return lines or ["none: no PR, branch, plan or master commit references this issue"]


def existing_work_summary(ref: IssueRef) -> str:
    """ One paragraph for the Bash hook: what already exists for the issue, and how to see the whole picture """
    work = IssueWork(ref=ref, title="", state="", author="", created="", labels=[], body="", is_pull_request=False,
                     linked=linked_work(ref), branches=issue_branches(ref.number), plans=issue_plans(ref.number),
                     commits=issue_commits(ref.number))
    if not work.has_existing_work:
        summary = f"No PR, branch, plan or master commit references {ref} yet."
        if unreferenced := unreferenced_open_prs(ref.repo if ref.public else DEFAULT_REPO):
            summary += (" Open PRs citing no issue at all, one of which may be its fix: "
                        + "; ".join(str(pull) for pull in unreferenced) + ".")
        return summary
    items = "; ".join(render_existing_work(work))
    return (f"Existing work for {ref}: {items}. Read it before proposing or writing any "
            f"(scripts/vg issue {ref.number} shows the whole picture).")


def issue_work_json(work: IssueWork) -> str:
    data = asdict(work)
    data["ref"] = {"repo": work.ref.repo, "number": work.ref.number, "url": work.ref.url, "public": work.ref.public}
    return json.dumps(data, indent=1)
