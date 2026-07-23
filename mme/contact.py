""" Resolve the MatchMaker Exchange `contact` block for a classification.

The MME `contact` (name / href / institution) is the public identity of a shared
record. It is sourced from the classification's owning Lab, falling back to the
server-wide `settings.MME_CONTACT`. Both the outbound submission builder and the
inbound serve path share this single resolver so every record is attributed to the
lab that produced it (not to the node operator).

A lab contact is only "valid" when it yields both `name` and `href` (matching MME's
existing name+href requirement); otherwise the resolver drops to the next tier.
"""
from django.conf import settings


def _mailto(email: str) -> str:
    """ Single point where a contact email becomes a `mailto:` href. """
    return f"mailto:{email.strip()}"


def _looks_like_url(value: str) -> bool:
    return value.startswith("http://") or value.startswith("https://")


def lab_mme_contact(lab) -> dict:
    """ Build {name, href, institution} from a Lab's contact fields.
        Returns {} unless it yields both `name` and `href`. """
    if lab is None:
        return {}

    name = (lab.contact_name or "").strip() or (lab.name or "").strip()

    email = (lab.contact_email or "").strip()
    url = (lab.url or "").strip()
    if email:
        href = _mailto(email)
    elif _looks_like_url(url):
        href = url
    else:
        href = ""

    if not (name and href):
        return {}

    organization = getattr(lab, "organization", None)
    org_name = (organization.name or "").strip() if organization else ""
    institution = org_name or (lab.name or "").strip()

    return {"name": name, "href": href, "institution": institution}


def settings_mme_contact() -> dict:
    """ settings.MME_CONTACT when it carries both name and href, else {}. """
    contact = settings.MME_CONTACT or {}
    if contact.get("name") and contact.get("href"):
        return dict(contact)
    return {}


def mme_contact_for_classification(classification) -> dict:
    """ Ordered fallback chain: lab contact -> server contact -> {}.
        An Organization tier can be inserted here later without touching call sites. """
    lab = getattr(classification, "lab", None)
    return lab_mme_contact(lab) or settings_mme_contact()


# Generic name for the shared resolver: the owning-lab {name, href, institution} of a
# shared classification record. Beacon's record-tier attribution (§5.5) uses this same
# chain so every record is attributed to the lab that produced it.
lab_contact = mme_contact_for_classification
