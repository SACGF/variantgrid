"""
Settings as `manage.py deployment_check` sees them: the Django system checks (`snpdb/checks.py`, which need
only settings) rerun and reported in the same shape as the other deployment checks, plus the settings that
name a database row, which only the deployed database can confirm.

Entry points: `check_system_checks`, `check_settings_references`.
"""
from django.conf import settings
from django.contrib.auth.models import Group
from django.contrib.sites.models import Site
from django.core.checks import run_checks

from genes.models.models_gene_coverage import CanonicalTranscriptCollection
from genes.models.models_gene_list import GeneList
from snpdb.models.models import Lab


def check_system_checks() -> dict:
    """ Each unsilenced system check Error is an invalid result keyed by its id; a Warning is a warning on
        an otherwise valid one. An Error here stops `upgrade.sh --quick` the way any other invalid check does """
    results = {}
    for message in run_checks():
        if message.is_silenced():
            continue
        key = message.id or message.msg
        text = message.msg + (f" Fix: {message.hint}" if message.hint else "")
        if message.is_serious():
            results[key] = {"valid": False, "fix": text}
        else:
            results.setdefault(key, {"valid": True})["warning"] = text
    return results or {"system checks": {"valid": True}}


def check_settings_references() -> dict:
    """ Settings that hold the pk or group name of a row: a restore from another deployment, or a renamed lab,
        leaves them pointing at nothing """
    results = {}

    def exists(name, qs, fix):
        results[name] = {"valid": qs.exists(), "fix": fix}

    exists("SITE_ID", Site.objects.filter(pk=settings.SITE_ID),
           f"Create django_site row {settings.SITE_ID} (admin, Sites) or change SITE_ID")
    for group_name in settings.LAB_EXTERNAL_MEMBERSHIP:
        exists(f"LAB_EXTERNAL_MEMBERSHIP[{group_name!r}]", Lab.objects.filter(group_name=group_name),
               f"No Lab has group_name {group_name!r} - use the lab's group_name as the key")
    if collection_id := settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID:
        exists("GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID",
               CanonicalTranscriptCollection.objects.filter(pk=collection_id),
               f"No CanonicalTranscriptCollection {collection_id} - import one or set the setting to None")
    if gene_list_id := settings.PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID:
        exists("PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID", GeneList.objects.filter(pk=gene_list_id),
               f"No GeneList {gene_list_id}")
    if settings.SEQAUTO_ENABLED:
        exists("SEQAUTO_API_WRITE_GROUP", Group.objects.filter(name=settings.SEQAUTO_API_WRITE_GROUP),
               f"Create group {settings.SEQAUTO_API_WRITE_GROUP!r} and add the user whose token the pipeline posts with")
    return results
