""" Prebuilds the phenotype ambiguous-acronym denylist into Redis, so a page render never pays for it: the build
    scans every ontology term and relation (~230MB), which took 100s on a cold database disk.

    Entry points: build_ambiguous_acronym_denylist_task - enqueued when an OntologyVersion is created
    (annotation.signals.ambiguous_acronym_denylist). """
import celery

from annotation.phenotype_matcher import get_ambiguous_acronym_denylist


@celery.shared_task(queue="db_workers")
def build_ambiguous_acronym_denylist_task():
    get_ambiguous_acronym_denylist()
