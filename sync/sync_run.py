from sync.alissa import *  # to get decorators to register
from sync.models import SyncStatus
from sync.shariant import *  # to get decorators to register
from sync.sync_runner import sync_runner_for_destination


def run_sync(sync_destination: SyncDestination, full_sync: bool = False, max_rows: Optional[int] = None) -> SyncRun:
    """
    :param sync_destination: which sync are we running
    :param full_sync: if True all eligible records will be synced regardless of if they've been synced before, otherwise only changes
    :param max_rows: If not None, then the number of rows uploaded will be limited to this (doesn't apply to download)

    Only allowed value currently
    {
        "type":"shariant",
        "direction":"upload",
        "mapping": {
            "labs": {
                # "<lab_group_name>": true | <mapped_group_name>,
                "sa_pathology/national_referral_lab_wch4": true, // map this lab as is
                "sa_pathology/legacy_national_referral_lab_wch4": "sa_pathology/national_referral_lab_wch4" // map legacy to non legacy
                # omitted labs will not be sent
            },
            "share_levels": {
                (logged_in_users|public): (lab|institution|logged_in_users|public)
                "logged_in_users": "lab",
                "public": "institution"
                # only logged_in_users and public will be considered, if omitted they will be sent as is
            },
            "users": {
                "local_username": "remove_username"
                # usernames not provided will be sent as is
            }
        }
    }
    """

    sync_run_instance = SyncRunInstance(sync_destination=sync_destination, full_sync=full_sync, max_rows=max_rows)
    try:
        sync_runner_for_destination(sync_destination).sync(sync_run_instance)
    finally:
        if sync_run_instance.sync_run.status == SyncStatus.IN_PROGRESS:
            sync_run_instance.run_failed()
    return sync_run_instance.sync_run
