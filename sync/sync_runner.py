from library.log_utils import report_exc_info
from sync.models.models import SyncDestination, SyncRun
from sync.shariant.shariant_download import sync_shariant_download
from sync.shariant.shariant_upload import sync_shariant_upload


def run_sync(sync_destination: SyncDestination, full_sync: bool = False) -> SyncRun:
    """
    :param sync_destination: which sync are we running
    :param full_sync: if True all eligible records will be synced regardless of if they've been synced before, otherwise only changes

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
    try:
        config = sync_destination.config
        sync_type = config.get('type')
        if sync_type in ('shariant', 'variantgrid'):
            direction = config.get('direction')
            if direction == 'upload':
                return sync_shariant_upload(sync_destination=sync_destination, full_sync=full_sync)
            if direction == 'download':
                return sync_shariant_download(sync_destination=sync_destination, full_sync=full_sync)
            raise ValueError('config.direction must be upload or download')
        raise ValueError(f'unknown sync_type: {sync_type}')
    except:
        report_exc_info(extra_data={"sync_name": sync_destination.name})
        raise
