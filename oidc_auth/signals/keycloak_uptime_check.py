from typing import Optional

from django.dispatch import receiver
from django.conf import settings

from library.keycloak import Keycloak
from library.uptime_check import uptime_check_signal, UptimeCheckResponse, UptimeCheckStatus


@receiver(signal=uptime_check_signal)
def keycloak_uptime_check(sender, **kwargs):
    if not settings.KEYCLOAK_SYNC_DETAILS:
        return None

    status = UptimeCheckStatus.OKAY
    note: Optional[str] = None
    try:
        Keycloak().ping()
    except Exception as ex:
        print(str(ex))
        status = UptimeCheckStatus.CRITICAL_FAILURE
        note = str(ex)

    return UptimeCheckResponse(
        name="KeyCloak",
        status=status,
        note=note
    )
