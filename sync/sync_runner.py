from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property
from typing import Optional, Callable, Dict, List

from dateutil import tz

from library.log_utils import report_exc_info
from library.oauth import ServerAuth
from sync.models import SyncStatus
from sync.models.models import SyncDestination, SyncRun


class SyncRunInstance:

    def __init__(self,
                 sync_destination: SyncDestination,
                 full_sync: bool = False,
                 max_rows: Optional[int] = None):
        self.sync_destination = sync_destination
        self.full_sync = full_sync
        self.max_rows = max_rows

    def get_config(self, param: str, mandatory: bool = True):
        value = self.sync_destination.config.get(param)
        if value is None:
            raise ValueError(f"Expected SyncDestination config property of {param}")
        return value

    def server_auth(self) -> ServerAuth:
        return ServerAuth.for_sync_details(self.sync_destination.sync_details)

    def last_success(self) -> Optional[SyncRun]:
        return SyncRun.objects.filter(destination=self.sync_destination, status=SyncStatus.SUCCESS).order_by('-created').first()

    def last_success_server_date(self, meta_key='server_date') -> Optional[datetime]:
        if success := self.last_success():
            if meta := success.meta:
                if server_date_str := meta.get(meta_key):
                    return datetime.strptime(server_date_str, "%a, %d %b %Y %H:%M:%S %Z").replace(tzinfo=tz.UTC)

    @cached_property
    def sync_run(self) -> SyncRun:
        return SyncRun.objects.create(destination=self.sync_destination, status=SyncStatus.IN_PROGRESS)

    def run_start(self):
        _ = self.sync_run

    def run_failed(self):
        self.sync_run.status = SyncStatus.FAILED
        self.sync_run.save()

    def run_completed(self, had_records: bool, meta: Optional[Dict] = None):
        self.sync_run.status = SyncStatus.SUCCESS if had_records else SyncStatus.NO_RECORDS
        self.sync_run.meta = meta
        self.sync_run.save()


class SyncRunner(ABC):

    @abstractmethod
    def sync(self, sync_run_instance: SyncRunInstance):
        pass


SyncRunnerFactory = Callable[[], SyncRunner]


@dataclass
class SyncRunnerFactoryRequirements:
    config: Dict
    factory: SyncRunnerFactory

    def matches(self, sync_destination: SyncDestination) -> bool:
        for key, value in self.config.items():
            config_value = sync_destination.config.get(key)
            # TODO allow fuzzier matches
            if isinstance(value, (set, dict, tuple)):
                if config_value not in value:
                    return False
            elif value != config_value:
                return False

        return True


_sync_runner_registry: List[SyncRunnerFactoryRequirements] = []


def register_sync_runner(config: Dict):
    def _wrapper(cls):
        _sync_runner_registry.append(SyncRunnerFactoryRequirements(config, cls))
        return cls
    return _wrapper


def sync_runner_for_destination(sync_destination: SyncDestination) -> SyncRunner:
    config = sync_destination.config
    for factory_requirements in _sync_runner_registry:
        if factory_requirements.matches(sync_destination):
            return factory_requirements.factory()

    raise ValueError(f"No SyncRunner is configured for the config of {sync_destination}")

