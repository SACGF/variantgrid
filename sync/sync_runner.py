from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Callable, Dict, List
from library.log_utils import report_exc_info
from sync.models.models import SyncDestination, SyncRun

class SyncRunner(ABC):

    def __init__(self, sync_destination: SyncDestination):
        self.sync_destination = sync_destination

    def get_config(self, param: str, mandatory: bool = True):
        value = self.sync_destination.config.get(param)
        if value is None:
            raise ValueError(f"Expected SyncDestination config property of {param}")
        return value

    @abstractmethod
    def sync(self, full_sync: bool = False, max_rows: Optional[int] = None) -> SyncRun:
        pass


SyncRunnerFactory = Callable[[SyncDestination], SyncRunner]


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
            return factory_requirements.factory(sync_destination)

    raise ValueError(f"No SyncRunner is configured for the config of {sync_destination}")

