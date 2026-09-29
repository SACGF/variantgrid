"""
The registry behind 'manage.py create_fake_data' (@see snpdb.management.commands.create_fake_data).

Each app makes the fake data for its own models in '<app>/fake_data.py' (or a 'fake_data' package): a FakeData
subclass per step, decorated with @register. discover() imports every installed app's module so they register
themselves - nothing here imports an app. Steps name the steps they build on in 'requires', and pass what they made
to later steps through the FakeDataContext rather than through imports.
"""
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

from django.contrib.auth.models import User
from django.core.management.base import OutputWrapper
from django.utils.module_loading import autodiscover_modules

if TYPE_CHECKING:
    from snpdb.models import GenomeBuild, Lab, Trio


@dataclass
class FakeDataContext:
    """ Shared by every step in one run - each step fills in its part for the steps after it """
    genome_build: 'GenomeBuild'
    seed: int
    stdout: OutputWrapper
    genes: list[str] = field(default_factory=list)
    """ Gene symbols the fake variants are placed in ("genes") """
    variant_ids_by_gene: dict[str, list[int]] = field(default_factory=dict)
    """ Annotated variants in those genes ("variants") """
    users: list[User] = field(default_factory=list)
    """ The fake users ("people") """
    labs: list['Lab'] = field(default_factory=list)
    """ The fake labs, germline first ("people") """
    trio: Optional['Trio'] = None
    """ ("trio") """


class FakeData:
    """ One step of fake data. create() is idempotent - running it again, or running a step whose data already
        exists, adds nothing - and delete() removes what create() made """
    name: str = ""
    help: str = ""
    requires: tuple[str, ...] = ()
    """ Names of the steps this builds on, which run first """

    @classmethod
    def add_arguments(cls, parser):
        """ Options for this step, beyond the shared --genome-build / --seed / --delete """

    def create(self, context: FakeDataContext, **options):
        raise NotImplementedError()

    def delete(self, context: FakeDataContext, **options):
        raise NotImplementedError()


_FAKE_DATA: dict[str, type[FakeData]] = {}


def register(klass: type[FakeData]) -> type[FakeData]:
    _FAKE_DATA[klass.name] = klass
    return klass


def discover() -> dict[str, type[FakeData]]:
    autodiscover_modules("fake_data")
    return dict(_FAKE_DATA)


def in_dependency_order(names: Iterable[str], fake_data: dict[str, type[FakeData]]) -> list[str]:
    """ The named steps and everything they require, each after its requirements """
    ordered = []

    def visit(name: str, path: tuple[str, ...]):
        if name in path:
            raise ValueError(f"Fake data steps require each other: {' -> '.join((*path, name))}")
        if name not in ordered:
            for required in fake_data[name].requires:
                visit(required, (*path, name))
            ordered.append(name)

    for step_name in names:
        visit(step_name, ())
    return ordered


def zipf_weight(rank: int) -> float:
    """ Long tail - a few genes/records get most of the action, the rest get a little """
    return 1 / (rank + 1) ** 0.8


def in_preferred_order(preferred: list[str], available: list[str]) -> list[str]:
    """ The available genes a step would rather use, in its order - or all of them when it has no preference
        among them (eg the three fake genes on a database with no real annotation) """
    return [gene for gene in preferred if gene in available] or list(available)
