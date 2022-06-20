from dataclasses import dataclass
from datetime import timedelta, datetime
from enum import Enum
from typing import Optional, Dict, List, TypeVar, Generic, Iterable, Type

from django.db.models import Model
from django.utils import timezone
from lazy import lazy
from model_utils.models import now

from ontology.models import OntologyTermRelation, OntologyTerm, OntologyImport


class OntologyBuilderDataUpToDateException(Exception):
    pass


@dataclass
class OperationCounter:
    inserts = 0
    updates = 0
    deletes = 0

    def count_op(self, created: bool):
        if created:
            self.inserts += 1
            return self.inserts
        self.updates += 1
        return self.updates

    def __str__(self) -> str:
        return f"Inserts {self.inserts}, Updates {self.updates}, Stale {self.deletes}"


@dataclass(frozen=True)
class RelationKey:
    source: str
    dest: str
    relation: str

    @staticmethod
    def from_existing(ot: OntologyTermRelation) -> 'RelationKey':
        return RelationKey(ot.source_term_id, ot.dest_term_id, ot.relation)


class ModifiedStatus(int, Enum):
    EXISTING = 0
    MODIFIED = 1
    CREATED = 2


T = TypeVar("T")


@dataclass
class CachedObj(Generic[T]):
    obj: T
    status: ModifiedStatus = ModifiedStatus.EXISTING

    def modify(self, ontology_import: OntologyImport):
        self.obj.from_import = ontology_import
        if self.status == ModifiedStatus.EXISTING:
            self.status = ModifiedStatus.MODIFIED
            # set modified date so we can use bulk update
            self.obj.modified = timezone.now()

    @staticmethod
    def bulk_apply(model: Type[Model], cache: Iterable['CachedObj'], fields: List[str], verbose=False):
        created = [c.obj for c in cache if c.status == ModifiedStatus.CREATED]
        modified = [c.obj for c in cache if c.status == ModifiedStatus.MODIFIED]
        if verbose:
            print(f"{model} Inserting {len(created):,}")
            print(f"{model} Modifying {len(modified):,}")
        batch_size = 10000
        if created:
            model.objects.bulk_create(created, batch_size=batch_size)
        if modified:
            model.objects.bulk_update(modified, fields, batch_size=batch_size)


class OntologyBuilder:

    def __init__(self, filename: str, context: str, import_source: str, processor_version: int = 1, force_update: bool = False):
        """
        :filename: Name of the resource used, only used for logging purposes
        :context: In what context are we inserting records (e.g. full file context, or partial panel app)
        :ontology_service: Service that sources this data, but import is not restricted to this service
        :force_update: If true, the ensure methods don't raise exceptions
        """
        self.start = datetime.now()
        self.filename = filename if ':' in filename else filename.split("/")[-1]
        self.import_source = import_source
        self.context = context
        self.processor_version = processor_version
        self.data_hash = None
        self.force_update = force_update

        self.previous_import: Optional[OntologyImport] = OntologyImport.objects.filter(import_source=import_source, context=context, completed=True).order_by('-modified').first()
        if self.previous_import and self.previous_import.processor_version != self.processor_version:
            self.previous_import = None  # if previous import was done with an older version of the import code, don't count it

        self.full_cache = False
        self.terms: Dict[str, CachedObj[OntologyTerm]] = dict()
        self.relations: Dict[RelationKey, CachedObj[OntologyTermRelation]] = dict()

    def ensure_old(self, max_age: timedelta):
        """
        Raises OntologyBuilderDataUpToDateException if data is up to date
        """
        if self.previous_import:
            if not self.force_update and self.previous_import.processed_date > datetime.now(tz=timezone.get_default_timezone()) - max_age:
                raise OntologyBuilderDataUpToDateException()

    def ensure_hash_changed(self, data_hash: str):
        self.data_hash = data_hash
        if not self.force_update and self.previous_import and self.previous_import.hash == data_hash:
            self.previous_import.processed_date = now
            self.previous_import.save()
            raise OntologyBuilderDataUpToDateException()

    def cache_everything(self):
        """
        Call this to prefetch all data from the previous import, should greatly
        """
        print("About to pre-cache all Ontology")
        for t in OntologyTerm.objects.all():
            self.terms[t.id] = CachedObj(t)

        for tr in OntologyTermRelation.objects.all():
            self.relations[RelationKey.from_existing(tr)] = CachedObj(tr)

        self.full_cache = True
        print("Cache complete")

    def _fetch_term(self, term_id: str) -> CachedObj[OntologyTerm]:
        if pre_cached := self.terms.get(term_id):
            return pre_cached
        if not self.full_cache and (in_db := OntologyTerm.objects.filter(id=term_id).first()):
            cached = CachedObj(obj=in_db)
            self.terms[term_id] = cached
            return cached

        parts = term_id.split(":")
        ontology_service = parts[0]
        ontology_index = int(parts[1])

        cached = CachedObj(OntologyTerm(
            id=term_id,
            ontology_service=ontology_service,
            index=ontology_index,
            from_import=self._ontology_import
        ), ModifiedStatus.CREATED)
        self.terms[term_id] = cached

        return cached

    def _fetch_relation(self, rk: RelationKey) -> CachedObj[OntologyTermRelation]:
        if pre_cached := self.relations.get(rk):
            return pre_cached
        if not self.full_cache and (in_db := OntologyTermRelation.objects.filter(
                source_term_id=rk.source,
                dest_term_id=rk.dest,
                relation=rk.relation
        ).first()):
            cached = CachedObj(in_db)
            self.relations[rk] = cached
            return cached

        self._fetch_term(rk.source)
        self._fetch_term(rk.dest)
        new_r = CachedObj(
            OntologyTermRelation(
                source_term_id=rk.source,
                dest_term_id=rk.dest,
                relation=rk.relation,
                from_import=self._ontology_import
            ),
            ModifiedStatus.CREATED
        )
        self.relations[rk] = new_r
        return new_r

    def add_ontology_relation(self, source_term_id: str, dest_term_id: str, relation: str, extra: Optional[Dict] = None):
        rk = RelationKey(source=source_term_id, dest=dest_term_id, relation=relation)
        cached = self._fetch_relation(rk)
        cached.modify(self._ontology_import)
        cached.obj.extra = extra

    def add_term(self,
                 term_id: str,
                 name: str,
                 definition: Optional[str] = None,
                 extra: Optional[Dict] = None,
                 aliases: Optional[List[str]] = None,
                 primary_source: bool = True):

        cached = self._fetch_term(term_id)
        if not primary_source and cached.status == ModifiedStatus.EXISTING:
            # record was imported by a different process, so leave it alone
            # e.g. it was imported via an OMIM import but we're just referencing a stub value
            # now from a MONDO import.
            return

        cached.modify(self._ontology_import)
        term = cached.obj

        if aliases:
            # want to maintain order (so don't convert to a set)
            unique_aliases = list()
            for alias in aliases:
                if alias not in unique_aliases and alias != name:
                    unique_aliases.append(alias)
            aliases = unique_aliases

        # if primary source of data, overwrite existing record
        # otherwise, only fill in provided data
        if name or primary_source:
            term.name = name

        if definition or primary_source:
            term.definition = definition

        if extra is not None or primary_source:
            term.extra = extra

        if aliases is not None or primary_source:
            term.aliases = aliases or list()

    def complete(self, purge_old_relationships=True, purge_old_terms=False, verbose=True):
        """
        :purge_old If True will mark OntologyTermGeneRelations not included in this import (but included in a previous import
        with the same context and ontology service) as deleted
        Will also complain about other records not included in this import but wont delete them
        """

        CachedObj.bulk_apply(OntologyTerm, self.terms.values(), ["name", "definition", "extra", "aliases", "from_import", "modified"], verbose=verbose)
        CachedObj.bulk_apply(OntologyTermRelation, self.relations.values(), ["extra", "from_import", "modified"], verbose=verbose)

        # Now to find previous imports - and their terms that weren't updated by this import (and purge them if requested)
        old_imports = set(OntologyImport.objects.filter(context=self.context, import_source=self.import_source).values_list("pk", flat=True))
        if self._ontology_import.pk in old_imports:
            old_imports.remove(self._ontology_import.pk)

        for model in [OntologyTermRelation, OntologyTerm]:

            olds = model.objects.filter(from_import__in=old_imports)
            count = olds.count()
            if (model == OntologyTermRelation and purge_old_relationships) or \
                    (model == OntologyTerm and purge_old_terms):  # we only delete relations, assume terms are going to persist forever or at work be marked deprecated
                if count:
                    if verbose:
                        print(f"{model} Deleting {count:,}")
                    olds.delete()
        if verbose:
            time_taken = datetime.now() - self.start
            print(f"Bulk complete in {time_taken}")

        self._ontology_import.completed = True
        self._ontology_import.save()

    @lazy
    def _ontology_import(self) -> OntologyImport:
        return OntologyImport.objects.create(
            import_source=self.import_source,
            context=self.context,
            filename=self.filename,
            processed_date=now,
            processor_version=self.processor_version,
            hash=self.data_hash)
