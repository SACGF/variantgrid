import enum
from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import Optional, Dict, Any

from django.utils import timezone
from lazy import lazy
from model_utils.models import now

from genes.models import GeneSymbol
from ontology.models import OntologyTermRelation, OntologyTerm, OntologyImport, OntologyService


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
        else:
            self.updates += 1
            return self.updates

    def __str__(self) -> str:
        return f"Inserts {self.inserts}, Updates {self.updates}, Stale {self.deletes}"


class OntologyBuilder:
    """
    Not the most efficient way of getting data in (no bulk inserts are used)
    But simplest way to update existing values and create new ones
    """
    class CreatedState(enum.IntEnum):
        """
        Bit overkill but returning a boolean from a few methods results in unclear
        """
        STUB = 1
        DETAILED = 2

    def __init__(self, filename: str, context: str, import_source: str, processor_version: int = 1, force_update: bool = False):
        """
        :filename: Name of the resource used, only used for logging purposes
        :context: In what context are we inserting records (e.g. full file context, or partial panel app)
        :ontology_service: Service that sources this data, but import is not restricted to this service
        :force_update: If true, the ensure methods don't raise exceptions
        """
        self.start = datetime.now()
        self.duration = None
        self.filename = filename if ':' in filename else filename.split("/")[-1]
        self.import_source = import_source
        self.context = context
        self.processor_version = processor_version
        self.data_hash = None
        self.force_update = force_update

        self.previous_import: Optional[OntologyImport] = OntologyImport.objects.filter(import_source=import_source, context=context, completed=True).order_by('-modified').first()
        if self.previous_import and self.previous_import.processor_version != self.processor_version:
            self.previous_import = None  # if previous import was done with an older version of the import code, don't count it

        self.total_count = 0
        self.counters: Dict[Any, OperationCounter] = defaultdict(OperationCounter)
        self.created_cache: Dict[str, int] = dict()
        self.gene_symbols: Dict[str, GeneSymbol] = dict()
        for gene_symbol in GeneSymbol.objects.all():
            self.gene_symbols[gene_symbol.symbol] = gene_symbol

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

    def complete(self, purge_old=True):
        """
        :purge_old If True will mark OntologyTermGeneRelations not included in this import (but included in a previous import
        with the same context and ontology service) as deleted
        Will also complain about other records not included in this import but wont delete them
        """
        if purge_old:
            for model in [OntologyTermRelation, OntologyTerm]:
                olds = model.objects.filter(from_import__context=self.context, from_import__import_source=self.import_source).exclude(from_import=self._ontology_import)
                for old in olds[0:3]:
                    print(f"This appears stale: {old}")

                count = olds.count()
                self.counters[model].deletes += count
                if model == OntologyTermRelation:  # we only delete relations, assume terms are going to persist forever or at work be marked deprecated
                    if count:
                        print("Deleting {count} stale relationships")
                    olds.delete()

        self._ontology_import.completed = True
        self._ontology_import.save()

    def report(self):
        self.duration = datetime.now() - self.start
        print(f"Time taken {self.duration}")
        for model, counter in self.counters.items():
            print(f"{model} - {counter}")

    @lazy
    def _ontology_import(self) -> OntologyImport:
        return OntologyImport.objects.create(
            import_source=self.import_source,
            context=self.context,
            filename=self.filename,
            processed_date=now,
            processor_version=self.processor_version,
            hash=self.data_hash)

    def _count(self, model, created: bool):
        value = self.counters[model].count_op(created)
        self.total_count += 1
        if self.total_count % 1000 == 0:
            print(f"Handled {self.total_count} records")

    def _add_term_stub(self, term_id: str):
        """
        Used internally to ensure relationships can link to existing terms
        """
        if term_id in self.created_cache:
            return term_id
        parts = term_id.split(":")
        ontology_service = parts[0]
        ontology_index = int(parts[1])

        OntologyTerm.objects.get_or_create(
            id=term_id,
            ontology_service=ontology_service,
            index=ontology_index,
            defaults={"from_import": self._ontology_import}
        )
        self.created_cache[term_id] = OntologyBuilder.CreatedState.STUB
        return term_id

    def add_term(self, term_id: str, name: str, definition: Optional[str], extra: Dict = None, primary_source: bool = True):
        """
        Returns the term_id
        """
        if self.created_cache.get(term_id) == OntologyBuilder.CreatedState.DETAILED:
            return

        term_id = term_id.replace("_", ":")
        if not primary_source:
            if existing := OntologyTerm.objects.filter(id=term_id).first():
                existing_import = existing.from_import
                if existing_import.context != self.context or existing_import.import_source != self.import_source:
                    # record was imported by a different process, so leave it alone
                    # e.g. it was imported via an OMIM import but we're just referencing a stub value
                    # now from a MONDO import.
                    self.created_cache[term_id] = OntologyBuilder.CreatedState.DETAILED
                    return

        parts = term_id.split(":")
        ontology_service = parts[0]
        ontology_index = int(parts[1])

        defaults = {
            "ontology_service": ontology_service,
            "index": ontology_index,
            "from_import": self._ontology_import
        }
        # if there is an existing record, and it has values for name, definition or extra and we don't
        # don't want to overwrite them
        if name:
            defaults["name"] = name
        if definition:
            defaults["definition"] = definition
        if extra is not None:
            defaults["extra"] = extra

        term, created = OntologyTerm.objects.update_or_create(
            id=term_id,
            defaults=defaults
        )
        self._count(OntologyTerm, created)
        self.created_cache[term_id] = OntologyBuilder.CreatedState.DETAILED
        return term

    def add_ontology_relation(self, source_term_id: str, dest_term_id: str, relation: str, extra: Dict = None):
        self._add_term_stub(source_term_id)
        self._add_term_stub(dest_term_id)

        defaults = {
            "from_import": self._ontology_import
        }
        if extra is not None:
            defaults["extra"] = extra

        relation, created = OntologyTermRelation.objects.update_or_create(
            source_term_id=source_term_id,
            dest_term_id=dest_term_id,
            relation=relation,
            defaults=defaults
        )
        self._count(OntologyTermRelation, created)
        return relation
