import enum
from datetime import timedelta, datetime
from typing import Optional, Dict

from django.utils import timezone
from lazy import lazy
from model_utils.models import now

from genes.models import GeneSymbol
from ontology.models import OntologyTermRelation, OntologyTermGeneRelation, OntologyTerm, OntologyImport, OntologySet


class OntologyBuilderDataUpToDateException(Exception):
    pass


class OntologyBuilder:
    """
    Not the most efficient way of getting data in (no bulk inserts are used)
    But simplest way to update existing values and create new ones
    """
    class BuildState(enum.IntEnum):
        """
        Bit overkill but returning a boolean from a few methods results in unclear
        """
        NO_UPDATE_REQUIRED = 0
        UPDATE_REQUIRED = 1

    def __init__(self, filename: str, context: str, ontology_set: OntologySet):
        self.filename = filename.split("/")[-1]
        self.ontology_set = ontology_set
        self.context = context
        self.data_hash = None
        self.previous_import: Optional[OntologyImport] = OntologyImport.objects.filter(ontology_set=ontology_set, context=context, completed=True).order_by('-modified').first()
        self.insert_count = 0
        self.update_count = 0
        self.delete_count = 0
        self.created_cache: Dict[str, OntologyTerm] = set()

    def ensure_old(self, max_age: timedelta):
        """
        Raises OntologyBuilderDataUpToDateException if data is up to date
        """
        if self.previous_import:
            if self.previous_import.processed_date > datetime.now(tz=timezone.get_default_timezone()) - max_age:
                raise OntologyBuilderDataUpToDateException()

    def ensure_hash_changed(self, data_hash: str):
        self.data_hash = data_hash
        if self.previous_import and self.previous_import.hash == data_hash:
            self.previous_import.processed_date = now
            self.previous_import.save()
            raise OntologyBuilderDataUpToDateException()

    def complete(self, purge_old=True):
        """
        Returns how many old relationships were deleted from a previous run
        """
        if purge_old:
            for model in [OntologyTermGeneRelation, OntologyTermRelation, OntologyTerm]:
                olds = model.objects.filter(from_import__context=self.context, from_import__ontology_set=self.ontology_set).exclude(from_import=self._ontology_import)
                for old in olds[0:3]:
                    print(f"This appears old {old}")

                if model == OntologyTermGeneRelation: #there are the only things that should actually get marked as deleted
                    olds.update(deletion_date=now, modified=now)
                    self.delete_count += olds.count()

        self._ontology_import.completed = True
        self._ontology_import.save()

    @lazy
    def _ontology_import(self) -> OntologyImport:
        return OntologyImport.objects.create(ontology_set=self.ontology_set, context=self.context, filename=self.filename, processed_date=now, hash=self.data_hash)

    def _count(self, created: bool):
        if created:
            self.insert_count += 1
        else:
            self.update_count += 1

    def add_term(self, term_id: str, name: str = None, definition: str = None, extra: Dict = None, stub: bool = False) -> OntologyTerm:
        term_id = term_id.replace("_", ":")
        if stub:
            if existing := OntologyTerm.objects.filter(id=term_id).first():
                existing_import = existing.from_import
                if existing_import.context != self.context or existing_import.ontology_set != self.ontology_set:
                    # record was imported by a different process, so leave it alone
                    # e.g. it was imported via an OMIM miport but we're just referencing a stub value
                    # now from a MONDO import.
                    return existing

                # nothing else to do to this record, skip for an optimisation
                if not name and not definition and not extra and existing.from_import == self._ontology_import:
                    return existing

        parts = term_id.split(":")
        ontology_set = parts[0]
        ontology_index = int(parts[1])

        defaults = {
            "ontology_set": ontology_set,
            "index": ontology_index,
            "from_import": self._ontology_import
        }
        # if there is an existing record, and it has values for name, definition or extra and we don't
        # don't want to overwrite them
        if name:
            defaults["name"] = name
        if definition:
            defaults["definition"] = definition
        if extra:
            defaults["extra"] = extra

        term, created = OntologyTerm.objects.update_or_create(
            id=term_id,
            defaults=defaults
        )
        self._count(created)
        return term

    def add_gene_relation(self, term_id: str, gene_symbol_str: str, relation: str, extra = None) -> Optional[OntologyTermGeneRelation]:
        """
        Overrides any previous relationship - only one term -> gene relationship can exist regardless of the relation
        If the gene symbol doesn't exist, nothing will happen
        """
        if gene_symbol := GeneSymbol.objects.filter(symbol=gene_symbol_str).first():
            key = f"{term_id}/{gene_symbol_str}"
            relation_obj, created = OntologyTermGeneRelation.objects.update_or_create(
                term=self.add_term(term_id=term_id, stub=True),
                gene_symbol=gene_symbol,
                relation=relation,
                defaults={
                    "extra": extra,
                    "from_import": self._ontology_import
                }
            )
            self._count(created)
            return relation_obj

        print(f"Could not find gene symbol {gene_symbol_str}")
        return None

    def add_ontology_relation(self, source_term_id: str, dest_term_id: str, relation: str):
        relation, created = OntologyTermRelation.objects.get_or_create(
            source_term=self.add_term(term_id=source_term_id, stub=True),
            dest_term=self.add_term(term_id=dest_term_id, stub=True),
            relation=relation,
            defaults={
                "from_import": self._ontology_import
            }
        )
        self._count(created)
        return relation
