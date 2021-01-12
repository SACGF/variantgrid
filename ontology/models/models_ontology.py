import re
from typing import Optional, List, Dict

from django.db import models
from django.db.models import PROTECT, CASCADE, QuerySet
from django.urls import reverse
from model_utils.models import TimeStampedModel

from genes.models import GeneSymbol
from library.utils import Constant

"""
A series of models that currently stores the combination of MONDO and OMIM (with support for HPO if we add an importer to it).
MONDO is given preferrential treatment (e.g. if there's an established relationship between an OMIM term and a gene symbol,
the relationship is actually stored against the MONDO equivilant of the OMIM term).
"""


class OntologyService(models.TextChoices):
    MONDO = "MONDO", "MONDO"
    OMIM = "OMIM", "OMIM"
    HPO = "HPO", "HPO"

    EXPECTED_LENGTHS: Dict[str, int] = Constant({
        MONDO[0]: 7,
        OMIM[0]: 6,
        HPO[0]: 7
    })

    URLS: Dict[str, str] = Constant({
        MONDO[0]: "https://vm-monitor.monarchinitiative.org/disease/MONDO:${1}",
        OMIM[0]: "http://www.omim.org/entry/${1}",
        HPO[0]: "https://hpo.jax.org/app/browse/term/HP:${1}"
    })

    @staticmethod
    def index_to_id(ontology_service: 'OntologyService', index: int):
        expected_length = OntologyService.EXPECTED_LENGTHS.get(ontology_service, 0)
        num_part = str(index).rjust(expected_length, '0')
        return f"{ontology_service}:{num_part}"

    PANEL_APP_AU = "PAAU", "PanelApp AU"


class OntologyRelation(models.TextChoices):
    IS_A = "is_a", "is a"
    """
    Alias relationship is our non-preferred term to preferred term
    If it's between ontology sets, the dest term will always be MONDO
    """
    EXACT = "exact", "exact"
    CLOSE = "close", "close"
    BROAD = "broad", "broad"
    NARROW = "narrow", "narrow"


class OntologyImport(TimeStampedModel):
    ontology_service = models.CharField(max_length=5, choices=OntologyService.choices)
    filename = models.TextField()
    context = models.TextField()
    hash = models.TextField()
    processed_date = models.DateTimeField(auto_created=True)
    completed = models.BooleanField(default=False)


class OntologyTerm(TimeStampedModel):

    """
    id is Term as it should be referenced <prefix>:<zero padded index> e.g.
    MONDO:0000043, OMIM:0000343
    """
    id = models.TextField(primary_key=True)
    ontology_service = models.CharField(max_length=5, choices=OntologyService.choices)
    index = models.IntegerField()
    name = models.TextField(null=True, blank=True)  # should only be null if we're using it as a placeholder reference
    definition = models.TextField(null=True, blank=True)
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)

    def __str__(self):
        return self.id

    class Meta:
        unique_together = ("ontology_service", "index")

    def get_absolute_url(self):
        return reverse('ontology_term', kwargs={"term": self.url_safe_id})

    @property
    def is_stub(self):
        return self._state.adding

    @property
    def url_safe_id(self):
        return self.id.replace(":", "_")

    @staticmethod
    def get_or_stub(id_str: str) -> 'OntologyTerm':
        parts = re.split("[:|_]", id_str)
        if len(parts) != 2:
            raise ValueError(f"Can not convert {id_str} to a proper id")

        prefix = OntologyService(parts[0])
        num_part = int(parts[1])
        clean_id = OntologyService.index_to_id(prefix, num_part)

        if existing := OntologyTerm.objects.filter(id=clean_id).first():
            return existing

        return OntologyTerm(
            id=clean_id,
            ontology_service=prefix,
            index=num_part,
            name="Not stored locally"
        )

    @property
    def padded_index(self) -> str:
        expected_length = OntologyService.EXPECTED_LENGTHS[self.ontology_service]
        return str(self.index).rjust(expected_length, '0')

    @property
    def url(self):
        return OntologyService.URLS[self.ontology_service].replace("${1}", self.padded_index)


class OntologyTermRelation(TimeStampedModel):
    """
    I haven't elected to use django_dag node_factory here as it only allows one kind of relationship
    and we have aliases, replacement for obsolete term and is_a - is_a being the only one that is a true
    relationship.
    Maybe the is_a should still be done via django_dag??
    """
    source_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE, related_name="subject")
    dest_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    relation = models.CharField(max_length=10, choices=OntologyRelation.choices)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)

    def __str__(self):
        return f"{self.source_term} -> ({self.relation}) -> {self.dest_term}"

    @property
    def relation_display(self):
        return self.relation.replace("_", " ")

    @staticmethod
    def parents_of(source_term: OntologyTerm) -> QuerySet:
        parent_ids = OntologyTermRelation.objects.filter(source_term=source_term, relation=OntologyRelation.IS_A).values_list("dest_term", flat=True)
        return OntologyTerm.objects.filter(id__in=parent_ids).order_by("-id")

    @staticmethod
    def children_of(term: OntologyTerm) -> QuerySet:
        """
        QuerySet of OntologyTerms
        """
        children_ids = OntologyTermRelation.objects.filter(dest_term=term, relation=OntologyRelation.IS_A).values_list("source_term", flat=True)
        return OntologyTerm.objects.filter(id__in=children_ids).order_by("-id")

    @staticmethod
    def all_relationships(term: OntologyTerm) -> List['OntologyTermRelation']:
        """
        Need to convert QuerySet to list to ensure ordering of outgoing vs incoming relationships
        """
        outgoing = OntologyTermRelation.objects.filter(source_term=term).order_by('-dest_term__id').all()
        incoming = OntologyTermRelation.objects.filter(dest_term=term).order_by('-source_term__id').all()
        return list(outgoing) + list(incoming)

    @staticmethod
    def mondo_version_of(term: OntologyTerm) -> Optional[OntologyTerm]:
        if term.ontology_service == OntologyService.MONDO:
            return term
        # Mondo will always be the source term
        if relationship := OntologyTermRelation.objects.filter(dest_term=term, source_term__ontology_service=OntologyService.MONDO, relation=OntologyRelation.EXACT).first():
            return relationship.source_term
        return None

    @staticmethod
    def omim_version_of(term: OntologyTerm) -> Optional[OntologyTerm]:
        if term.ontology_service == OntologyService.OMIM:
            return term
        # Mondo will always be the source term
        if relationship := OntologyTermRelation.objects.filter(dest_term__ontology_service=OntologyService.OMIM, source_term=term, relation=OntologyRelation.EXACT).first():
            return relationship.dest_term
        return None


class OntologyTermGeneRelation(TimeStampedModel):
    term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    relation = models.TextField()
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)
    deletion_date = models.DateTimeField(null=True, blank=True)

    @staticmethod
    def relationship_for(term: OntologyTerm) -> QuerySet:
        return OntologyTermGeneRelation.objects.filter(term=term, deletion_date__isnull=True)
