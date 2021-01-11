from typing import Optional

from django.db import models
from django.db.models import PROTECT, CASCADE
from model_utils.models import TimeStampedModel

from genes.models import GeneSymbol
from library.utils import Constant


class OntologySet(models.TextChoices):
    MONDO = "MONDO"
    OMIM = "OMIM"
    HPO = "HPO"

    EXPECTED_LENGTHS = Constant({
        MONDO: 7,
        OMIM: 6,
        HPO: 7
    })

    URLS = Constant({
        MONDO: "https://vm-monitor.monarchinitiative.org/disease/MONDO:${1}",
        OMIM: "http://www.omim.org/entry/${1}",
        HPO: "https://hpo.jax.org/app/browse/term/HP:${1}"
    })

    @staticmethod
    def index_to_id(ontology_set: 'OntologySet', index: int):
        expected_length = OntologySet.EXPECTED_LENGTHS.get(ontology_set, 0)
        num_part = str(index).rjust(expected_length, '0')
        return f"{ontology_set}:{num_part}"

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
    #TODO maybe work out better name than ontology_set
    ontology_set = models.CharField(max_length=5, choices=OntologySet.choices)
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
    ontology_set = models.CharField(max_length=5, choices=OntologySet.choices)
    index = models.IntegerField()
    name = models.TextField(null=True, blank=True) # should only be null if we're using it as a placeholder reference
    definition = models.TextField(null=True, blank=True)
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)

    class Meta:
        unique_together = ("ontology_set", "index")

    @property
    def is_stub(self):
        return self._state.adding

    @staticmethod
    def get_or_stub(id_str: str) -> 'OntologyTerm':
        parts = id_str.split(":")
        if len(parts) != 2:
            raise ValueError(f"Can not convert {id_str} to a proper id")

        prefix = OntologySet(parts[0])
        num_part = int(parts[1])
        clean_id = OntologySet.index_to_id(prefix, num_part)

        if existing := OntologyTerm.objects.filter(id=clean_id).first():
            return existing

        return OntologyTerm(
            id=clean_id,
            ontology_set=prefix,
            index=num_part,
            name="Not stored locally"
        )

    @property
    def padded_index(self) -> str:
        expected_length = OntologySet.EXPECTED_LENGTHS[self.ontology_set]
        return str(self.index).rjust(expected_length, '0')

    @property
    def url(self):
        return OntologySet.URLS[self.ontology_set].replace("${1}", self.padded_index)


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

    @staticmethod
    def parent_of(source_term: OntologyTerm) -> Optional[OntologyTerm]:
        if relationship := OntologyTermRelation.objects.filter(source_term=source_term, relation=OntologyRelation.IS_A).first():
            return relationship.dest_term
        return None

    @staticmethod
    def mondo_version_of(term: OntologyTerm) -> Optional[OntologyTerm]:
        if term.ontology_set == OntologySet.MONDO:
            return term
        # Mondo will always be the source term
        if relationship := OntologyTermRelation.objects.filter(dest_term=term, source_term__ontology_set=OntologySet.MONDO, relation=OntologyRelation.EXACT).first():
            return relationship.dest_term
        return None

    @staticmethod
    def omim_version_of(term: OntologyTerm) -> Optional[OntologyTerm]:
        if term.ontology_set == OntologySet.OMIM:
            return term
        # Mondo will always be the source term
        if relationship := OntologyTermRelation.objects.filter(dest_term__ontology_set=OntologySet.OMIM, source_term=term, relation=OntologyRelation.EXACT).first():
            return relationship.dest_term
        return None


class OntologyTermGeneRelation(TimeStampedModel):
    term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    relation = models.TextField()
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)