from django.db import models
from django.db.models import PROTECT, CASCADE
from model_utils.models import TimeStampedModel

from genes.models import GeneSymbol


class OntologySet(models.TextChoices):
    MONDO = "MONDO", "MONDO"
    OMIM = "OMIM", "OMIM"
    HPO = "HPO", "HPO"

    @staticmethod
    def index_to_id(ontology_set: 'OntologySet', index: int):
        EXPECTED_LENGTHS = {
            OntologySet.MONDO: 7,
            OntologySet.OMIM: 6,
            OntologySet.HPO: 7
        }
        expected_length = EXPECTED_LENGTHS.get(ontology_set, 0)
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
    ontology_set = models.CharField(max_length=5, choices=OntologySet.choices)
    filename = models.TextField()
    hash = models.TextField()
    notes = models.TextField()


class OntologyTerm(TimeStampedModel):
    """
    Term as it should be referenced <prefix>:<zero padded index> e.g.
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


class OntologyTermGeneRelation(TimeStampedModel):
    term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    relation = models.TextField()
    extra = models.JSONField(null=True, blank=True)
    from_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)