from typing import List

from django.contrib.auth.models import User
from django.db import models
from django.db.models import QuerySet
from django.db.models.deletion import CASCADE, SET_NULL

from genes.models import GeneSymbol, Gene
from ontology.models import OntologyTerm, OntologySnake
from patients.models import Patient


PATIENT_TPM_PATH = "patient_text_phenotype__phenotype_description__textphenotypesentence__text_phenotype__textphenotypematch"
PATIENT_ONTOLOGY_TERM_PATH = PATIENT_TPM_PATH + "__ontology_term"
PATIENT_GENE_SYMBOL_PATH = PATIENT_TPM_PATH + "__gene_symbol"


class DescriptionProcessingStatus(models.Model):
    CREATED = 'C'
    TOKENIZED = 'T'
    ERROR = 'E'
    SUCCESS = 'S'

    CHOICES = [
        (CREATED, 'Created'),
        (TOKENIZED, 'Tokenized'),
        (ERROR, 'Error'),
        (SUCCESS, 'Success'),
    ]


class PhenotypeDescription(models.Model):
    """ This represents a whole phenotype description - body of text which is broken up into multiple sentences
        If an object uses PhenotypeDescription - @see HasPhenotypeDescriptionMixin
    """

    original_text = models.TextField()
    status = models.CharField(max_length=1, choices=DescriptionProcessingStatus.CHOICES)

    def get_results(self):
        results = []
        for sentence in self.textphenotypesentence_set.all():
            results.extend(sentence.get_results())
        return results

    def get_ontology_term_ids(self) -> List[OntologyTerm]:
        ot_qs = self.textphenotypesentence_set.filter(text_phenotype__textphenotypematch__ontology_term__isnull=False)
        return ot_qs.values_list("text_phenotype__textphenotypematch__ontology_term_id", flat=True)

    def get_gene_symbols(self) -> QuerySet:
        ontology_term_ids = self.get_ontology_term_ids()
        return OntologySnake.special_case_gene_symbols_for_terms(ontology_term_ids)

    def __str__(self):
        text = self.original_text[:50]
        status = self.get_status_display()
        name = f"Text: {text} ({status})"
        name += '\n'.join(map(str, self.get_results()))
        return name


class TextPhenotype(models.Model):
    """ This is a sentence, that has matches hanging off it """

    text = models.TextField(primary_key=True)
    processed = models.BooleanField(default=False)

    def __str__(self):
        return f"{self.text} (processed: {self.processed})"


class TextPhenotypeSentence(models.Model):
    """ A description broken up into sentences which are represented by TextPhenotypes and matches """
    phenotype_description = models.ForeignKey(PhenotypeDescription, on_delete=CASCADE)
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    sentence_offset = models.IntegerField()

    def get_results(self):
        results = []

        for r in self.text_phenotype.textphenotypematch_set.all():
            r.offset_start += self.sentence_offset
            r.offset_end += self.sentence_offset
            results.append(r.to_dict())

        return results


class PhenotypeMatchTypes(models.Model):
    HPO = 'H'
    OMIM = 'O'
    GENE = 'G'
    CHOICES = [
        (HPO, 'HPO'),
        (OMIM, 'OMIM'),
        (GENE, 'Gene'),
    ]


class TextPhenotypeMatch(models.Model):
    FIELDS = {PhenotypeMatchTypes.HPO: 'ontology_term',
              PhenotypeMatchTypes.OMIM: 'ontology_term',
              PhenotypeMatchTypes.GENE: 'gene_symbol'}

    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    match_type = models.CharField(max_length=1, choices=PhenotypeMatchTypes.CHOICES)
    ontology_term = models.ForeignKey(OntologyTerm, null=True, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    offset_start = models.IntegerField()
    offset_end = models.IntegerField()

    def to_dict(self):
        """ This is what's sent as JSON back to client for highlighting and grids """
        if self.match_type == PhenotypeMatchTypes.GENE:
            gene_symbols = [self.gene_symbol_id]
        else:
            gene_symbols_qs = OntologySnake.special_case_gene_symbols_for_terms([self.ontology_term_id])
            gene_symbols = list(gene_symbols_qs.values_list("symbol", flat=True))

        string = str(self.record)

        return {"accession": string,
                "gene_symbols": gene_symbols,
                "match": string,
                "match_type": self.match_type,
                "name": self.record.name,
                "offset_start": self.offset_start,
                "offset_end": self.offset_end,
                "pk": self.record.pk}

    @property
    def record(self):
        f = TextPhenotypeMatch.FIELDS[self.match_type]
        return getattr(self, f)

    def __str__(self):
        return f"{self.record} from ({self.offset_start}-{self.offset_end})"


class PatientTextPhenotype(models.Model):
    """ Used to link patient & phenotype """
    patient = models.OneToOneField(Patient, related_name='patient_text_phenotype', on_delete=CASCADE)
    phenotype_description = models.OneToOneField(PhenotypeDescription, on_delete=CASCADE)
    approved_by = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    def __str__(self):
        return f"{self.patient}: {self.phenotype_description}"
