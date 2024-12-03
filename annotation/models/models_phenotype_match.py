from cache_memoize import cache_memoize
from django.contrib.auth.models import User
from django.db import models
from django.db.models import QuerySet, OuterRef, Count, Subquery
from django.db.models.deletion import CASCADE, SET_NULL

from library.constants import DAY_SECS
from ontology.models import OntologyTerm, OntologyVersion
from patients.models import Patient

PATIENT_TPM_PATH = "patient_text_phenotype__phenotype_description__textphenotypesentence__text_phenotype__textphenotypematch"
PATIENT_ONTOLOGY_TERM_PATH = PATIENT_TPM_PATH + "__ontology_term"


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

    @cache_memoize(timeout=DAY_SECS, args_rewrite=lambda s: (s.pk, ))
    def get_ontology_term_ids(self) -> list[int]:
        ot_qs = self.textphenotypesentence_set.filter(text_phenotype__textphenotypematch__ontology_term__isnull=False)
        # Sort so can be cached
        ot_qs = ot_qs.order_by("text_phenotype__textphenotypematch__ontology_term_id")
        return ot_qs.values_list("text_phenotype__textphenotypematch__ontology_term_id", flat=True)

    def get_gene_symbols(self, ontology_version) -> QuerySet:
        terms = tuple(self.get_ontology_term_ids())
        gene_symbols = ontology_version.cached_gene_symbols_for_terms_tuple(terms)
        return gene_symbols

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

    def get_ambiguous_matches(self):
        """ Where an Ontology Service has multiple matches to the exact same text """

        tpm_qs = self.textphenotypematch_set.all()
        outer = {
            "ontology_term__ontology_service": OuterRef('ontology_term__ontology_service'),
            "offset_start": OuterRef('offset_start'),
            "offset_end": OuterRef('offset_end'),
        }
        sub_query = tpm_qs.filter(**outer).values(*outer.keys()).annotate(num_terms=Count("pk")).filter(num_terms__gte=2)
        return tpm_qs.annotate(num_terms=Subquery(sub_query.values("num_terms"))).filter(num_terms__gte=2)



class TextPhenotypeSentence(models.Model):
    """ A description broken up into sentences which are represented by TextPhenotypes and matches """
    phenotype_description = models.ForeignKey(PhenotypeDescription, on_delete=CASCADE)
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    sentence_offset = models.IntegerField()

    def get_results(self):
        results = []

        ambiguous = set(self.text_phenotype.get_ambiguous_matches())
        for tpm in self.text_phenotype.textphenotypematch_set.all():
            tpm.offset_start += self.sentence_offset
            tpm.offset_end += self.sentence_offset
            data = tpm.to_dict()
            if tpm in ambiguous:
                data["ambiguous"] = tpm.match_text
            results.append(data)
        return results


class TextPhenotypeMatch(models.Model):
    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)
    offset_start = models.IntegerField()
    offset_end = models.IntegerField()

    @property
    def match_text(self) -> str:
        txt = self.text_phenotype.text
        return txt[self.offset_start:self.offset_end]

    def to_dict(self):
        """ This is what's sent as JSON back to client for highlighting and grids """

        ontology_version = OntologyVersion.latest()
        gene_symbols_qs = ontology_version.cached_gene_symbols_for_terms_tuple((self.ontology_term.pk,))
        gene_symbols = list(gene_symbols_qs.values_list("symbol", flat=True))
        accession = str(self.ontology_term)
        return {
            "accession": accession,
            "gene_symbols": gene_symbols,
            "match": accession,
            "ontology_service": self.ontology_term.get_ontology_service_display(),
            "name": self.ontology_term.name,
            "offset_start": self.offset_start,
            "offset_end": self.offset_end,
            "pk": self.ontology_term.pk,
        }

    def __str__(self):
        return f"{self.ontology_term} from ({self.offset_start}-{self.offset_end})"


class PatientTextPhenotype(models.Model):
    """ Used to link patient & phenotype """
    patient = models.OneToOneField(Patient, related_name='patient_text_phenotype', on_delete=CASCADE)
    phenotype_description = models.OneToOneField(PhenotypeDescription, on_delete=CASCADE)
    approved_by = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    def __str__(self):
        return f"{self.patient}: {self.phenotype_description}"


def patients_qs_for_ontology_term(user, ontology_term):
    return Patient.filter_for_user(user).filter(**{PATIENT_ONTOLOGY_TERM_PATH: ontology_term}).order_by("id")
