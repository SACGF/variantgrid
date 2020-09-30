from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL

from annotation.models.models_mim_hpo import HumanPhenotypeOntology, MIMMorbidAlias, MIMMorbid
from genes.models import GeneSymbol, Gene
from patients.models import Patient


PATIENT_TPM_PATH = "patient_text_phenotype__phenotype_description__textphenotypesentence__text_phenotype__textphenotypematch"
PATIENT_OMIM_PATH = PATIENT_TPM_PATH + "__omim_alias__mim_morbid"
PATIENT_HPO_PATH = PATIENT_TPM_PATH + "__hpo"
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

    def get_hpo_qs(self):
        hpo_path = 'text_phenotype__textphenotypematch__hpo'
        pheno_qs = self.textphenotypesentence_set.filter(**{hpo_path + "__isnull": False})
        pheno_ids = pheno_qs.values_list(hpo_path, flat=True)
        return HumanPhenotypeOntology.objects.filter(pk__in=pheno_ids)  # @UndefinedVariable

    def get_mim_qs(self):
        mim_morbid_alias_qs = self.textphenotypesentence_set.filter(text_phenotype__textphenotypematch__omim_alias__isnull=False)
        omim_ids = mim_morbid_alias_qs.values_list("text_phenotype__textphenotypematch__omim_alias__mim_morbid_id", flat=True)
        return MIMMorbid.objects.filter(pk__in=omim_ids)

    def get_mim_and_pheno_mim_qs(self):
        tps = self.textphenotypesentence_set
        mim_morbid_qs = tps.filter(text_phenotype__textphenotypematch__omim_alias__isnull=False)
        mim_morbid_values = mim_morbid_qs.values_list("text_phenotype__textphenotypematch__omim_alias__mim_morbid_id", flat=True)
        pheno_morbid_qs = tps.filter(text_phenotype__textphenotypematch__hpo__phenotypemim__mim_morbid__isnull=False)
        pheno_morbid_values = pheno_morbid_qs.values_list("text_phenotype__textphenotypematch__hpo__phenotypemim__mim_morbid", flat=True)
        omim_ids = set(pheno_morbid_values) | set(mim_morbid_values)
        return MIMMorbid.objects.filter(pk__in=omim_ids)

    def get_gene_qs(self):
        mimmorbid_qs = self.get_mim_and_pheno_mim_qs()
        return Gene.objects.filter(mimgene__mim_morbid__in=mimmorbid_qs)

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
    FIELDS = {PhenotypeMatchTypes.HPO: 'hpo',
              PhenotypeMatchTypes.OMIM: 'omim_alias',
              PhenotypeMatchTypes.GENE: 'gene_symbol'}

    text_phenotype = models.ForeignKey(TextPhenotype, on_delete=CASCADE)
    match_type = models.CharField(max_length=1, choices=PhenotypeMatchTypes.CHOICES)
    hpo = models.ForeignKey(HumanPhenotypeOntology, null=True, on_delete=CASCADE)
    omim_alias = models.ForeignKey(MIMMorbidAlias, null=True, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    offset_start = models.IntegerField()
    offset_end = models.IntegerField()

    def to_dict(self):
        """ This is what's sent as JSON back to client for highlighting and grids """
        gene_symbols_qs = GeneSymbol.objects.filter(geneversion__gene__in=self.record.get_genes()).order_by("pk")
        gene_symbols = list(gene_symbols_qs.values_list("symbol", flat=True))

        return {"accession": self.record.accession,
                "gene_symbols": gene_symbols,
                "match": str(self.record),
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
