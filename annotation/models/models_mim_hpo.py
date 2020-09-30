from django.db import models
from django.db.models.deletion import CASCADE
from django_dag.models import edge_factory, node_factory

from annotation.models.models_enums import HPOSynonymScope


# MIM gives multiple names for disorders, separated by ";;". We store the 1st one as the name, and aliases in a 2nd table
class MIMMorbid(models.Model):
    PREFIX = "OMIM:"
    accession = models.IntegerField(primary_key=True)
    description = models.TextField()

    def get_genes(self):
        from genes.models import Gene
        return Gene.objects.filter(mimgene__mim_morbid=self).distinct()

    def __str__(self):
        return f"{self.accession} {self.description}"


# OMIM also contains aliases (ie multiple terms for the same underlying disease)
# Save these as MIMMorbidAlias (to get all the terms) we can auto-complete on.
# This also includes an entry from MIMMorbid (the default term)
class MIMMorbidAlias(models.Model):
    mim_morbid = models.ForeignKey(MIMMorbid, on_delete=CASCADE)
    description = models.TextField()

    @property
    def accession(self):
        return f"OMIM: {self.mim_morbid.pk}"

    @property
    def name(self):
        description = self.mim_morbid.description
        is_alias = self.mim_morbid.description != self.description
        if is_alias:
            description += f" (aka {self.description})"
        return description

    def get_genes(self):
        return self.mim_morbid.get_genes()

    def __str__(self):
        return f"{self.accession} {self.name}"


class HumanPhenotypeOntology(node_factory('HPOEdge')):
    """ Inserted with PK = human phenotype ontology ID """
    PREFIX = "HP:"
    name = models.TextField(null=True)
    definition = models.TextField(null=True)

    @property
    def accession(self):
        return self.PREFIX + "%07d" % int(self.pk)

    def get_genes(self):
        from genes.models import Gene
        return Gene.objects.filter(mimgene__mim_morbid__phenotypemim__hpo=self).distinct()

    def __str__(self):
        return f"{self.accession}: {self.name}"


class HPOEdge(edge_factory(HumanPhenotypeOntology, concrete=False)):
    pass


# This also contains HPO entries so we can auto-complete all HPO/Alias terms
# In that case name = hpo.name and you won't see AKA in repr()
class HPOSynonym(models.Model):
    hpo = models.ForeignKey(HumanPhenotypeOntology, on_delete=CASCADE)
    name = models.TextField()
    scope = models.CharField(max_length=1, choices=HPOSynonymScope.CHOICES)

    class Meta:
        unique_together = ('hpo', 'name')

    def __str__(self):
        description = str(self.hpo)
        if self.name != self.hpo.name:
            description += f" (aka {self.name})"
        return description


class PhenotypeMIM(models.Model):
    hpo = models.ForeignKey(HumanPhenotypeOntology, on_delete=CASCADE)
    mim_morbid = models.ForeignKey(MIMMorbid, on_delete=CASCADE)

    class Meta:
        unique_together = ('hpo', 'mim_morbid')

    def __str__(self):
        return f"{self.hpo}/{self.mim_morbid}"
