from django.db import models
from django.db.models import CASCADE
from django.urls import reverse
from django_dag.models import node_factory, edge_factory

from snpdb.models import SoftwareVersion


class VariantCaller(SoftwareVersion):
    run_params = models.TextField(null=True, blank=True)

    def __str__(self):
        description = f"{self.name} (v. {self.version})"
        if self.run_params:
            description += f": {self.run_params}"
        return description

    def get_absolute_url(self):
        return reverse('view_variant_caller', kwargs={'pk': self.pk})


class Aligner(SoftwareVersion):

    def __str__(self):
        return f"{self.name} (v. {self.version})"

    def get_absolute_url(self):
        return reverse('view_aligner', kwargs={'pk': self.pk})


class VariantCallingPipeline(SoftwareVersion):
    description = models.TextField(null=True, blank=True)
    aligner = models.ForeignKey(Aligner, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, on_delete=CASCADE)
    other_details = models.TextField(null=True, blank=True)

    def get_absolute_url(self):
        return reverse('view_variant_calling_pipeline', kwargs={'pk': self.pk})

    def __str__(self):
        return f"{self.description} (Algn: {self.aligner}, V. Caller: {self.variant_caller})"

# Software Pipeline Models replace VariantCallingPipeline and use the acyclic
# graph module from django_dag


class SoftwarePipeline(SoftwareVersion):
    description = models.TextField(null=True, blank=True)


class SoftwarePipelineNode(node_factory('SoftwarePipelineEdge')):
    softwarepipeline = models.ForeignKey(SoftwarePipeline, on_delete=CASCADE)
    name = models.TextField()
    version = models.TextField()
    parameters = models.TextField()
    description = models.TextField()


class SoftwarePipelineEdge(edge_factory(SoftwarePipelineNode, concrete=False)):
    name = models.CharField(max_length=32, blank=True, null=True)
